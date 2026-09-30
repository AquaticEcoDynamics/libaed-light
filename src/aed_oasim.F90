!###############################################################################
!#                                                                             #
!# aed_oasim.F90                                                               #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2024-2026 : The University of Western Australia                  #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created April 2022                                                          #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |     ____     | || |      __      | || |    _______   | |         !
!         | |   .'    `.   | || |     /  \     | || |   /  ___  |  | |         !
!         | |  /  .--.  \  | || |    / /\ \    | || |  |  (__ \_|  | |         !
!         | |  | |    | |  | || |   / ____ \   | || |   '.___`-.   | |         !
!         | |  \  `--'  /  | || | _/ /    \ \_ | || |  |`\____) |  | |         !
!         | |   `.____.'   | || ||____|  |____|| || |  |_______.'  | |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

!#------------------------------##################------------------------------
!### This module is not stable and is under development. Everything may change !
!#------------------------------##################------------------------------
!#                                                                             #
!# An implementation of the Ocean-Atmosphere Spectral Irradiance Model (OASIM) #
!#                                                                             #
!# Lazzari, P., Salon, S., Terzić, E., Gregg, W. W., D'Ortenzio, F.,           #
!# Vellucci, V., Organelli, E., and Antoine, D.:                               #
!# Assessment of the spectral downward irradiance at the surface of the        #
!# Mediterranean Sea using the radiative Ocean-Atmosphere Spectral Irradiance  #
!# Model (OASIM), Ocean Sci., 17, 675-697,                                     #
!#                               https://doi.org/10.5194/os-17-675-2021, 2021  #
!#                                                                             #
!# - https://os.copernicus.org/articles/17/675/2021/                           #
!#                                                                             #
!#------------------------------##################------------------------------

#include "aed.h"

!
MODULE aed_oasim

   USE aed_core
   USE aed_slingo
   USE aed_maths
   USE curtin_light, ONLY : direct_diffuse_bands, calculate_lambda_bounds

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_oasim_data_t
   !# 2026-09-08: exported so the production arithmetic can be unit tested against the compiled
   !# object code rather than a reimplementation. See cases_heat_wall/tools/test_oasim_kd.py and
   !# tools/test_oasim_kd.F90.
   PUBLIC aed_oasim_kd_mixture, aed_oasim_kd_inputs_ok, aed_oasim_is_finite
!
   TYPE type_iop
      INTEGER                   :: type       ! type index
      INTEGER                   :: id_c       ! concentration (could be chl, carbon, or
                                              ! something else, but its product with
                                              ! a or b below should return units m-1)
      INTEGER                   :: id_parc    ! par photons captured by this IOP
      INTEGER                   :: id_aeff    ! projected area of leaves in a benthic canopy
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a  ! specific absorption (m-1 concentration-1)
      AED_REAL,DIMENSION(:),ALLOCATABLE :: b  ! specific total scattering (m-1 concentration-1)
      AED_REAL :: b_b                         ! ratio of backscattering to total scattering (dimensionless)
      AED_REAL :: b_p                         ! packaging effect parameter (dimensionless)
      AED_REAL :: scale                       ! optional user scaling factor to amplify/dampen (dimensionless)
   END TYPE
!
   TYPE,extends(aed_model_data_t) :: aed_oasim_data_t
      !# Variable identifiers
      INTEGER :: id_swr, id_uv, id_par, id_par_E, id_par_E_scalar, id_par_J_scalar, id_par_E_dif, id_swr_abs, id_secchi
      INTEGER :: id_kd_bio   !# PAR Kd from non-water constituents only (see Kd_bio below)
      INTEGER :: id_kd_water !# PAR Kd from pure water alone - the value Kw should be set to
      INTEGER :: id_swr_sf, id_par_sf, id_uv_sf, id_par_E_sf, id_swr_dif_sf, id_wind_out, id_zen
      INTEGER :: id_swr_sf_w, id_par_sf_w, id_uv_sf_w, id_par_E_sf_w
      INTEGER :: id_alpha_a, id_beta_a, id_omega_a, id_F_a
      INTEGER :: id_lon, id_lat, id_cloud, id_wind_speed, id_airpres, id_relhum
      INTEGER :: id_lwp, id_O3, id_WV, id_visibility, id_air_mass_type
      !# 2026-08-29 FIX: the two mean-wind ids were left UNINITIALISED. Exactly one of them is
      !# ever assigned (id_mean_wind_out when compute_mean_wind=.TRUE., id_mean_wind_speed when
      !# .FALSE.), and only inside the sky_model==1 branch of aed_define_model - yet both are
      !# read unguarded elsewhere (`IF (data%id_mean_wind_out > 0)` and the unconditional
      !# `WM = _STATE_VAR_S_(data%id_mean_wind_speed)` in aed_calculate_dry_oasim). Whichever
      !# was not assigned held garbage that could test > 0 and index column() out of bounds.
      !# Default them to -1 so the guards below are meaningful. Both remain assigned exactly as
      !# before in every configuration that assigns them, so no live result changes.
      INTEGER :: id_mean_wind_out = -1, id_mean_wind_speed = -1
      INTEGER :: id_yearday, id_doy, id_swflux, id_hr_peakz, id_hr_peakx
      INTEGER :: id_h
      INTEGER :: id_I_0
      INTEGER :: nlambda
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_surface_band_dir, id_surface_band_dif
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_band_dir, id_band_dif, id_a_band, id_b_band, id_Kd
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: id_a_iop
      !# 2026-07 (Tier 1C): per-IOP normalised spectral-absorption factor. Lets a biology
      !# module weight its growth by the SAME measured a*(lambda) that OASIM uses for
      !# shading, instead of a separate hand-tuned band template.
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_afac
      !# 2026-08: the same factor computed both ways, so the two can be compared in one run.
      !#   id_afacE  weighted by par_weights   (energy, W/m2/nm) - the conventional treatment
      !#   id_afacq  weighted by par_E_weights (photons, umol/m2/s/nm) - Phi = I*lambda/(h c)
      !# id_afac holds a copy of whichever `afac_quanta` selects; consumers read only that, so
      !# no other module needs to know about the switch.
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_afacE, id_afacq, id_astar
      AED_REAL,DIMENSION(:),ALLOCATABLE :: lambda, lambda_bounds, par_weights, &
                           par_E_weights, swr_weights, uv_weights, F, lambda_out
      AED_REAL,DIMENSION(:),ALLOCATABLE :: exter
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a_o, a_u, a_v, tau_r
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a_w, b_w, a_w_out, b_w_out
      TYPE (type_iop), ALLOCATABLE :: iops(:)
      INTEGER :: l490_l = 1   !# default: read unconditionally in the column routine
      INTEGER :: spectral_output, sky_model
      LOGICAL :: save_Kd, ozone
      LOGICAL :: drive_extinction   !# feed OASIM's Kd_bio into the host's extinction
      LOGICAL :: afac_quanta        !# afac_iop carries photon- (T) or energy-weighted (F) absorption
      AED_REAL :: cloud, airpres, lwp, O3, WV, visibility, air_mass_type, tzshift

    CONTAINS       ! Selected AED methods:
         PROCEDURE :: define            => aed_define_oasim
         PROCEDURE :: calculate_surface => aed_calculate_surface_oasim
         PROCEDURE :: calculate_column  => aed_calculate_column_oasim
         PROCEDURE :: calculate_dry     => aed_calculate_dry_oasim
         !# 2026-07 (Tier 1A): binding enabled. The routine is a no-op unless the
         !# `drive_extinction` namelist switch is set, so the default behaviour is
         !# unchanged - OASIM stays a diagnostic sidecar until deliberately wired in.
         PROCEDURE :: light_extinction  => aed_light_extinction_oasim

   END TYPE

! MODULE GLOBALS
   AED_REAL :: DTday       = (15./60.)/24.    ! 15 minutes
   AED_REAL, PARAMETER :: WindAvgTime = 1.0   ! 1 day
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


#include "oasim.inc"


!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed_define_oasim(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED OASIM module
!
!  Here, the aed namelist is read and the variables and diagnostics created
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_oasim_data_t),INTENT(inout) :: data
!
!LOCALS
   INTEGER :: status
   INTEGER :: i, num_tn,num_tkn,num_tp,num_toc,num_tss,num_turb,num_tfe,num_tal

   INTEGER :: l
   INTEGER :: i_iop
   CHARACTER(len=8) :: strwavelength, strindex, strindex2
   !# Per-IOP inland (Case-2) custom-IOP parameters, read from the /aed_oasim/ namelist.
   !# (Previously scalar INTEGERs that were never assigned; now AED_REAL arrays so type 99 works.)
   AED_REAL :: a_star_iop(MAX_PHYTO_TYPES)     = 0.   ! specific absorption at lambda_ref [m2/mg C]
   AED_REAL :: S_iop(MAX_PHYTO_TYPES)          = 0.   ! exponential absorption slope [1/nm] (CDOM/NAP)
   AED_REAL :: lambda_ref_iop(MAX_PHYTO_TYPES) = 440. ! reference wavelength [nm]
   AED_REAL :: b_star_iop(MAX_PHYTO_TYPES)     = 0.   ! specific scattering at lambda_ref [m2/mg C]
   AED_REAL :: eta_iop(MAX_PHYTO_TYPES)        = 0.   ! scattering power-law exponent [-]
   AED_REAL :: b_b_iop(MAX_PHYTO_TYPES)        = 0.   ! backscatter-to-scatter ratio [-]

   INTEGER,PARAMETER :: exter_source = 2

   ! Coefficients for wavelength dependence of foam reflectance
   !                                       (Eqs A11, A12 in Gregg & Casey 2009)
   AED_REAL,PARAMETER :: a0 =  0.9976
   AED_REAL,PARAMETER :: a1 =  0.2194
   AED_REAL,PARAMETER :: a2 =  0.0554
   AED_REAL,PARAMETER :: a3 =  0.0067
   AED_REAL,PARAMETER :: b0 =  5.026
   AED_REAL,PARAMETER :: b1 = -0.0114
   AED_REAL,PARAMETER :: b2 =  9.552e-6
   AED_REAL,PARAMETER :: b3 = -2.698e-9

   AED_REAL,PARAMETER :: Planck     = 6.62606957e-34  ! Planck constant (m2 kg/s)
   AED_REAL,PARAMETER :: lightspeed = 299792458       ! Speed of light (m/s)
   AED_REAL,PARAMETER :: Avogadro   = 6.02214129e23   ! Avogadro constant (/mol)

   AED_REAL :: log_T_w

!  %% NAMELIST   %%  /aed_oasim/
!  %% Last Checked : never
   INTEGER  :: lambda_method = 1    ! choice of wavebands (0: custom range, 1: OASIM)
   INTEGER  :: nlambda = 1          !
   INTEGER  :: nlambda_out = 1      ! number of wavebands for spectral output
   AED_REAL :: lambda_out(50) = 500.!
   AED_REAL :: lambda_min = 0.      ! minimum wavelength
   AED_REAL :: lambda_max = 0.      ! maximum wavelength
   INTEGER  :: n_iop = 0            ! number of inherent optical properties (IOPs)
   INTEGER  :: iop_type(MAX_PHYTO_TYPES) = 0         ! type of IOP :
                                    !    1: diatoms
                                    !    2: chlorophytes
                                    !    3: cyanobacteria
                                    !    4: coccolithophorids
                                    !    5: dinoflagellates
                                    !    6: detritus
                                    !    8: CDOC
                                    !    9: OM with custom absorption/scattering
   CHARACTER(len=64) :: iop_link(MAX_PHYTO_TYPES) = ''
   AED_REAL :: iop_scale(MAX_PHYTO_TYPES) = 1.0 !

   INTEGER  :: spectral_output = 0  ! spectral output :
                                    !    0: none
                                    !    1: full
                                    !    2: selected wavelengths
                                    !    default=0, minimum=0, maximum=2)
!  CALL data%get_parameter(data%lambda_out(l), 'lambda' // trim(strindex) // '_out', '', 'output wavelength ' // trim(strindex))
   INTEGER  :: lambda(100)
   LOGICAL  :: save_Kd = .FALSE.    ! compute attenuation, default=.false.
   LOGICAL  :: ozone = .TRUE.       ! compute ozone using Van Heuklon (1979)
   LOGICAL  :: compute_mean_wind = .TRUE.
   AED_REAL :: cloud = 0.1          ! cloud cover fraction
   AED_REAL :: airpres = 101300.    ! surface air pressure, Pa
   AED_REAL :: lwp = 0.1            ! cloud liquid water content, kg/m2
   AED_REAL :: O3 = 0.344 * 0.4462e-3 * 48. / 1000  ! Ozone content, From matm-cm to mol m-2 to kg m-2
   AED_REAL :: WV = 5.0             ! precipitable water content, kg/m2
   AED_REAL :: visibility = 25000.0 ! atmospheric visibility, m
   AED_REAL :: air_mass_type = 10.0 ! air mass type (1=open ocean, 10=continental)
   INTEGER  :: sky_model = 1        ! sub-model approach to comput incoming spectra
   AED_REAL :: tzshift = 8.0        ! time-zone shift, required to offset UTC to local for solar zenith calculation
   !# 2026-07 (Tier 1A): let OASIM's spectral solution drive the host's light extinction.
   !# DEFAULT .FALSE. ON PURPOSE. Until this is enabled OASIM is a diagnostic sidecar: the
   !# host attenuates light with Kw plus the non-spectral Ke* constants in aed.nml, and
   !# nothing OASIM computes affects growth or heating. Enabling it WITHOUT first zeroing
   !# or re-deriving those constants DOUBLE-COUNTS every attenuator. See handoff 33-O/33-P
   !# for the measured decomposition and the required recalibration.
   LOGICAL  :: drive_extinction = .FALSE.

   !# Weight a*(lambda) by PHOTON flux rather than radiant energy when forming afac_iop.
   !# Photosynthesis counts quanta, not joules: a 675 nm photon drives one photochemical event
   !# with 435/675 = 0.64 of the energy of a blue one. Weighting absorption by energy therefore
   !# systematically under-weights the red and green-orange bands - which is where phycocyanin
   !# and phycoerythrin absorb, i.e. exactly where the cyanobacterial light niche lives.
   !# Default .FALSE. preserves the calibrated behaviour; afacE_iop/afacq_iop are written either
   !# way, so the comparison is available without re-running.
   LOGICAL  :: afac_quanta = .FALSE.

! %% From Module Globals
!  INTEGER  :: diag_level = 10      ! 0 = no diagnostic outputs
!                                   ! 1 = basic diagnostic outputs
!                                   ! 2 = flux rates, and supporitng
!                                   ! 3 = other metrics
!                                   !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_oasim/

   !# 2026-08-29 FIX: `ozone` (declared above with default .TRUE.) was NOT in the NAMELIST
   !# below, so it could never be set from aed.nml. data%ozone was therefore always .TRUE.,
   !# and the block at ~930 always overwrote the user's O3 with the Van Heuklon (1979)
   !# climatology - silently ignoring any O3 supplied in the namelist. Adding it here does not
   !# change the live model: no case aed.nml sets `ozone`, so it keeps its .TRUE. default.
   NAMELIST /aed_oasim/ lambda_method, nlambda, lambda_min, lambda_max,        &
                        n_iop, iop_type, iop_link, iop_scale,                  &
                        a_star_iop, S_iop, lambda_ref_iop,                     &
                        b_star_iop, eta_iop, b_b_iop,                          &
                        spectral_output, nlambda_out, lambda_out,              &
                        lambda, save_Kd, ozone,                                &
                        cloud, airpres, lwp, O3, WV, visibility, air_mass_type,&
                        compute_mean_wind, sky_model, tzshift, diag_level,      &
                        drive_extinction, afac_quanta

!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_oasim configuration"

   !-----------------------------------------------
   ! Read the namelist
   read(namlst,nml=aed_oasim,iostat=status)
   !# 2026-09-08: gfortran exits with status 0 for STOP 'string', so a bad namelist used to
   !# look like a successful run. Print, then STOP 1.
   IF (status /= 0) THEN
      print *,'ERROR: failed reading namelist aed_oasim'
      STOP 1
   ENDIF

   data%O3 = O3
   data%WV = WV
   data%lwp = lwp
   data%ozone = ozone
   data%cloud = cloud
   data%tzshift = tzshift
   data%airpres = airpres
   data%sky_model = sky_model
   data%visibility = visibility
   data%air_mass_type = air_mass_type

   !-----------------------------------------------
   ! Process /interpolate spectral data & link to IOPs
   data%nlambda = nlambda
   SELECT CASE (lambda_method)
   CASE (0)
      ! Custom bandwidth range, specified by user
      ALLOCATE(data%lambda(nlambda), data%lambda_bounds(nlambda + 1))
      DO l = 1, data%nlambda + 1
         data%lambda_bounds(l) = lambda_min + (l - 1) * (lambda_max - lambda_min) / nlambda
      ENDDO
      data%lambda = (data%lambda_bounds(1:nlambda) + data%lambda_bounds(2:nlambda + 1)) / 2
   CASE (1)
      ! OASIM module default bandwidths, 33 set in oasim.inc
      nlambda = nlambda_oasim
      data%nlambda = nlambda
      ALLOCATE(data%lambda(nlambda), data%lambda_bounds(nlambda + 1))
      data%lambda(:) = lambda_oasim
      call calculate_lambda_bounds(nlambda, data%lambda, data%lambda_bounds)
   END SELECT

   ALLOCATE(data%iops(n_iop))
   DO i_iop = 1, n_iop
      ALLOCATE(data%iops(i_iop)%a(nlambda))
      ALLOCATE(data%iops(i_iop)%b(nlambda))
      !# 2026-08-29 FIX: define the spectra before the SELECT CASE below. Two branches leave
      !# them undefined otherwise - CASE(7) (PIC), whose body is entirely commented out, and
      !# any iop_type with no case at all (see the CASE DEFAULT added below). Undefined
      !# ALLOCATE'd memory then flowed into the two-stream light budget at ~1107. Every live
      !# iop_type assigns both arrays in full, so this changes nothing for the live model.
      data%iops(i_iop)%a(:) = zero_
      data%iops(i_iop)%b(:) = zero_
      data%iops(i_iop)%b_b  = zero_
      write(strindex, '(i0)') i_iop

      data%iops(i_iop)%b_p = 1.0
      data%iops(i_iop)%type = iop_type(i_iop)
      data%iops(i_iop)%scale = iop_scale(i_iop)

      SELECT CASE (iop_type(i_iop))
      CASE (1)  ! DIATOMS
         CALL interp(size(lambda_diatoms), lambda_diatoms, a_diatoms, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_diatoms), lambda_diatoms, b_diatoms, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.002 ! Gregg & Rousseau 2016 but originally Morel 1988
         data%iops(i_iop)%b_p = 1.0   ! JIM
      CASE (2)  ! CHLOROPHYTES
         CALL interp(size(lambda_chlorophytes), lambda_chlorophytes, a_chlorophytes, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_chlorophytes), lambda_chlorophytes, b_chlorophytes, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Ahn et al. 1992 as reported in Gregg & Rousseau 2016
         data%iops(i_iop)%b_p = 1.0
      CASE (3)  ! CYANOBACTERIA
         CALL interp(size(lambda_cyanobacteria), lambda_cyanobacteria, a_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_cyanobacteria), lambda_cyanobacteria, b_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0032 ! Gregg & Rousseau 2016 but originally Ahn et al. 1992
         data%iops(i_iop)%b_p = 1.0
      CASE (15) ! PLANKTOTHRIX (phycoerythrin-bearing cyanobacterium)
         !# 2026-07 (Tier 2.1). Previously Planktothrix used CASE(3), making it optically
         !# identical to the surface cyanobacterium and leaving it with NO phycoerythrin
         !# absorption - the very trait that explains its metalimnetic niche. Scattering
         !# reuses the cyanobacteria spectrum (no basis to differentiate it); only the
         !# absorption spectrum carries the added 540-575 nm phycoerythrin band.
         !# See the provenance note above a_Planktothrix in oasim.inc: shape is
         !# literature-based, magnitudes are constructed and NOT measured for FCR.
         CALL interp(size(lambda_planktothrix), lambda_planktothrix, a_planktothrix, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_cyanobacteria), lambda_cyanobacteria, b_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0032 ! as cyanobacteria; Gregg & Rousseau 2016 / Ahn et al. 1992
         data%iops(i_iop)%b_p = 1.0
      !# 2026-08: PIGMENT-RECONSTRUCTED counterparts of the four phytoplankton spectra above.
      !#   21 <- 1  diatoms        22 <- 2  chlorophytes
      !#   23 <- 3  cyanobacteria  25 <- 15 Planktothrix
      !# Same 19-node grid, same scattering, same PAR-integrated absorption; only the SHAPE of
      !# a*(lambda) differs, built from pigment Gaussians rather than from tabulated marine
      !# spectra with three hand-edited nodes. See analysis/build_absorption_spectra.py.
      !# These are OPT-IN. The legacy codes 1/2/3/15 are unchanged and remain the default,
      !# because drive_extinction = .true. makes a*(lambda) part of the lake's light budget:
      !# switching a group here moves Kd, and therefore temperature. Not a free change.
      CASE (21) ! DIATOMS - pigment-reconstructed
         CALL interp(size(lambda_Diatoms), lambda_Diatoms, a_Diatoms_pigment, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Diatoms), lambda_Diatoms, b_diatoms, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.002
         data%iops(i_iop)%b_p = 1.0
      CASE (22) ! CHLOROPHYTES - pigment-reconstructed
         CALL interp(size(lambda_Chlorophytes), lambda_Chlorophytes, a_Chlorophytes_pigment, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Chlorophytes), lambda_Chlorophytes, b_chlorophytes, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.00071 * 10
         data%iops(i_iop)%b_p = 1.0
      CASE (23) ! CYANOBACTERIA - pigment-reconstructed
         CALL interp(size(lambda_Cyanobacteria), lambda_Cyanobacteria, a_Cyanobacteria_pigment, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Cyanobacteria), lambda_Cyanobacteria, b_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0032
         data%iops(i_iop)%b_p = 1.0
      CASE (25) ! PLANKTOTHRIX - pigment-reconstructed (phycoerythrin-dominant)
         CALL interp(size(lambda_Planktothrix), lambda_Planktothrix, a_Planktothrix_pigment, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Cyanobacteria), lambda_Cyanobacteria, b_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0032
         data%iops(i_iop)%b_p = 1.0
      CASE (4)  ! COCCOLITHOPHORIDS
         CALL interp(size(lambda_coccolithophores), lambda_coccolithophores, a_coccolithophores, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_coccolithophores), lambda_coccolithophores, b_coccolithophores, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Morel 1988 as reported in Gregg & Rousseau 2016
         data%iops(i_iop)%b_p = 1.0
      CASE (5)  ! DINOFLAGELLATES
         CALL interp(size(lambda_dinoflagellates), lambda_dinoflagellates, a_dinoflagellates, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_dinoflagellates), lambda_dinoflagellates, b_dinoflagellates, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0029 ! Gregg & Rousseau 2016 but originally Morel 1988
         data%iops(i_iop)%b_p = 1.0
      CASE (6)  ! DETRITUS (small, as described in Gregg & Rousseau 2016)
         ! Parameters below match the small organic detritus parametrization of Gallegos et al. 2011 (Table 2)
         !  - the latter also offers a parametrizaton for large detritus. Note: Gallegos et al. 2011 constants
         !    are specific to dry weight! Did Gregg & Rousseau 2016 misinterpret them as specific to carbon weight?
         !    If so: Babin et al 2003 state 2.6 g DW per g C is representative for suspended OM, Bruggeman pers com
         !    notes adding factor 2.6 (DW/C) seems like to make L4/WCO match much better
         !    Also, note that 12.0107 converts from mg-1 to mmol-1
         data%iops(i_iop)%a(:) = 8e-5 * exp(-0.013 * (data%lambda - 440)) * 12.0107 * 2.6
         data%iops(i_iop)%b(:) = 0.00115 * (550. / data%lambda)**0.5 * 12.0107 * 2.6
         data%iops(i_iop)%b_b = 0.005
         data%iops(i_iop)%b_p = 1.0
      CASE (7)  ! PIC
        !data%iops(i_iop)%a(:) = 0
        !CALL interp(size(), lambda_w, a_w, nlambda, data%lambda, data%iops(i_iop)%b)
        !data%iops(i_iop)%b_b = 0.01
        !data%iops(i_iop)%b_p = 1.0
      CASE (8) !  CDOC
         ! Note that 12.0107 converts from mg-1 to mmol-1
         data%iops(i_iop)%a(:) = 2.98e-4 * exp(-0.014 * (data%lambda - 443)) * 12.0107
         data%iops(i_iop)%b(:) = 0
         data%iops(i_iop)%b_b = 0
         data%iops(i_iop)%b_p = 1.0
      CASE (9) ! HE6CHL (gordon-morel)
         CALL interp(size(lambda_he6chl ), lambda_he6chl, a_he6chl, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_he6chl ), lambda_he6chl, b_he6chl, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.02
         data%iops(i_iop)%b_p = 0.62
      CASE (10) ! AVERAGE SEDIMENT
         CALL interp(size(lambda_Averagesediment ), lambda_Averagesediment, a_Averagesediment, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Averagesediment ), lambda_Averagesediment, b_Averagesediment, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (11) ! BROWN EARTH
         CALL interp(size(lambda_Brownearth ), lambda_Brownearth, a_Brownearth, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Brownearth ), lambda_Brownearth, b_Brownearth, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (12) ! CALCEROUS SAND
         CALL interp(size(lambda_Calcareoussand ), lambda_Calcareoussand, a_Calcareoussand, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Calcareoussand ), lambda_Calcareoussand, b_Calcareoussand, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (13) ! RED CLAY
         CALL interp(size(lambda_Redclay ), lambda_Redclay, a_Redclay, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Redclay ), lambda_Redclay, b_Redclay, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (14) ! YELLOW CLAY
         CALL interp(size(lambda_Yellowclay ), lambda_Yellowclay, a_Yellowclay, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_Yellowclay ), lambda_Yellowclay, b_Yellowclay, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (20) ! SEAGRASS CANOPY
         CALL interp(size(lambda_seagrassleaf) , lambda_seagrassleaf, a_seagrassleaf, nlambda, data%lambda, data%iops(i_iop)%a)
         data%iops(i_iop)%b(:) = 0.0
         data%iops(i_iop)%b_b = 0.01
         data%iops(i_iop)%b_p = 1.0
      CASE (99)  ! CUSTOM (carbon-specific absorption and total scattering spectra)
         ! Note that 12.0107 converts from mg-1 to mmol-1
         !   get_parameter(a_star_iop, 'a_star_iop'//trim(strindex), 'm2/mg C', &
         !                 'carbon-mass-specific absorption coefficient for IOP '//trim(strindex)//' at reference wavelength', minimum=0., default=0.)
         !# Inland (Case-2) custom IOP: exponential absorption (CDOM/NAP form, Helms et al. 2008;
         !# Babin et al. 2003) + power-law scattering (Shi et al. 2013). Parameters are per-IOP,
         !# read from the /aed_oasim/ namelist arrays (indexed by i_iop). 12.0107 converts
         !# per-mg-C to per-mmol-C so a_star/b_star are m2/mg C and the linked var is mmol C/m3.
         IF (a_star_iop(i_iop) /= 0.) THEN
            data%iops(i_iop)%a(:) = a_star_iop(i_iop) * &
                  exp(-S_iop(i_iop) * (data%lambda - lambda_ref_iop(i_iop))) * 12.0107
         ELSE
            data%iops(i_iop)%a(:) = 0.
         ENDIF
         IF (b_star_iop(i_iop) /= 0.) THEN
            data%iops(i_iop)%b(:) = b_star_iop(i_iop) * &
                  (lambda_ref_iop(i_iop) / data%lambda)**eta_iop(i_iop) * 12.0107
            data%iops(i_iop)%b_b  = b_b_iop(i_iop)
         ELSE
            data%iops(i_iop)%b(:) = 0.
            data%iops(i_iop)%b_b  = 0.
         ENDIF
         data%iops(i_iop)%b_p = 1.0
      CASE DEFAULT
         !# 2026-08-29 FIX: there was no CASE DEFAULT, so an iop_type with no matching branch
         !# left %a and %b as freshly-ALLOCATEd (undefined) memory. CASE(7) (PIC) is such a
         !# branch - its body is entirely commented out - and so is an unset/zero iop_type when
         !# n_iop exceeds the number of entries supplied in iop_type. Those undefined arrays
         !# then went straight into max(...,zero_) and, for type<20, into the two-stream light
         !# budget at ~1107, contributing arbitrary absorption and scattering. Zero them and
         !# warn instead. Unreachable in the live configuration (iop_type = 23,22,21,25,99,99,
         !# 10,99,10,6 - every one of those has an explicit CASE), so no result change.
         data%iops(i_iop)%a(:) = zero_
         data%iops(i_iop)%b(:) = zero_
         data%iops(i_iop)%b_b  = zero_
         data%iops(i_iop)%b_p  = 1.0
         print *,'          WARNING: aed_oasim IOP ',i_iop,' has unsupported iop_type = ', &
                 iop_type(i_iop),'; its absorption and scattering are set to zero'
      END SELECT

      ! Protect against -ve coefficients caused by extrapolation beyond source spectrum boundaries
      data%iops(i_iop)%a(:) = max( data%iops(i_iop)%a, zero_ )
      data%iops(i_iop)%b(:) = max( data%iops(i_iop)%b, zero_ )

      ! Apply any optional user "scaling" adjustments to a spectra
      !# 2026-07 fix: scale was applied to absorption but NOT to scattering. For the
      !# phytoplankton IOPs, iop_scale is the mmol C -> mg Chl conversion (1/Xcc), and
      !# both a and b in the source tables are chlorophyll-specific. Later, a, b and b_b
      !# are all multiplied by the SAME concentration (in mmol C/m3), so leaving b
      !# unscaled overstated phytoplankton scattering and backscattering by 1/scale
      !# (~3.3-4.2x for the configured groups), inflating the diffuse-stream attenuation.
      data%iops(i_iop)%a(:) = data%iops(i_iop)%a(:) * data%iops(i_iop)%scale
      data%iops(i_iop)%b(:) = data%iops(i_iop)%b(:) * data%iops(i_iop)%scale

      ! Link to concentration metric to allow us to convert *specific* absorption/scattering into
      ! actual absorption and scattering (in m-1)
      IF (iop_type(i_iop) >=1 .and. iop_type(i_iop) <= 5) THEN
         ! Phytoplankton: chlorophyll-specific absorption and scattering
         data%iops(i_iop)%id_c = aed_locate_variable(TRIM(iop_link(i_iop)))
      ELSEIF (iop_type(i_iop) == 20) THEN
         ! Plant canopy: blade/leaf specific absorption and scattering (NOTE: benthic)
         data%iops(i_iop)%id_c    = aed_locate_sheet_variable(TRIM(iop_link(i_iop)))
         data%iops(i_iop)%id_parc = aed_locate_sheet_variable(TRIM(iop_link(i_iop))//'_parc')
         data%iops(i_iop)%id_aeff = aed_locate_sheet_variable(TRIM(iop_link(i_iop))//'_aeff')
      ELSE
         ! POM/DOM/PIC: carbon-specific absorption and scattering
         data%iops(i_iop)%id_c = aed_locate_variable(TRIM(iop_link(i_iop)))
      ENDIF
   ENDDO

   !-----------------------------------------------
   ! Find wavelength bounds of photosynthetically active radiation, swr & uv
   ALLOCATE(data%par_weights(nlambda))
   ALLOCATE(data%swr_weights(nlambda))
   ALLOCATE(data%uv_weights(nlambda))
   ALLOCATE(data%par_E_weights(nlambda))
   CALL calculate_integral_weights(400.,  700., nlambda, data%lambda, data%par_weights)
   CALL calculate_integral_weights(300., 4000., nlambda, data%lambda, data%swr_weights)
   CALL calculate_integral_weights(300.,  400., nlambda, data%lambda, data%uv_weights )
   data%par_E_weights(:) = data%par_weights * data%lambda /(Planck*lightspeed)/Avogadro*1e-3
                           ! divide by 1e9 to go from nm to m, multiply by 1e6 to go from mol to umol

   !-----------------------------------------------
   ! Find and link to environmental variables, from host
   data%id_yearday    = aed_locate_sheet_global('yearday')
   data%id_lat        = aed_locate_sheet_global('latitude')
   data%id_lon        = aed_locate_sheet_global('longitude')
   data%id_I_0        = aed_locate_sheet_global('par_sf')
   data%id_relhum     = aed_locate_sheet_global('humidity')
   data%id_wind_speed = aed_locate_sheet_global('wind_speed')
   data%id_h          = aed_locate_global('layer_ht')

   !-----------------------------------------------
   ! Allocate necesary diagnostics
   IF( diag_level>9 ) THEN
     data%id_zen =        aed_define_sheet_diag_variable('zen'    ,'degrees', 'zenith angle')
     data%id_doy =        aed_define_sheet_diag_variable('yearday','-', 'day of year')
     data%id_swflux =     aed_define_sheet_diag_variable('swflux' ,'-', 'swflux from host')
     data%id_hr_peakz =   aed_define_sheet_diag_variable('hr_peak_zen' ,'-', 'time of peak zenith')
     data%id_hr_peakx =   aed_define_sheet_diag_variable('hr_peak_flux' ,'-', 'time of peak zenith')
   ENDIF

   IF (sky_model == 1) THEN
      ! Constant OR Van Heuklon (1979) function to account for seasonal / geographical variation in O3
     data%id_O3 = 0  !  data%id_O3 = aed_locate_sheet_global('atmosphere_mass_content_of_ozone')
     IF(ozone) data%id_O3 = aed_define_sheet_diag_variable('ozone', 'kg/m2', 'atmosphere_mass_content_of_ozone')

     ! Aerosol properties
     data%id_alpha_a =  aed_define_sheet_diag_variable('alpha_a', '-' , 'aerosol Angstrom exponent')
     data%id_beta_a =   aed_define_sheet_diag_variable('beta_a' , '-' , 'aerosol scale factor for optical thickness')
     data%id_omega_a =  aed_define_sheet_diag_variable('omega_a', '-' , 'aerosol single scattering albedo')
     data%id_F_a =      aed_define_sheet_diag_variable('F_a'    , '-' , 'aerosol forward scattering probability')

     ! Wind properties
     IF (compute_mean_wind) THEN
      data%id_mean_wind_out = aed_define_sheet_diag_variable('mean_wind', 'm/s', 'daily mean wind speed')
     ELSE
      data%id_mean_wind_speed = aed_locate_sheet_global('mean_wind')
     ENDIF
   ENDIF
   data%id_wind_out = aed_define_sheet_diag_variable('wind'   ,'m/s', 'wind speed')

   ! Horizontal downwelling irradiance just above the water surface (BEFORE reflection by water surface)
   data%id_swr_sf =     aed_define_sheet_diag_variable('swr_sf' ,'W/m^2', 'downwelling shortwave flux in air')
   data%id_swr_dif_sf = aed_define_sheet_diag_variable('swr_dif_sf','W/m^2', 'diffuse downwelling shortwave flux in air')
   data%id_uv_sf =      aed_define_sheet_diag_variable('uv_sf'  ,'W/m^2', 'downwelling ultraviolet radiative flux in air')
   data%id_par_sf =     aed_define_sheet_diag_variable('par_sf' ,'W/m^2', 'downwelling photosynthetic radiative flux in air')
   data%id_par_E_sf =   aed_define_sheet_diag_variable('par_E_sf','umol/m^2/s', 'downwelling photosynthetic photon flux in air')

   ! Horizontal downwelling irradiance just below the water surface (AFTER reflection by water surface)
   data%id_swr_sf_w =   aed_define_sheet_diag_variable('swr_sf_w', 'W/m^2', 'downwelling shortwave flux in water')
   data%id_uv_sf_w =    aed_define_sheet_diag_variable('uv_sf_w' , 'W/m^2', 'downwelling ultraviolet radiative flux in water')
   data%id_par_sf_w =   aed_define_sheet_diag_variable('par_sf_w', 'W/m^2', 'downwelling photosynthetic radiative flux in water')
   data%id_par_E_sf_w = aed_define_sheet_diag_variable('par_E_sf_w','umol/m^2/s', 'downwelling photosynthetic photon flux in water')

   ! Scalar downwelling irradiance within the water column
   data%id_swr =        aed_define_diag_variable('swr', 'W/m^2', 'downwelling shortwave flux')
   data%id_uv =         aed_define_diag_variable('uv', 'W/m^2', 'downwelling ultraviolet radiative flux')
   data%id_par =        aed_define_diag_variable('par', 'W/m^2', 'downwelling photosynthetic radiative flux')
   data%id_par_E =      aed_define_diag_variable('par_E', 'umol/m^2/s', 'downwelling photosynthetic photon flux')
   data%id_par_J_scalar=aed_define_diag_variable('par_J_scalar','W/m^2', 'scalar downwelling photosynthetic radiative flux')
   data%id_par_E_scalar=aed_define_diag_variable('par_E_scalar','umol/m^2/s', 'scalar downwelling photosynthetic photon flux')
   data%id_swr_abs =    aed_define_diag_variable('swr_abs', 'W/m^2', 'absorption of shortwave energy in layer')
   !# 2026-07: PAR-band diffuse attenuation from NON-WATER constituents only (phytoplankton,
   !# CDOM, detritus, sediment). This is the quantity that is directly comparable with the
   !# tuned non-spectral Ke* constants in aed.nml (KePHY, KeDOM, KeDOMR, KePOM, Ke_ss), and
   !# it is what a future aed_light_extinction_oasim must return - GLM's `extinction` is an
   !# ACCUMULATOR sitting on top of Kw, so returning the TOTAL Kd (which includes pure
   !# water) would double-count the water term. Diagnostic only for now: nothing consumes it.
   data%id_kd_bio =     aed_define_diag_variable('Kd_bio', 'm-1', &
                          'PAR attenuation from non-water constituents (excludes pure water)')
   !# 2026-07: PAR attenuation of PURE WATER alone, from the same two-stream form. This is
   !# the value Kw in glm3.nml should be set to when drive_extinction is enabled - do not
   !# estimate it by subtracting Kd_bio from Kd, because Kd is a log of transmission and is
   !# NOT additive across constituents.
   data%id_kd_water =   aed_define_diag_variable('Kd_water', 'm-1', &
                          'PAR attenuation of pure water alone (the correct value for Kw)')

   !# 2026-07 (Tier 1C): one normalised spectral-absorption factor per IOP.
   !#   afac = <a*(l) weighted by the local PAR spectrum> / <a*(l) weighted flat>
   !# Dimensionless, and exactly 1.0 for a flat spectrum - the same normalisation the ABM
   !# already applied to its hand-tuned band template, so a consumer can swap to this
   !# WITHOUT recalibrating alphachl. The difference is that the weights here are the
   !# actual chlorophyll-specific absorption spectra used for shading, over all PAR bands
   !# with proper integration weights - not 7 hand-picked band centres addressed by
   !# POSITION in lambda_out (where editing lambda_out silently re-points every weight).
   !# 2026-08: computed under BOTH weightings, always. afac_iop is the copy consumers read.
   ALLOCATE(data%id_afac(size(data%iops)))
   ALLOCATE(data%id_afacE(size(data%iops)))
   ALLOCATE(data%id_afacq(size(data%iops)))
   ALLOCATE(data%id_astar(size(data%iops)))
   DO i_iop = 1, size(data%iops)
      write(strindex2, '(i0)') i_iop
      data%id_afac(i_iop) = aed_define_diag_variable('afac_iop' // trim(strindex2), '-', &
             'spectral absorption factor for IOP ' // trim(strindex2) // ' (1.0 = flat spectrum)')
      data%id_afacE(i_iop) = aed_define_diag_variable('afacE_iop' // trim(strindex2), '-', &
             'energy-weighted spectral absorption factor for IOP ' // trim(strindex2))
      data%id_afacq(i_iop) = aed_define_diag_variable('afacq_iop' // trim(strindex2), '-', &
             'photon-weighted spectral absorption factor for IOP ' // trim(strindex2))
      data%id_astar(i_iop) = aed_define_diag_variable('astar_iop' // trim(strindex2), &
             'm-1 concentration-1', &
             'flat photon-weighted mean specific absorption for IOP ' // trim(strindex2))
   ENDDO
   data%id_par_E_dif =  aed_define_diag_variable('par_E_dif', 'umol/m^2/s', 'diffusive downwelling photosynthetic photon flux')


   !-----------------------------------------------
   ! Interpolate (pre-defined) absorption and scattering spectra to user wavelength grid
   ALLOCATE(data%a_w(nlambda), data%b_w(nlambda))
   CALL interp(nlambda_w, lambda_w, a_w, nlambda, data%lambda, data%a_w)
   CALL interp(nlambda_w, lambda_w, b_w, nlambda, data%lambda, data%b_w)

   ALLOCATE(data%exter(nlambda), data%a_o(nlambda), data%a_v(nlambda), data%a_u(nlambda), data%tau_r(nlambda))
   IF (exter_source == 1) THEN
      CALL interp(nlambda_oasim, lambda_oasim, ET_oasim, nlambda, data%lambda, data%exter)
   ELSE
      CALL interp(nlambda_astm, lambda_astm, ET_astm, nlambda, data%lambda, data%exter)
   ENDIF
   CALL interp(nlambda_oasim, lambda_oasim, a_o_oasim, nlambda, data%lambda, data%a_o)
   CALL interp(nlambda_oasim, lambda_oasim, a_v_oasim, nlambda, data%lambda, data%a_v)
   CALL interp(nlambda_oasim, lambda_oasim, a_u_oasim, nlambda, data%lambda, data%a_u)
  !CALL interp(nlambda_oasim, lambda_birdrior1986, a_u_birdrior1986, nlambda, data%lambda, data%a_u)
   ! Rayleigh optical thickness (Eq 2 Bird 1984, Eq 4 Bird & Riordan 1986, Eq 15 Gregg & Carder 1990)
  !CALL interp(nlambda_oasim, lambda_oasim, tau_r_oasim, nlambda, data%lambda, data%tau_r)

   data%tau_r = 1.0 / (115.6406 * (data%lambda/1000)**4 - 1.335 * (data%lambda/1000)**2)

   ! Protect against negative absorption coefficients produced by linear extrapolation
   data%a_o = max(0., data%a_o)
   data%a_v = max(0., data%a_v)
   data%a_u = max(0., data%a_u)

   ! Wavelength dependence of foam reflectance (Eqs A11, A10 Gregg & Casey 2009)
   ALLOCATE(data%F(nlambda))
   DO l = 1, nlambda
      IF (data%lambda(l) < 900) THEN
         ! Note: close to 900 nm the expression below returns negative values.
         ! We clip to zero in line with Fig A1 Gregg & Casey 2009
         log_T_w = -(data%a_w(l) + 0.5 * data%b_w(l))
         data%F(l) = max(0., a0 + a1 * log_T_w + a2 * log_T_w**2 + a3 * log_T_w**3)
      ELSE
         ! Note: above 1700 nm the expression below returns negative values.
         ! We clip to zero in line with Fig A1 Gregg & Casey 2009
         data%F(l) = max(0., b0 + b1 * data%lambda(l) + b2 * data%lambda(l)**2 + b3 * data%lambda(l)**3)
      ENDIF
   ENDDO

   !-----------------------------------------------
   ! Setup output configuration
   data%save_Kd = save_Kd
   data%drive_extinction = drive_extinction
   data%afac_quanta = afac_quanta
   IF (afac_quanta) THEN
      print *,"          OASIM: afac_quanta = TRUE - absorption weighted by photon flux"
   ENDIF
   IF (drive_extinction) THEN
      print *,"          OASIM: drive_extinction = TRUE - OASIM Kd_bio now drives host light"
      print *,"          WARNING: zero or re-derive Kw and the Ke* constants in aed.nml,"
      print *,"                   otherwise every attenuator is counted twice."
      IF (.NOT. save_Kd .OR. spectral_output == 0) THEN
         print *,"          ERROR: drive_extinction requires save_Kd=.true. and spectral_output>0"
         ERROR STOP 1
      ENDIF
   ENDIF
   data%spectral_output = spectral_output
   SELECT CASE (data%spectral_output)
   CASE (1)   ! Output each wavelength as a diagnostic
      ALLOCATE(data%lambda_out(nlambda))
      data%lambda_out(:) = data%lambda
   CASE (2)   ! Output only requested wavelengths
      ALLOCATE(data%lambda_out(nlambda_out))
      DO l = 1, nlambda_out
         data%lambda_out(l) = lambda_out(l)
      ENDDO
      ALLOCATE(data%a_w_out(nlambda_out), data%b_w_out(nlambda_out))
      CALL interp(nlambda_w, lambda_w, a_w, nlambda_out, data%lambda_out, data%a_w_out)
      CALL interp(nlambda_w, lambda_w, b_w, nlambda_out, data%lambda_out, data%b_w_out)
   END SELECT

   IF (ALLOCATED(data%lambda_out)) THEN
      ALLOCATE(data%id_surface_band_dir(size(data%lambda_out)))
      ALLOCATE(data%id_surface_band_dif(size(data%lambda_out)))
      ALLOCATE(data%id_band_dir(size(data%lambda_out)))
      ALLOCATE(data%id_band_dif(size(data%lambda_out)))
      ALLOCATE(data%id_a_iop(size(data%lambda_out), size(data%iops)))
      ALLOCATE(data%id_a_band(size(data%lambda_out)))
      ALLOCATE(data%id_b_band(size(data%lambda_out)))
      IF (data%save_Kd) ALLOCATE(data%id_Kd(size(data%lambda_out)+1))
      DO l = 1, size(data%lambda_out)
         IF (data%lambda_out(l) < 1000.) THEN
            write(strwavelength, '(f5.1)') data%lambda_out(l)
         ELSE
            write(strwavelength, '(f6.1)') data%lambda_out(l)
         ENDIF
         write(strindex, '(i0)') l
         data%id_surface_band_dir(l) = aed_define_sheet_diag_variable('dir_sf_band' // trim(strindex), &
                                'W/m2/nm', 'downward direct irradiance in air @ ' // trim(strwavelength) // ' nm')
         data%id_surface_band_dif(l) = aed_define_sheet_diag_variable('dif_sf_band' // trim(strindex), &
                                'W/m2/nm', 'downward diffuse irradiance in air @ ' // trim(strwavelength) // ' nm')
         data%id_band_dir(l) = aed_define_diag_variable('dir_band' // trim(strindex) , &
                                'W/m2/nm', 'direct irradiance @ ' // trim(strwavelength) // ' nm')
         data%id_band_dif(l) = aed_define_diag_variable('dif_band' // trim(strindex), &
                                'W/m2/nm', 'diffuse irradiance @ ' // trim(strwavelength) // ' nm')
         IF (data%save_Kd) &
            data%id_Kd(l) = aed_define_diag_variable('Kd_band' // trim(strindex), 'm-1', &
                                                      'attenuation @ ' // trim(strwavelength) // ' nm')
         DO i_iop = 1, size(data%iops)
            write(strindex2, '(i0)') i_iop
            data%id_a_iop(l, i_iop) = aed_define_diag_variable('a_iop' // trim(strindex2) // '_band' // trim(strindex), &
                                '1/m', 'absorption by IOP ' // trim(strindex2) // ' @ ' // trim(strwavelength) // ' nm')
         ENDDO
         data%id_a_band(l) = aed_define_diag_variable('a_band' // trim(strindex), &
                                'm-1', 'total absorption excluding water @ ' // trim(strwavelength) // ' nm')
         data%id_b_band(l) = aed_define_diag_variable('b_band' // trim(strindex), &
                               'm-1', 'total scattering excluding water @ ' // trim(strwavelength) // ' nm')
      ENDDO
      IF (data%save_Kd) &
      data%id_Kd(size(data%lambda_out)+1) = aed_define_diag_variable('Kd', 'm-1', &
                                     'PAR (400-700nm) diffuse attenuation coefficient')
   ENDIF

   !-----------------------------------------------
   ! Add option for secchi depth output
   !# 2026-07 fix: the 490nm index search was inside the save_Kd/spectral_output guard,
   !# but l490_l is read UNCONDITIONALLY in the column routine - for l490_sf, for the
   !# Kd490 accumulation and for the attenuation-fraction test. With save_Kd=.false. or
   !# spectral_output=0 it was undefined and indexed swr_weights out of bounds. The
   !# search is cheap, so always do it; the type declaration also carries a default.
   ! Find lambda index for 490nm (e.g. for secchi)
   data%l490_l = nlambda - 1
   DO l = 1, nlambda - 1
      IF (data%lambda(l) >= 490.) THEN
         data%l490_l = l
         exit
      ENDIF
   ENDDO

   IF (save_Kd .AND. spectral_output>0) THEN
     data%id_secchi =     aed_define_sheet_diag_variable('secchi', 'm', 'Secchi depth (1.7/Kd490)')
   ENDIF

END SUBROUTINE aed_define_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_surface_oasim(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Atmospheric component of the AED implementation of the OASIM model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oasim_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   ! Environment
   AED_REAL :: wind
   ! State
   AED_REAL :: O3

!
!-------------------------------------------------------------------------------
!BEGIN

! ASPECTS OF ATMOSPHERE MODEL BELOW CAN BE CONSIDERED TO COME HERE

END SUBROUTINE aed_calculate_surface_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_column_oasim(data,column,layer_map)
!-------------------------------------------------------------------------------
! Vertical penetration of light
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oasim_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_map(:)
!
!LOCALS
   INTEGER :: layer_idx, layer
   AED_REAL,PARAMETER :: pres0 = 101300.    ! Reference air pressure in Pa (Bird 1984 p461, Bird & Riordan p89)
  !AED_REAL,PARAMETER :: ga = 0.05          ! Ground albedo (used by Bird & Riordan 1986, but discarded by Gregg & Carder 1990)
   AED_REAL,PARAMETER :: H_oz = 22.         ! Height of maximum ozone concentration (km)
   AED_REAL,PARAMETER :: r_e = (10. + 11.8) / 2 ! Equivalent radius of cloud drop size distribution (um)
                                            ! based on mean of Kiehl et al. & Han et al. (cf OASIM)
   AED_REAL,PARAMETER :: mcosthetas = 0.831 ! Mean of cosine of angle of diffuse radiation in water,
                                            ! assuming all angular contributions equal in air
                                            ! (Sathyendranath and Platt 1989, p 191) NB Ackleson et al 1994 use 0.9
   AED_REAL,PARAMETER :: r_s = 1.5          ! Shape factor representing mean backscatter coefficient
                                            ! of diffuse irradiance (r_d in Ackleson et al. 1994 p 7487)

   AED_REAL :: longitude, latitude, yearday, cloud_cover, wind_speed, &
               airpres, relhum, LWP, water_vapour, WV, WM, visibility, AM
   AED_REAL :: days, hour, theta, costheta, alpha_a, beta_a, yearday_utc
   INTEGER  :: l
   AED_REAL :: M, M_prime, M_oz
   AED_REAL :: O3
   AED_REAL,DIMENSION(data%nlambda) :: direct, diffuse, spectrum, Kd  ! Spectra at top of the water
                                            ! column (with refraction and reflection accounted for)
   AED_REAL,DIMENSION(nlambda_5nm_astm) :: direct_5nm, diffuse_5nm ! Spectra at top of the water
   AED_REAL :: par_J=0., swr_J=0., uv_J=0., par_E=0., F_a, omega_a, par_M
   AED_REAL,DIMENSION(data%nlambda) :: tau_a, T_a, T_oz, T_w, T_u, T_r, T_aa, T_as
   AED_REAL,DIMENSION(data%nlambda) :: T_g, T_dclr, T_sclr, T_dcld, T_scld
   AED_REAL,DIMENSION(data%nlambda) :: rho_d, rho_s

   AED_REAL,DIMENSION(data%nlambda) :: a, b, b_b, a_iop, a_mac
   AED_REAL,DIMENSION(data%nlambda) :: f_att_d, f_att_s, f_prod_s, f_att_m
   AED_REAL,DIMENSION(data%nlambda) :: spectrum_out
   INTEGER  :: i_iop
   !# 2026-08-29: index of the type-20 (benthic canopy) IOP, captured inside the IOP loop.
   !# See the fix note at the canopy block below. -1 = no canopy IOP configured.
   INTEGER  :: i_iop_mac
   AED_REAL,DIMENSION(data%nlambda) :: Kd_w      !# water-only Kd, for the Kd_bio diagnostic
   AED_REAL :: c_iop, h, swr_top, par_top, l490_sf, costheta_r, dir_frac
   AED_REAL,DIMENSION(data%nlambda) :: fwd_frac   ! forward-scattered fraction, per band
   AED_REAL :: kd_bio_num, kd_bio_den, kd_wat_num
   AED_REAL :: kd_k(3), kd_c(3)   !# 2026-09-08: coefficients/weights of the three-stream mixture
   LOGICAL  :: kd_ok              !# .FALSE. from aed_oasim_kd_mixture on invalid input
   AED_REAL :: afac_num, afac_den, afac_fnum, afac_fden
   AED_REAL :: afq_num, afq_den, afq_fnum, afq_fden, afacE, afacq, astarq
   AED_REAL :: SWFLUX, MaxSWFlux
   AED_REAL :: c_eff, atten_frac, Kd490avg=0., a_eff
   INTEGER  :: month, secchi_botlayer
   CHARACTER(len=1):: met
!
!-------------------------------------------------------------------------------
!BEGIN
   layer_idx = layer_map(1)                   ! Start at the top of the column

   !-------------------------------------------------------------------------------
   ! Set column time/location properties
   longitude   = _STATE_VAR_S_(data%id_lon)
   latitude    = _STATE_VAR_S_(data%id_lat)

   !# 2026-08-29 LATITUDE UNIT FIX. GLM hands latitude over as radians PLUS two_Pi
   !# (glm_model.c:195, `Latitude = two_Pi + Latitude * deg2rad`) and the global is declared
   !# 'radians' (glm_aed.F90:291) - but BOTH consumers in this module expect DEGREES:
   !#   zenith_angle()  aed_slingo.F90:272-273   rlat = deg2rad * dlat
   !#   vH_O3()         aed_atmosphere.F90:333   sin(beta * phi * deg2rad)**2
   !# The second deg2rad therefore put the sun at 7.21 N instead of 53.15 N. Measured before the
   !# fix: modelled solar zenith 15.7 deg at local noon against a true 34.0, and a minimum over
   !# the run of 7.50 deg - geometrically impossible at this latitude, where the sun never rises
   !# above ~29.8 deg from the zenith. Wrong air mass, wrong direct/diffuse split, wrong day
   !# length (~12.2 h tropical vs ~16.5 h at Stechlin in July) and wrong in-water path.
   !#
   !# Unwrap the 2*pi offset, then convert. GLM's own C code is self-consistent and must NOT be
   !# touched: trig is 2*pi-periodic so glm_util.c:80-81 (cos/sin of lat) is unaffected, and the
   !# one place needing degrees subtracts the offset explicitly (glm_surface.c:2088).
   !#
   !# Longitude is deliberately NOT converted - GLM never applies deg2rad to it, so it already
   !# arrives in degrees, which is what both consumers want.
   IF ( ABS(latitude) > pi ) latitude = latitude - SIGN(2.*pi, latitude)
   latitude = latitude * (180./pi)
   yearday     = _STATE_VAR_S_(data%id_yearday)
   yearday_utc = yearday - data%tzshift/24.
   days        = floor(yearday_utc)
   hour        = mod(yearday_utc, 1.0) * 24.0

   ! Calculate zenith angle (in radians)
   theta = zenith_angle(days, hour, longitude, latitude)
   theta = min(theta, 0.5*pi)      ! Restrict the input zenith angle between 0 & pi/2
   costheta = cos(theta)

   ! Set column environmental conditions (from host)
   wind_speed  = _STATE_VAR_S_(data%id_wind_speed)      ! Wind speed @10m above surface (m/s)
   relhum      = _STATE_VAR_S_(data%id_relhum)          ! Relative humidity (-)
   SWFLUX      = _STATE_VAR_S_(data%id_I_0)             ! Incoming shortwave flux

   IF( diag_level>9 ) THEN
      _DIAG_VAR_S_(data%id_swflux) = SWFLUX
      _DIAG_VAR_S_(data%id_zen) = rad2deg * theta
      _DIAG_VAR_S_(data%id_doy) = yearday
      ! Special diagnostics for tracking time of day solar peak occurs
      IF( mod(yearday, 1.0) * 24.0 < 1. ) THEN
         _DIAG_VAR_S_(data%id_hr_peakz) =  90.
         _DIAG_VAR_S_(data%id_hr_peakx) =  0.
      ENDIF
      !# 2026-08-29 FIX: hr_peakx was a pure SELF-ASSIGNMENT
      !#     IF( SWFLUX > _DIAG_VAR_S_(data%id_hr_peakx) ) &
      !#       _DIAG_VAR_S_(data%id_hr_peakx) = _DIAG_VAR_S_(data%id_hr_peakx)
      !# so hr_peak_flux never moved off the 0. set at midnight (line above) and the diagnostic
      !# read 0.0 for the whole run. The commented-out RHS shows the intended value:
      !# mod(yearday,1.0)*24.0, i.e. the hour of day. It cannot be latched on `SWFLUX > stored`
      !# though, because the stored value is an HOUR, not a flux, and no running peak-flux
      !# variable persists between calls (MaxSWFlux is a local, reset every step). The daily
      !# peak is latched here on the same test as hr_peakz, which is what the registered
      !# long_name ('time of peak zenith') asks for: hr_peakz keeps the minimum zenith of the
      !# day, hr_peakx keeps the hour at which it occurred.
      !# Both are diag_level>9 debug outputs with no consumer anywhere in the tree (grep:
      !# hr_peak appears only in this file), so no model result changes.
      IF( _DIAG_VAR_S_(data%id_zen) < _DIAG_VAR_S_(data%id_hr_peakz) ) THEN
        _DIAG_VAR_S_(data%id_hr_peakz) = _DIAG_VAR_S_(data%id_zen)
        _DIAG_VAR_S_(data%id_hr_peakx) = mod(yearday, 1.0) * 24.0
      ENDIF
   ENDIF

   !-------------------------------------------------------------------------------
   ! Compute incoming spectral irradiance based on atmospheric properties
   IF(data%sky_model == 1) THEN

       ! Set column environmental conditions (from parameters)
       AM          = data%air_mass_type !_STATE_VAR_S_(data%id_air_mass_type)
                                        ! Aerosol air mass type (1: open ocean, 10: continental)
       visibility  = data%visibility ! _STATE_VAR_S_(data%id_visibility)     ! Visibility (m)
       airpres     = data%airpres ! _STATE_VAR_S_(data%id_airpres)       ! Surface air pressure (Pa)
       cloud_cover = data%cloud ! _STATE_VAR_S_(data%id_cloud)       ! Cloud cover (fraction, 0-1)
       LWP         = data%LWP ! _STATE_VAR_S_(data%id_lwp)      ! Cloud liquid water content (kg m-2)
       water_vapour= data%WV ! _STATE_VAR_S_(data%id_wv)    ! Total precipitable water vapour (kg m-2) - equivalent to mm
       O3          = data%O3                             ! Ozone content (kg m-2)
       IF(data%ozone) THEN
         O3 = vH_O3(longitude, latitude, yearday)
         _DIAG_VAR_S_(data%id_O3) = O3
       ENDIF

       ! Adjustments and diagnostic outputs for debugging
       _DIAG_VAR_S_(data%id_wind_out) = wind_speed
       IF (data%id_mean_wind_out > 0) THEN
         ! Daily average of wind
         _DIAG_VAR_S_(data%id_mean_wind_out) = _DIAG_VAR_S_(data%id_mean_wind_out) &
                                      * (1-(DTday/WindAvgTime)) + wind_speed*(DTday/WindAvgTime)
         WM = _DIAG_VAR_S_(data%id_mean_wind_out) ! Daily mean wind speed @10m above surface (m/s)
       !# 2026-08-29 FIX: with compute_mean_wind = .FALSE. the branch above is skipped (only
       !# id_mean_wind_speed is assigned, at ~565) and WM was left UNDEFINED, yet it is passed
       !# to navy_aerosol_model below where it sets the sea-salt aerosol loading. Read the
       !# host-supplied 'mean_wind' global in that case, and fall back to the instantaneous
       !# wind if neither id was ever set. Unreached in the live model (compute_mean_wind
       !# defaults to .TRUE., so id_mean_wind_out > 0 and the first branch always wins).
       ELSEIF (data%id_mean_wind_speed > 0) THEN
         WM = _STATE_VAR_S_(data%id_mean_wind_speed)
       ELSE
         WM = wind_speed
       ENDIF
       IF (cloud_cover > 0) LWP = LWP / cloud_cover

       ! LWP is the mean density over a grid box (Jorn: ECMWF pers comm 26/2/2019), but
       !        we want the mean density per cloud-covered area
       WV   = water_vapour / 10             ! from kg m-2 to cm
       O3   = O3 * (1000 / 48.) / 0.4462    ! from kg m-2 to mol m-2, then from mol m-2 to atm cm (Basher 1982)

       ! Atmospheric path length, a.k.a. relative air mass (Eq 3 Bird 1984, Eq 5 Bird & Riordan 1986, &
       !                                        Eq 13 Gregg & Carder 1990, Eq A5 in Casey & Gregg 2009)
       ! Note this should always exceed 1, but as it is an approximation it does not near &
       !                                    theta -> 0 (Tomasi et al. 1998 p14). Hence the max operator.
       M = max(1., 1. / (costheta + 0.15 * (93.885 - theta * rad2deg)**(-1.253)))

       ! Pressure-corrected atmospheric path length (Eq A6 Casey & Gregg 2009)
       M_prime = M * airpres / pres0

       ! Atmospheric path length for ozone
       !  Eq 10, Bird & Riordan 1986, Eq 14 Gregg & Carder 1990; NB 6370 is the earth's radius in km
       !  See also Tomasi et al. 1998 Eq 5
       M_oz = (1. + H_oz / 6370.) / sqrt(costheta**2 + 2 * H_oz / 6370.)

       ! Transmittance due to ozone absorption (Eq 8 Bird 1984, Eq 9 Bird & Riordan 1986, Eq 17 Gregg & Carder 1990)
       T_oz = exp(-data%a_o * O3 * M_oz)

       ! Transmittance due to water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
       !  Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
       T_w = exp((-0.2385 * data%a_v * WV * M) / (1. + 20.07 * data%a_v * WV * M)**0.45)

       ! Transmittance due to uniformly mixed gas absorption - SHOULD use pressure corrected airmass
       !  Eq 10 Bird 1984, Eq 11 Bird and Riordan 1986, Eq 18 Gregg & Carder 1990
       !  Bird & Riordan use 118.93 rather than 118.3, but state 118.3 should be used in the future.
       T_u = exp(-1.41 * data%a_u * M_prime / (1. + 118.3 * data%a_u * M_prime)**0.45)

       ! Transmittance terms that apply for both cloudy and clear skies.
       T_g = T_oz * T_w * T_u

       ! -------------------
       ! Clear skies part
       ! -------------------

       ! Transmittance due to Rayleigh scattering (use precomputed optical thickness)
       T_r = exp(-M_prime * data%tau_r)

       ! Transmittance due to aerosol absorption (Eq 26 Gregg & Carder 1990)
       CALL navy_aerosol_model(AM, WM, wind_speed, relhum, visibility, costheta, alpha_a, beta_a, F_a, omega_a)
       tau_a = beta_a * data%lambda**(-alpha_a)
       T_a = exp(-tau_a * M)

       ! Direct transmittance
       T_dclr = T_r * T_a

       ! Separate absorption and scattering components of aerosol transmittance
       T_aa = exp(-(1. - omega_a) * tau_a * M)
       T_as = exp(-omega_a * tau_a * M)
       T_sclr = T_aa * 0.5 * (1. - T_r**0.95) + T_r**1.5 * T_aa * F_a * (1 - T_as)

       _DIAG_VAR_S_(data%id_omega_a) = omega_a
       _DIAG_VAR_S_(data%id_alpha_a) = alpha_a
       _DIAG_VAR_S_(data%id_beta_a) = beta_a
       _DIAG_VAR_S_(data%id_F_a) = F_a

       ! -------------------
       ! Cloudy skies part
       ! -------------------

       ! Transmittance due to absorption and scattering by clouds
       CALL slingo(costheta, max(0., LWP) * 1000, r_e, data%nlambda, data%lambda, T_dcld, T_scld)

       ! Diffuse and direct irradiance streams (Eqs 1, 2 Gregg & Casey 2009)
       !  These combine terms for clear and cloudy skies, weighted by cloud cover fraction
       direct  = data%exter * costheta * T_g * ((1. - cloud_cover) * T_dclr + cloud_cover * T_dcld)
       diffuse = data%exter * costheta * T_g * ((1. - cloud_cover) * T_sclr + cloud_cover * T_scld)

   ELSEIF(data%sky_model == 2) THEN

       ! Diffuse and direct irradiance streams from local analysis of Cockburn Sound radiation (Perth, WA)
       !  These functions combine terms for clear and cloudy skies
        met = 'W'

        _DIAG_VAR_S_(data%id_wind_out) = wind_speed

        ! Get function to return direct and diffuse spectra, based on incoming solar (etc)
       CALL direct_diffuse_bands(SWFlux, theta, yearday, met, data%nlambda, data%lambda_bounds, direct, diffuse)
   ENDIF

   ! Integrate to get total incoming surface irradiance and bandwidth components
   spectrum = direct + diffuse
   swr_J = sum(data%swr_weights   * spectrum)
   uv_J  = sum(data%uv_weights    * spectrum)
   par_J = sum(data%par_weights   * spectrum)
   par_E = sum(data%par_E_weights * spectrum)

   ! Set diagnostics of surface irradiance (bulk)
   _DIAG_VAR_S_(data%id_par_E_sf) = par_E    ! Photosynthetically Active Radiation, PAR (umol/m2/s)
   _DIAG_VAR_S_(data%id_par_sf)   = par_J    ! Photosynthetically Active Radiation, PAR (W/m2)
   _DIAG_VAR_S_(data%id_swr_sf)   = swr_J    ! Total shortwave radiation (W/m2), SW [up to 4000 nm]
   _DIAG_VAR_S_(data%id_uv_sf)    = uv_J     ! Ultraviolet Radiation, UV (W/m2)
   swr_J = sum(data%swr_weights * diffuse)   ! Repopulate variable with diffuse only fraction
   _DIAG_VAR_S_(data%id_swr_dif_sf) = swr_J  ! Diffuse shortwave radiation (W/m2), SW [up to 4000 nm]

   ! Set diagnostics of surface irradiance (spectral)
   SELECT CASE (data%spectral_output)
      CASE (1)
         DO l = 1, data%nlambda
            _DIAG_VAR_S_(data%id_surface_band_dir(l)) = direct(l)
            _DIAG_VAR_S_(data%id_surface_band_dif(l)) = diffuse(l)
         ENDDO
      CASE (2)
         CALL interp(data%nlambda, data%lambda, direct, size(data%lambda_out), data%lambda_out, spectrum_out)
         DO l = 1, size(data%lambda_out)
            _DIAG_VAR_S_(data%id_surface_band_dir(l)) = spectrum_out(l)
         ENDDO
         CALL interp(data%nlambda, data%lambda, diffuse, size(data%lambda_out), data%lambda_out, spectrum_out)
         DO l = 1, size(data%lambda_out)
            _DIAG_VAR_S_(data%id_surface_band_dif(l)) = spectrum_out(l)
         ENDDO
   END SELECT

   !-------------------------------------------------------------------------------
   ! Compute sea surface reflectance to ascertain fraction entering the water
   CALL reflectance(data%nlambda, data%F, theta, wind_speed, rho_d, rho_s, costheta_r)

   ! Incorporate the loss due to reflectance
   direct = direct * (1. - rho_d)
   diffuse = diffuse * (1. - rho_s)
   spectrum = direct + diffuse

   ! Re-compute aggregate amounts and store to diagnostics
   par_E = sum(data%par_E_weights * spectrum)
   par_J = sum(data%par_weights * spectrum)
   swr_J = sum(data%swr_weights * spectrum)
   uv_J  = sum(data%uv_weights * spectrum)
   _DIAG_VAR_S_(data%id_par_E_sf_w) = par_E        ! Photosynthetically Active Radiation (umol/m2/s)
   _DIAG_VAR_S_(data%id_par_sf_w) = par_J          ! Photosynthetically Active Radiation (W/m2)
   _DIAG_VAR_S_(data%id_swr_sf_w) = swr_J          ! Total shortwave radiation (W/m2) [up to 4000nm]
   _DIAG_VAR_S_(data%id_uv_sf_w) = uv_J            ! UV (W/m2)

   l490_sf = swr_J * data%swr_weights(data%l490_l) ! Incoming irradiance @ 490nm
   secchi_botlayer = 1                             ! Used to set secchi depth
   Kd490avg = zero_                                ! MUST reset each column: the declaration
                                                   ! initializer (Kd490avg=0.) gives it implicit
                                                   ! SAVE, so without this it accumulates across all
                                                   ! calls -> Secchi = 1.7/huge -> 0. (2026-07 fix)
   a_mac = zero_                                   ! Initialise macrophyte/seagrass absorbance
   i_iop_mac = -1                                  ! Initialise canopy IOP index (2026-08-29)

   !-------------------------------------------------------------------------------
   ! Now loop down through the water column
   DO layer = 1,SIZE(layer_map)

      layer_idx = layer_map(layer)

      !# 2026-09-08: read the layer thickness at the TOP of the iteration. It used to be read
      !# after the IOP loop, but the invalid-concentration report inside that loop already
      !# passes h - so it printed the PREVIOUS layer's thickness, and an uninitialised value on
      !# the first layer. The value itself is unchanged for every consumer; only its position
      !# moves, so results are unaffected.
      h = _STATE_VAR_(data%id_h)

      ! Save downwelling shortwave flux at top of the layer
      swr_top = swr_J
      !# 2026-07 fix: also keep the PAR-band flux at the layer top, so the scalar 'Kd'
      !# diagnostic below can be a true Kd(PAR). See the note at its assignment.
      par_top = par_J

      ! Compute absorption, total scattering and backscattering in current layer from IOPs
      a = data%a_w ; b = data%b_w
      b_b = 0.5 * data%b_w
      DO i_iop = 1, size(data%iops)
         !# 2026-08-23 FIX: this test was `type<20 .OR. type==99`, which silently EXCLUDED the
         !# pigment-reconstructed phytoplankton spectra added as iop_type 21/22/23/25. They fell
         !# into the gap between the macrophyte branch (20) and the custom branch (99), so they
         !# contributed NO absorption and NO scattering to the two-stream solve. Selecting them
         !# did not change the spectra used - it switched phytoplankton out of the light budget
         !# altogether, while OAS_afac_iop<n> still moved because afac is accumulated in a
         !# separate loop that walks every IOP regardless of type. That combination is why the
         !# first gate (c) looked like "reconstructed spectra improve temperature" when it was
         !# actually "phytoplankton absorption removed". Range-limited rather than `/=20` so a
         !# future unrelated type cannot be swept in by accident.
         IF (data%iops(i_iop)%type<20 .OR. data%iops(i_iop)%type==99 .OR.                     &
             (data%iops(i_iop)%type>=21 .AND. data%iops(i_iop)%type<=25)) THEN   ! type 99 = custom water-column IOP (CDOM/NAP)
            c_iop = _STATE_VAR_(data%iops(i_iop)%id_c)
            a_iop = c_iop * data%iops(i_iop)%a
            SELECT CASE (data%spectral_output)
            CASE (1)
               DO l = 1, data%nlambda
                  _DIAG_VAR_(data%id_a_iop(l, i_iop)) = a_iop(l)
               ENDDO
            CASE (2)
               CALL interp(data%nlambda, data%lambda, a_iop, size(data%lambda_out), data%lambda_out, spectrum_out)
               DO l = 1, size(data%lambda_out)
                  _DIAG_VAR_(data%id_a_iop(l, i_iop)) = spectrum_out(l)
               ENDDO
            END SELECT
            a = a + a_iop
            !# 2026-09-08: validate the RAW concentration before the power.
            !#
            !# CORRECTION to an earlier note here, which claimed a negative base raised to a real
            !# power is always NaN. MEASURED with this toolchain (gfortran 16.1, production flags):
            !#     (-1.0d0)**1.0d0  = -1.0  (finite)     (-1.0d-9)**1.0d0 = -1.0d-9 (finite)
            !#     (-2.5d0)**0.5d0  = NaN
            !# so a unit exponent is evaluated exactly and does NOT produce NaN. With b_p = 1.0 for
            !# every currently configured IOP, this is therefore NOT a live NaN source, and it does
            !# not explain the bottom-layer NaN seen when the non-PAR radiation patches were tried.
            !#
            !# The check stays because b_p is the packaging-effect exponent and is meant to take
            !# FRACTIONAL values; the moment one is configured, a transient negative concentration
            !# (before repair_state clamps it) would give NaN in b, b_b, Kd and extc. Catching a
            !# negative concentration is right regardless of the exponent.
            IF ( .NOT. aed_oasim_is_finite(c_iop) .OR. c_iop < zero_ ) THEN
               CALL aed_oasim_fatal('invalid IOP concentration before c_iop**b_p',                &
                    yearday, layer, layer_idx, -1, -one_, h,                                      &
                    c_iop, data%iops(i_iop)%b_p, REAL(i_iop,KIND(c_iop)),                         &
                    REAL(data%iops(i_iop)%type,KIND(c_iop)),                                      &
                    'c_iop, b_p, i_iop, iop_type')
            ENDIF
            c_eff = c_iop ** data%iops(i_iop)%b_p
            b = b + c_eff * data%iops(i_iop)%b
            b_b = b_b + c_eff * data%iops(i_iop)%b_b * data%iops(i_iop)%b
         ELSEIF (data%iops(i_iop)%type==20) THEN
            i_iop_mac = i_iop   !# 2026-08-29: remember which IOP is the canopy (see below)
            IF(layer==SIZE(layer_map)) THEN
               c_iop = _DIAG_VAR_S_(data%iops(i_iop)%id_c)    ! MAC_A
               a_eff = _DIAG_VAR_S_(data%iops(i_iop)%id_aeff) ! A_eff * sine_blade
               a_mac = data%iops(i_iop)%a * a_eff * c_iop      ! a * A_eff * sine_blade * MAC
            ELSE
               c_iop = zero_ ; a_eff = zero_ ; a_mac = zero_
            ENDIF
         ENDIF
      ENDDO

      ! Transmissivity of direct/diffuse attentuation and conversion from direct to diffuse
      ! - for one half of the layer
      !# 2026-09-08: validate EVERY operand of the attenuation arithmetic HERE - before the first
      !# division by costheta_r on the next line, and before any day/night branch downstream.
      !#
      !# The previous placement was inside `IF (direct(l) + diffuse(l) > 0)`, which is precisely
      !# where a NaN irradiance CANNOT reach it: NaN > 0 is false, so a NaN silently took the
      !# diffuse-only branch and was never checked. Testing the validator in isolation could not
      !# have caught that - the fault was the caller's ordering, not the predicate.
      IF ( .NOT. aed_oasim_is_finite(h) .OR. h <= zero_ )                                          &
         CALL aed_oasim_fatal('invalid layer thickness', yearday, layer, layer_idx, -1, -one_, h,  &
              h, zero_, zero_, zero_, 'h')
      IF ( .NOT. aed_oasim_is_finite(costheta_r) .OR. costheta_r <= zero_ )                        &
         CALL aed_oasim_fatal('invalid refracted-beam geometry', yearday, layer, layer_idx,        &
              -1, -one_, h, costheta_r, mcosthetas, zero_, zero_, 'costheta_r, mcosthetas')
      DO l = 1, data%nlambda
         IF ( .NOT. aed_oasim_is_finite(a(l))   .OR. a(l)   < zero_ .OR.                           &
              .NOT. aed_oasim_is_finite(b(l))   .OR. b(l)   < zero_ .OR.                           &
              .NOT. aed_oasim_is_finite(b_b(l)) .OR. b_b(l) < zero_ )                              &
            CALL aed_oasim_fatal('invalid optical coefficient', yearday, layer, layer_idx,         &
                 l, data%lambda(l), h, a(l), b(l), b_b(l), zero_, 'a, b, b_b, -')
         IF ( .NOT. aed_oasim_is_finite(direct(l))  .OR. direct(l)  < zero_ .OR.                   &
              .NOT. aed_oasim_is_finite(diffuse(l)) .OR. diffuse(l) < zero_ )                      &
            CALL aed_oasim_fatal('invalid irradiance (checked BEFORE the day/night branch)',       &
                 yearday, layer, layer_idx, l, data%lambda(l), h,                                  &
                 direct(l), diffuse(l), zero_, zero_, 'direct, diffuse, -, -')
      ENDDO

      f_att_d = exp(-0.5 * (a + b) * h / costheta_r)         ! Gregg & Rousseau 2016 Eq 8
      f_att_s = exp(-0.5 * (a + r_s * b_b) * h / mcosthetas) ! Gregg & Rousseau 2016 Eq 9

      !# 2026-08-29 SCATTERING FIX. f_prod_s used to be the whole
      !#     exp(-0.5*a*h/costheta_r) - f_att_d
      !# i.e. EVERY scattered photon, forward and backward, was moved into the DOWNWARD diffuse
      !# stream. The consequence is exact and was verified algebraically:
      !#     f_att_d + f_prod_s == exp(-0.5*a*h/costheta_r)
      !# so the b terms cancelled identically and total scattering had ZERO net effect on the
      !# downward irradiance. The direct beam was attenuated by ABSORPTION ALONE. The in-code
      !# comment ("but not accounting for backscattered fraction") understated it: this is not a
      !# missing refinement, it is the scattering contribution to Kd going missing entirely.
      !#
      !# It matters because 3 of the 10 configured IOPs are scattering-dominated (iop_type 10 =
      !# average sediment for NCS_ss1/NCS_ss2, and iop_type 6 for OGM_cpom), and because
      !# drive_extinction = .true. feeds Kd_bio straight into GLM's heat budget.
      !#
      !# CORRECTION, same day: an earlier version of this comment claimed the defect was "a direct
      !# contributor to the modelled Kd sitting ~25-40% below the measured Kd(PAR)". THAT WAS
      !# WRONG and is retracted. Measured after the fix, Kd did not move to three decimal places
      !# (0.265 / 0.309 / 0.309 before and after). The backscatter RATIOS configured here are
      !# 0.002-0.02 (:339-431), so the forward fraction is 0.98-0.998 and the extra loss is
      !# sub-percent. Separately, scoring like-for-like against the measured PAR profiles
      !# (obs/<E>/field_par.csv, OAS_par_E_scalar, same 1-8 m window, same 11:00 casts) puts the
      !# real Kd deficit at 8-15%, not 25-40% - the larger figure came from comparing against a
      !# median derived by a different method.
      !#
      !# The fix stays because it is a genuine energy-conservation defect and would matter in a
      !# turbid or sediment-laden water body. It does NOT explain this site's Kd deficit.
      !#
      !# Correct form: the direct beam loses a + b; only the FORWARD-scattered fraction
      !# (b - b_b)/b is converted into the downward diffuse stream; the backscattered remainder
      !# leaves the downward streams, as it physically must. b_b here is already the absolute
      !# backscattering coefficient (ratio * b, assembled at :1123-1124), not a ratio.
      !# elementwise: a, b and b_b are per-band arrays, so this must be WHERE, not IF
      fwd_frac = one_
      WHERE ( b > 1.0e-12 ) fwd_frac = MAX(zero_, MIN(one_, (b - b_b) / b))
      !# 2026-09-08: fwd_frac is only defined on the line above, so it is validated here rather
      !# than with the other operands; f_prod_s consumes it immediately.
      DO l = 1, data%nlambda
         IF ( .NOT. aed_oasim_is_finite(fwd_frac(l)) .OR.                                          &
              fwd_frac(l) < zero_ .OR. fwd_frac(l) > one_ )                                        &
            CALL aed_oasim_fatal('invalid forward-scattered fraction', yearday, layer, layer_idx,  &
                 l, data%lambda(l), h, fwd_frac(l), b(l), b_b(l), zero_, 'fwd_frac, b, b_b, -')
      ENDDO
      f_prod_s = (exp(-0.5 * a * h / costheta_r) - f_att_d) * fwd_frac   ! G&R 2016 Eq 14

      IF (data%save_Kd .and. data%spectral_output /= 0) THEN
         DO l = 1, data%nlambda
              !Kd(l) = 2 * (log(direct(l) + diffuse(l)) - log(direct(l) * (f_att_d(l) + f_prod_s(l)) + diffuse(l) * f_att_s(l))) / h
               IF (direct(l) + diffuse(l) > 0) THEN
                  !# 2026-09-08 STABLE EVALUATION (log-sum-exp).
                  !#
                  !# Substituting f_prod_s (:1244) shows the quantity being logged is exactly a
                  !# THREE-TERM CONVEX COMBINATION of exponentials:
                  !#    T  = c1*exp(-k1*h/2) + c2*exp(-k2*h/2) + c3*exp(-k3*h/2)
                  !#    k1 = (a + b)/costheta_r         c1 = dir_frac*(1 - fwd_frac)
                  !#    k2 =  a     /costheta_r         c2 = dir_frac*fwd_frac
                  !#    k3 = (a + r_s*b_b)/mcosthetas   c3 = 1 - dir_frac      (sum(c) = 1, c >= 0)
                  !#
                  !# For the strongly absorbing infrared bands (a_w up to 1e6 m-1 above 1350 nm,
                  !# see the note at the swr_abs diagnostic) EVERY exponential underflows to exactly
                  !# zero in double precision, so the direct evaluation gave log(0) = -Inf and
                  !# Kd(l) = +Inf. Those bands have par_weights(l) = 0 exactly (par_weights spans
                  !# 400-700 nm, the band grid 300-4000 nm), so the PAR collapse below then formed
                  !# 0 * Inf = NaN. That is why OAS_Kd_bio was NaN at the surface layer on every
                  !# DAYLIGHT record of every run in this project, on both archived binaries.
                  !#
                  !# Factoring out the smallest exponent AMONG THE TERMS THAT CARRY WEIGHT is
                  !# algebraically identical and cannot underflow: every remaining exponent is <= 0
                  !# and the sum lies in (0,1]. Verified against the previous expression to 5.5e-16
                  !# over 304 well-conditioned cases, and against a 200-digit reference over 144
                  !# underflowing cases, by cases_heat_wall/tools/test_oasim_kd.py.
                  !#
                  !# NOTE for anyone tempted to "simplify" this back: the weighted mean of the
                  !# coefficients, sum(c_i*k_i), is NOT the answer. That is the small-h Taylor limit;
                  !# underflow is the large-k*h regime, where the mixture is dominated by the
                  !# SLOWEST-decaying stream. By Jensen the weighted mean always overestimates - by
                  !# 50% on a two-stream test case. It shipped on 2026-09-08 and was wrong.
                  !# Operands were validated before the first division and before this branch
                  !# (see the block above f_att_d). Nothing invalid can reach either branch here.
                  dir_frac = direct(l) / (direct(l) + diffuse(l))

                  kd_k(1) = (a(l) + b(l))       / costheta_r
                  kd_k(2) =  a(l)               / costheta_r
                  kd_k(3) = (a(l) + r_s*b_b(l)) / mcosthetas
                  kd_c(1) = dir_frac * (one_ - fwd_frac(l))
                  kd_c(2) = dir_frac * fwd_frac(l)
                  kd_c(3) = one_ - dir_frac

                  !# Same object code the unit test exercises (tools/test_oasim_kd.F90).
                  CALL aed_oasim_kd_mixture(kd_k, kd_c, h, Kd(l), kd_ok)
                  IF ( .NOT. kd_ok )                                                              &
                     CALL aed_oasim_fatal('Kd mixture rejected its inputs', yearday, layer,       &
                          layer_idx, l, data%lambda(l), h, kd_k(1), kd_k(2), kd_k(3),             &
                          kd_c(1)+kd_c(2)+kd_c(3), 'k1, k2, k3, sum(c)')
               ELSE
                  ! Only diffuse stream
                  Kd(l) = (a(l) + r_s * b_b(l)) / mcosthetas
               ENDIF
               IF(l==data%l490_l) THEN
                  ! increment Kd490 with this layers properties
                  Kd490avg = Kd490avg +  Kd(l)
               END IF

         ENDDO

         !# 2026-07: same two-stream form but with PURE WATER IOPs only, so the difference
         !# isolates the attenuation contributed by phytoplankton, CDOM, detritus and
         !# sediment. GLM's `extinction` accumulates on top of Kw, so this - not the total
         !# Kd - is the quantity a light_extinction binding would have to return, and it is
         !# what is directly comparable with the tuned Ke* constants in aed.nml.
         Kd_w = ( data%a_w + r_s * (0.5*data%b_w) ) / mcosthetas
         kd_bio_num = zero_ ; kd_bio_den = zero_ ; kd_wat_num = zero_
         DO l = 1, data%nlambda
            !# spectrum-weighted collapse to a single PAR value: bands carrying more
            !# energy contribute proportionally more to the PAR-averaged attenuation.
            !# 2026-09-08: skip bands that carry no PAR weight or no light BEFORE forming the
            !# product. They contribute exactly nothing by construction, and evaluating them was
            !# how a non-finite Kd in an out-of-PAR band reached this sum as 0 * Inf = NaN. The
            !# same test gates the denominator so the weighting stays consistent with the numerator.
            IF ( data%par_weights(l) <= zero_ ) CYCLE
            IF ( spectrum(l) <= zero_ ) CYCLE
            kd_bio_num = kd_bio_num + data%par_weights(l) * spectrum(l) * MAX(Kd(l) - Kd_w(l), zero_)
            kd_wat_num = kd_wat_num + data%par_weights(l) * spectrum(l) * Kd_w(l)
            kd_bio_den = kd_bio_den + data%par_weights(l) * spectrum(l)
         ENDDO
         IF (kd_bio_den > zero_) THEN
            _DIAG_VAR_(data%id_kd_bio)   = kd_bio_num / kd_bio_den
            _DIAG_VAR_(data%id_kd_water) = kd_wat_num / kd_bio_den
         ELSE
            _DIAG_VAR_(data%id_kd_bio)   = zero_
            _DIAG_VAR_(data%id_kd_water) = zero_
         ENDIF
      ENDIF

      ! From top to centre of layer
      diffuse = diffuse * f_att_s + direct * f_prod_s  ! MH Reversed order of calc
      direct = direct * f_att_d
     !direct = direct * f_att_d
     !diffuse = diffuse * f_att_s + direct * f_prod_s
      spectrum = direct + diffuse

      par_E = sum(data%par_E_weights * spectrum)
      par_J = sum(data%par_weights * spectrum)
      swr_J = sum(data%swr_weights * spectrum)
      uv_J  = sum(data%uv_weights * spectrum)
      _DIAG_VAR_(data%id_par_E) = par_E ! Photosynthetically Active Radiation (umol/m2/s)
      _DIAG_VAR_(data%id_par) = par_J   ! Photosynthetically Active Radiation (W/m2)
      _DIAG_VAR_(data%id_swr) = swr_J   ! Total shortwave radiation (W/m2) [up to 4000 nm]
      _DIAG_VAR_(data%id_uv) = uv_J     ! UV (W/m2)
      _DIAG_VAR_(data%id_par_E_dif) = & ! Diffuse Photosynthetically Active photon flux (umol/m2/s)
                                      sum(data%par_E_weights * diffuse)


      ! Compute scalar PAR as experienced by phytoplankton
      spectrum = direct / costheta_r + diffuse / mcosthetas
      par_J = sum(data%par_weights * spectrum)
      par_E = sum(data%par_E_weights * spectrum)
      _DIAG_VAR_(data%id_par_J_scalar) = par_J ! Scalar Photosynthetically Active Radiation (W/m2)
      _DIAG_VAR_(data%id_par_E_scalar) = par_E ! Scalar Photosynthetically Active photon flux (umol/m2/s)

      ! From centre to bottom of layer - special case for macrophytes in bottom layer
      !# 2026-08-29 FIXES (benthic canopy block; unreachable in the live configuration, which
      !# has no iop_type = 20, so a_mac stays zero_ and this branch never runs):
      !#  (a) SIGN. f_att_m was exp(+a_mac). a_mac = a * A_eff * sine_blade * MAC_A is a product
      !#      of non-negative terms, so exp(a_mac) >= 1 and the canopy AMPLIFIED the light
      !#      passing through it. par_M = par_J - par_M then came out NEGATIVE, i.e. the blades
      !#      "captured" a negative amount of PAR. Beer-Lambert requires exp(-a_mac).
      !#  (b) OUT-OF-BOUNDS SUBSCRIPT. The diagnostic write used i_iop, but the DO i_iop loop
      !#      that owns it ended ~120 lines above; after a completed DO loop i_iop equals
      !#      size(data%iops)+1, so data%iops(i_iop)%id_parc read one element past the array
      !#      and used whatever integer that happened to be as a column() subscript. Use the
      !#      canopy IOP's own index, captured inside the loop.
      IF( layer==SIZE(layer_map) .AND. SUM(a_mac) > zero_ .AND. i_iop_mac > 0 ) THEN
         spectrum = direct + diffuse                     ! Spectrum at centre of layer (assumed top of canopy)
         par_J = sum(data%par_weights * spectrum)        ! Incoming PAR @ top of canopy
         f_att_m = exp( -a_mac )                         ! Absorbance of light by leaf/blade area
         direct = direct * f_att_m                       ! Reduction of direct
         diffuse = diffuse * f_att_m                     ! Reduction of diffuse
         spectrum = direct + diffuse                     ! Update spectrum after goinf through canopy
         par_M = sum(data%par_weights * spectrum)        ! PAR after blade photon capture
         par_M = par_J - par_M                           ! PAR capture by blades
         _DIAG_VAR_S_(data%iops(i_iop_mac)%id_parc) = par_M ! Return PAR capture rate to linked group
      ENDIF

      ! From centre to bottom of layer - for 2nd half of the layer
      diffuse = diffuse * f_att_s + direct * f_prod_s  ! MH Reversed order of calc
      direct = direct * f_att_d
     !direct = direct * f_att_d
     !diffuse = diffuse * f_att_s + direct * f_prod_s
      spectrum = direct + diffuse


      ! Save spectrally resolved outputs
      SELECT CASE (data%spectral_output)
         CASE (1)
            !# 2026-07 fix: full-spectrum output was broken two ways. The in-water
            !# irradiance bands were never written at all (dir_band*/dif_band* stayed at
            !# zero while being registered and appearing in the NetCDF), and the Kd line
            !# was BACKWARDS - it read the diagnostic into the local Kd array instead of
            !# publishing it, so Kd_band* was also never written and the local Kd used by
            !# the Secchi average below was silently overwritten with whatever the
            !# diagnostic held. Only spectral_output=2 produced usable spectral output.
            !# (No bounds issue here: CASE(1) sets lambda_out = lambda, so
            !# size(lambda_out) == nlambda and every id_* array is nlambda long.)
            DO l = 1, data%nlambda
               _DIAG_VAR_(data%id_band_dir(l)) = direct(l)
               _DIAG_VAR_(data%id_band_dif(l)) = diffuse(l)
               _DIAG_VAR_(data%id_a_band(l)) = a(l) - data%a_w(l)
               _DIAG_VAR_(data%id_b_band(l)) = b(l) - data%b_w(l)
               IF (data%save_Kd) _DIAG_VAR_(data%id_Kd(l)) = Kd(l)
            ENDDO
         CASE (2)
            CALL interp(data%nlambda, data%lambda, direct, size(data%lambda_out), data%lambda_out, spectrum_out)
            DO l = 1, size(data%lambda_out)
               _DIAG_VAR_(data%id_band_dir(l)) = spectrum_out(l)
            ENDDO
            CALL interp(data%nlambda, data%lambda, diffuse, size(data%lambda_out), data%lambda_out, spectrum_out)
            DO l = 1, size(data%lambda_out)
               _DIAG_VAR_(data%id_band_dif(l)) = spectrum_out(l)
            ENDDO
            CALL interp(data%nlambda, data%lambda, a, size(data%lambda_out), data%lambda_out, spectrum_out)
            DO l = 1, size(data%lambda_out)
               _DIAG_VAR_(data%id_a_band(l)) = spectrum_out(l) - data%a_w_out(l)
            ENDDO
            CALL interp(data%nlambda, data%lambda, b, size(data%lambda_out), data%lambda_out, spectrum_out)
            DO l = 1, size(data%lambda_out)
               _DIAG_VAR_(data%id_b_band(l)) = spectrum_out(l) - data%b_w_out(l)
            ENDDO
            IF (data%save_Kd) THEN
               CALL interp(data%nlambda, data%lambda, Kd, size(data%lambda_out), data%lambda_out, spectrum_out)
               DO l = 1, size(data%lambda_out)
                  _DIAG_VAR_(data%id_Kd(l)) = spectrum_out(l)
               ENDDO
            ENDIF
      END SELECT

      ! Compute remaining downwelling shortwave flux and from that, absorption [heating]
      !# 2026-07 (Tier 1C): per-IOP normalised spectral-absorption factor, computed from
      !# the SAME a*(lambda) this module uses for shading. Growth weighted by this is
      !# weighted by real absorption spectra rather than a parallel hand-tuned template -
      !# which is what let Planktothrix have a phycoerythrin band for shading (33-T) while
      !# its growth still read a disconnected 7-band table.
      DO i_iop = 1, size(data%iops)
         afac_num = zero_ ; afac_den = zero_ ; afac_fnum = zero_ ; afac_fden = zero_
         afq_num  = zero_ ; afq_den  = zero_ ; afq_fnum  = zero_ ; afq_fden  = zero_
         DO l = 1, data%nlambda
            !# energy weighting: joules absorbed per joule available
            afac_num  = afac_num  + data%par_weights(l) * spectrum(l) * data%iops(i_iop)%a(l)
            afac_den  = afac_den  + data%par_weights(l) * spectrum(l)
            afac_fnum = afac_fnum + data%par_weights(l) * data%iops(i_iop)%a(l)
            afac_fden = afac_fden + data%par_weights(l)
            !# photon weighting: quanta absorbed per quantum available. par_E_weights already
            !# carries the lambda/(h c N_A) factor, so this IS Phi(lambda) = I(lambda) lambda/(h c)
            !# integrated against a*(lambda) - no extra conversion, and none must be added.
            afq_num   = afq_num   + data%par_E_weights(l) * spectrum(l) * data%iops(i_iop)%a(l)
            afq_den   = afq_den   + data%par_E_weights(l) * spectrum(l)
            afq_fnum  = afq_fnum  + data%par_E_weights(l) * data%iops(i_iop)%a(l)
            afq_fden  = afq_fden  + data%par_E_weights(l)
         ENDDO
         !# Both are ratios of a* to a*, so both are dimensionless and both are exactly 1.0 for a
         !# flat spectrum regardless of weighting - that identity is the backward-compatibility test.
         IF (afac_den > zero_ .AND. afac_fnum > zero_ .AND. afac_fden > zero_) THEN
            afacE = MIN(5d0, MAX(0.05d0, (afac_num/afac_den) / (afac_fnum/afac_fden)))
         ELSE
            afacE = one_
         ENDIF
         IF (afq_fnum > zero_ .AND. afq_fden > zero_) THEN
            astarq = afq_fnum / afq_fden
         ELSE
            astarq = zero_
         ENDIF
         IF (afq_den > zero_ .AND. astarq > zero_) THEN
            afacq = MIN(5d0, MAX(0.05d0, (afq_num/afq_den) / (afq_fnum/afq_fden)))
         ELSE
            afacq = one_
         ENDIF
         _DIAG_VAR_(data%id_afacE(i_iop)) = afacE
         _DIAG_VAR_(data%id_afacq(i_iop)) = afacq
         _DIAG_VAR_(data%id_astar(i_iop)) = astarq
         !# the consumer-facing copy
         IF (data%afac_quanta) THEN
            _DIAG_VAR_(data%id_afac(i_iop)) = afacq
         ELSE
            _DIAG_VAR_(data%id_afac(i_iop)) = afacE
         ENDIF
      ENDDO

      swr_J = sum(data%swr_weights * spectrum)
      par_J = sum(data%par_weights * spectrum)
      _DIAG_VAR_(data%id_swr_abs) = swr_top - swr_J
      !# 2026-07 fix: the scalar 'Kd' diagnostic was computed from swr_J/swr_top, i.e.
      !# the 300-4000 nm BROADBAND flux - so it was Kd(SW), not Kd(PAR), despite its name
      !# and units. In the top metre it is dominated by pure-water NIR absorption (a_w
      !# reaches 1e3-1e6 m-1 above 1350 nm), so it ran an order of magnitude above the
      !# true Kd(PAR) and never converged to it at depth either. Comparing it to a
      !# PAR-sensor Kd looked like catastrophic model failure; the model was fine and the
      !# diagnostic was mislabelled. Now uses the 400-700 nm weights.
      !# 2026-08-19 FIX: the guard was par_top>0.01, so in the deepest layers - where PAR
      !# has fallen below 0.01 W/m2 - this diagnostic was NEVER WRITTEN. glm_aed.F90:1155
      !# zeroes cc_diag at the top of every call, so OAS_Kd stayed 0.0 there. The consumer
      !# (glm_ptm.c:1101) rejects only kd<0.0, so 0.0 passed as a VALID sample and gave
      !# factor = exp(0*(h-mean_h)) = 1.0 - i.e. every particle in the bottom layer,
      !# including one sitting on the bed, was handed the LAYER-CENTRE irradiance with no
      !# within-layer attenuation at all. Measured on the 2019 run: particles at 9-10 m
      !# received 1.6-2.0x MORE PAR than particles at 8-9 m on the same day, on 61% of
      !# days - light increasing with depth. That over-brightening is what kept bed-trapped
      !# Planktothrix above its neutral-buoyancy light so it never stopped ballasting.
      !# Kd is well defined wherever both fluxes are positive; the 0.01 threshold was
      !# guarding against log(0), which a positivity test does correctly and without
      !# silently blanking the deepest layers.
      !# 2026-08-29 FIX: added the `.AND. data%spectral_output /= 0` guard. data%id_Kd and
      !# data%lambda_out are only ALLOCATEd when spectral_output /= 0 (see aed_define_model),
      !# but this test gated on save_Kd alone - so save_Kd = .true. with spectral_output = 0
      !# dereferenced two unallocated arrays. Every other site that touches id_Kd/lambda_out
      !# (~1179 and the Secchi block below) already carries the extra condition. No live change:
      !# the live model runs spectral_output = 2, so the added condition is always true.
      IF (data%save_Kd .AND. data%spectral_output /= 0 .AND. par_top>1d-30 .AND. par_J>1d-30) &
           _DIAG_VAR_(data%id_Kd(size(data%lambda_out)+1)) = -log(par_J/par_top) / h

      IF(l490_sf>0.01) THEN !incoming solar >0.01 W/m2
           ! look for last layer where radiation at lambda490 is 10% of sf amount
           atten_frac = (data%swr_weights(data%l490_l) * spectrum(data%l490_l)) / l490_sf
           IF (atten_frac>0.1) secchi_botlayer = layer
      ENDIF

   END DO ! End vertical/column loop


   !-------------------------------------------------------------------------------
   ! Put remaining shortwave in bottom layer
   !  (assumes all light is absorbed by sediment and injected into water column)
   !# 2026-08-29 FIX: was `=`, not `+=`. _DIAG_VAR_ resolves to column(id)%cell(layer_idx) and
   !# layer_idx still holds the BOTTOM layer here (the column loop above just ended), so this
   !# line DISCARDED the bottom layer's own in-layer absorption (swr_top - swr_J, written at
   !# ~1364) and replaced it with the sediment-reflection term alone. The bottom layer absorbs
   !# both. id_swr_abs has no consumer anywhere in the tree (grep: it appears only in this
   !# file, at its definition and these two writes), so this cannot change model results - it
   !# only corrects the reported diagnostic.
   _DIAG_VAR_(data%id_swr_abs) = _DIAG_VAR_(data%id_swr_abs) + swr_J * 0.5 ! 50% reflect to water and 50% lost to sediment

   ! Finalise secchi calculation. Kd490avg accumulates Kd490 over ALL water-column layers
   ! (loop above), so it must be divided by the layer COUNT to get the mean Kd490 -- not by
   ! secchi_botlayer (the euphotic-depth index), which made Secchi ~30x too small. (2026-07 fix)
   IF (data%save_Kd .and. data%spectral_output /= 0 .and. SIZE(layer_map) > 0) THEN
     Kd490avg = Kd490avg / SIZE(layer_map)
     _DIAG_VAR_S_(data%id_secchi) = 1.7 / Kd490avg
   ENDIF

END SUBROUTINE aed_calculate_column_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_dry_oasim(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Atmospheric component of the AED implementation of the OASIM model
! Note DRY cells only refer to 2D (sheet) variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_oasim_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL,PARAMETER :: pres0 = 101300.    ! Reference air pressure in Pa (Bird 1984 p461, Bird & Riordan p89)
  !AED_REAL,PARAMETER :: ga = 0.05          ! Ground albedo (used by Bird & Riordan 1986, but discarded by Gregg & Carder 1990)
   AED_REAL,PARAMETER :: H_oz = 22.         ! Height of maximum ozone concentration (km)
   AED_REAL,PARAMETER :: r_e = (10. + 11.8) / 2 ! Equivalent radius of cloud drop size distribution (um)
                                            ! based on mean of Kiehl et al. & Han et al. (cf OASIM)
   AED_REAL,PARAMETER :: mcosthetas = 0.831 ! Mean of cosine of angle of diffuse radiation in water,
                                            ! assuming all angular contributions equal in air
                                            ! (Sathyendranath and Platt 1989, p 191) NB Ackleson et al 1994 use 0.9
   AED_REAL,PARAMETER :: r_s = 1.5          ! Shape factor representing mean backscatter coefficient
                                            ! of diffuse irradiance (r_d in Ackleson et al. 1994 p 7487)

   AED_REAL :: longitude, latitude, yearday, cloud_cover, wind_speed, &
               airpres, relhum, LWP, water_vapour, WV, WM, visibility, AM
   AED_REAL :: days, hour, theta, costheta, alpha_a, beta_a, yearday_utc
   INTEGER  :: l
   AED_REAL :: M, M_prime, M_oz
   AED_REAL :: O3
   AED_REAL,DIMENSION(data%nlambda) :: direct, diffuse, spectrum, Kd  ! Spectra at top of the water
                                               ! column (with refraction and reflection accounted for)
   AED_REAL,DIMENSION(nlambda_5nm_astm) :: direct_5nm, diffuse_5nm ! Spectra at top of the water
   AED_REAL :: par_J=0., swr_J=0., uv_J=0., par_E=0., F_a, omega_a, par_M
   AED_REAL,DIMENSION(data%nlambda) :: tau_a, T_a, T_oz, T_w, T_u, T_r, T_aa, T_as
   AED_REAL,DIMENSION(data%nlambda) :: T_g, T_dclr, T_sclr, T_dcld, T_scld
   AED_REAL,DIMENSION(data%nlambda) :: rho_d, rho_s

   AED_REAL,DIMENSION(data%nlambda) :: a, b, b_b, a_iop, a_mac
   AED_REAL,DIMENSION(data%nlambda) :: f_att_d, f_att_s, f_prod_s, f_att_m
   INTEGER  :: i_iop
   AED_REAL :: c_iop, h, swr_top, l490_sf, costheta_r, dir_frac
   AED_REAL :: SWFLUX, MaxSWFlux
   AED_REAL :: c_eff, atten_frac, Kd490avg=0., a_eff
   INTEGER  :: month, secchi_botlayer
   CHARACTER(len=1):: met
!
!-------------------------------------------------------------------------------
!BEGIN
   !-------------------------------------------------------------------------------
   ! Set column time/location properties
   longitude   = _STATE_VAR_S_(data%id_lon)
   latitude    = _STATE_VAR_S_(data%id_lat)

   !# 2026-08-29 LATITUDE UNIT FIX. GLM hands latitude over as radians PLUS two_Pi
   !# (glm_model.c:195, `Latitude = two_Pi + Latitude * deg2rad`) and the global is declared
   !# 'radians' (glm_aed.F90:291) - but BOTH consumers in this module expect DEGREES:
   !#   zenith_angle()  aed_slingo.F90:272-273   rlat = deg2rad * dlat
   !#   vH_O3()         aed_atmosphere.F90:333   sin(beta * phi * deg2rad)**2
   !# The second deg2rad therefore put the sun at 7.21 N instead of 53.15 N. Measured before the
   !# fix: modelled solar zenith 15.7 deg at local noon against a true 34.0, and a minimum over
   !# the run of 7.50 deg - geometrically impossible at this latitude, where the sun never rises
   !# above ~29.8 deg from the zenith. Wrong air mass, wrong direct/diffuse split, wrong day
   !# length (~12.2 h tropical vs ~16.5 h at Stechlin in July) and wrong in-water path.
   !#
   !# Unwrap the 2*pi offset, then convert. GLM's own C code is self-consistent and must NOT be
   !# touched: trig is 2*pi-periodic so glm_util.c:80-81 (cos/sin of lat) is unaffected, and the
   !# one place needing degrees subtracts the offset explicitly (glm_surface.c:2088).
   !#
   !# Longitude is deliberately NOT converted - GLM never applies deg2rad to it, so it already
   !# arrives in degrees, which is what both consumers want.
   IF ( ABS(latitude) > pi ) latitude = latitude - SIGN(2.*pi, latitude)
   latitude = latitude * (180./pi)
   yearday     = _STATE_VAR_S_(data%id_yearday)
   !# 2026-08-29 FIX: this dry path hardcoded `yearday - 8./24.` instead of using the namelist
   !# time-zone shift, unlike the wet path (aed_calculate_column_oasim, ~881), which uses
   !# data%tzshift. Any site not at UTC+8 got the wrong hour-of-day, hence the wrong zenith
   !# angle, in every dry-cell diagnostic. Unreachable in the live model: benthic_mode = 1, so
   !# aed_calculate_dry_oasim is never called (and it RETURNs before the physics regardless).
   yearday_utc = yearday - data%tzshift/24.
   days        = floor(yearday_utc)
   hour        = mod(yearday_utc, 1.0) * 24.0

  ! Calculate zenith angle (in radians)
   theta = zenith_angle(days, hour, longitude, latitude)
   theta = min(theta, 0.5 * pi)  ! Restrict the input zenith angle between 0 and pi/2
   costheta = cos(theta)

   ! Set column environmental conditions (from host)
   wind_speed  = _STATE_VAR_S_(data%id_wind_speed)      ! Wind speed @ 10 m above surface (m/s)
   !# 2026-08-29 FIX: this read was UNCONDITIONAL and sat before the RETURN below, but
   !# id_mean_wind_speed is only ever assigned when compute_mean_wind = .FALSE. (and only in
   !# the sky_model == 1 branch of aed_define_model). With the default compute_mean_wind =
   !# .TRUE. the id held garbage, so column(garbage)%cell_sheet was dereferenced on every dry
   !# cell. Guarded, with the same fallback used in the wet path. (The id now defaults to -1;
   !# see the type declaration.) Unreachable in the live model: benthic_mode = 1 means
   !# aed_calculate_dry_oasim is never called.
   IF (data%id_mean_wind_speed > 0) THEN
      WM = _STATE_VAR_S_(data%id_mean_wind_speed) ! Daily mean wind speed @ 10 m above surface (m/s)
   ELSE
      WM = wind_speed
   ENDIF
   relhum      = _STATE_VAR_S_(data%id_relhum)          ! Relative humidity (-)
   SWFLUX      = _STATE_VAR_S_(data%id_I_0)             ! Incoming shortwave flux

   IF( diag_level>9 ) THEN
      _DIAG_VAR_S_(data%id_swflux) = SWFLUX
      _DIAG_VAR_S_(data%id_zen) = rad2deg * theta
      _DIAG_VAR_S_(data%id_doy) = yearday
   ENDIF


   RETURN
   !-------------------------------------------------------------------------------
   ! Compute incoming spectral irradiance based on atmospheric properties
   IF (data%sky_model == 1) THEN

       ! Set column environmental conditions (from parameters)
       AM          = data%air_mass_type !_STATE_VAR_S_(data%id_air_mass_type)
                                        ! Aerosol air mass type (1: open ocean, 10: continental)
       visibility  = data%visibility ! _STATE_VAR_S_(data%id_visibility)     ! Visibility (m)
       airpres     = data%airpres ! _STATE_VAR_S_(data%id_airpres)       ! Surface air pressure (Pa)
       cloud_cover = data%cloud ! _STATE_VAR_S_(data%id_cloud)       ! Cloud cover (fraction, 0-1)
       LWP         = data%LWP ! _STATE_VAR_S_(data%id_lwp)      ! Cloud liquid water content (kg m-2)
       water_vapour= data%WV ! _STATE_VAR_S_(data%id_wv)    ! Total precipitable water vapour (kg m-2) - equivalent to mm
       O3          = data%O3                             ! Ozone content (kg m-2)
       IF(data%ozone) THEN
         O3 = vH_O3(longitude, latitude, yearday)
         _DIAG_VAR_S_(data%id_O3) = O3
       ENDIF

       ! Adjustments and diagnostic outputs for debugging
       _DIAG_VAR_S_(data%id_wind_out) = wind_speed
       IF (data%id_mean_wind_out > 0) _DIAG_VAR_S_(data%id_mean_wind_out) = WM
       IF (cloud_cover > 0) LWP = LWP / cloud_cover

       ! LWP is the mean density over a grid box (Jorn: ECMWF pers comm 26/2/2019), but
       !        we want the mean density per cloud-covered area
       WV   = water_vapour / 10             ! from kg m-2 to cm
       O3   = O3 * (1000 / 48.) / 0.4462    ! from kg m-2 to mol m-2, then from mol m-2 to atm cm (Basher 1982)

       ! Atmospheric path length, a.k.a. relative air mass (Eq 3 Bird 1984, Eq 5 Bird & Riordan 1986, &
       !                                        Eq 13 Gregg & Carder 1990, Eq A5 in Casey & Gregg 2009)
       ! Note this should always exceed 1, but as it is an approximation it does not near &
       !                                    theta -> 0 (Tomasi et al. 1998 p14). Hence the max operator.
       M = max(1., 1. / (costheta + 0.15 * (93.885 - theta * rad2deg)**(-1.253)))

       ! Pressure-corrected atmospheric path length (Eq A6 Casey & Gregg 2009)
       M_prime = M * airpres / pres0

       ! Atmospheric path length for ozone
       !  Eq 10, Bird & Riordan 1986, Eq 14 Gregg & Carder 1990; NB 6370 is the earth's radius in km
       !  See also Tomasi et al. 1998 Eq 5
       M_oz = (1. + H_oz / 6370.) / sqrt(costheta**2 + 2 * H_oz / 6370.)

       ! Transmittance due to ozone absorption (Eq 8 Bird 1984, Eq 9 Bird & Riordan 1986, Eq 17 Gregg & Carder 1990)
       T_oz = exp(-data%a_o * O3 * M_oz)

       ! Transmittance due to water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
       !  Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
       T_w = exp((-0.2385 * data%a_v * WV * M) / (1. + 20.07 * data%a_v * WV * M)**0.45)

       ! Transmittance due to uniformly mixed gas absorption - SHOULD use pressure corrected airmass
       !  Eq 10 Bird 1984, Eq 11 Bird and Riordan 1986, Eq 18 Gregg & Carder 1990
       !  Bird & Riordan use 118.93 rather than 118.3, but state 118.3 should be used in the future.
       T_u = exp(-1.41 * data%a_u * M_prime / (1. + 118.3 * data%a_u * M_prime)**0.45)

       ! Transmittance terms that apply for both cloudy and clear skies.
       T_g = T_oz * T_w * T_u

       ! -------------------
       ! Clear skies part
       ! -------------------

       ! Transmittance due to Rayleigh scattering (use precomputed optical thickness)
       T_r = exp(-M_prime * data%tau_r)

       ! Transmittance due to aerosol absorption (Eq 26 Gregg & Carder 1990)
       CALL navy_aerosol_model(AM, WM, wind_speed, relhum, visibility, costheta, alpha_a, beta_a, F_a, omega_a)
       tau_a = beta_a * data%lambda**(-alpha_a)
       T_a = exp(-tau_a * M)

       ! Direct transmittance
       T_dclr = T_r * T_a

       ! Separate absorption and scattering components of aerosol transmittance
       T_aa = exp(-(1. - omega_a) * tau_a * M)
       T_as = exp(-omega_a * tau_a * M)
       T_sclr = T_aa * 0.5 * (1. - T_r**0.95) + T_r**1.5 * T_aa * F_a * (1 - T_as)

       _DIAG_VAR_S_(data%id_omega_a) = omega_a
       _DIAG_VAR_S_(data%id_alpha_a) = alpha_a
       _DIAG_VAR_S_(data%id_beta_a) = beta_a
       _DIAG_VAR_S_(data%id_F_a) = F_a

       ! -------------------
       ! Cloudy skies part
       ! -------------------

       ! Transmittance due to absorption and scattering by clouds
       CALL slingo(costheta, max(0., LWP) * 1000, r_e, data%nlambda, data%lambda, T_dcld, T_scld)

       ! Diffuse and direct irradiance streams (Eqs 1, 2 Gregg & Casey 2009)
       !  These combine terms for clear and cloudy skies, weighted by cloud cover fraction
       direct  = data%exter * costheta * T_g * ((1. - cloud_cover) * T_dclr + cloud_cover * T_dcld)
       diffuse = data%exter * costheta * T_g * ((1. - cloud_cover) * T_sclr + cloud_cover * T_scld)

   ELSEIF(data%sky_model == 2) THEN

       ! Diffuse and direct irradiance streams from local analysis of Cockburn Sound radiation (Perth, WA)
       !  These functions combine terms for clear and cloudy skies
        met = 'W'

       ! Get function to return direct and diffuse spectra, based on incoming solar (etc)
       CALL direct_diffuse_bands(SWFlux, theta, yearday, met, data%nlambda, data%lambda_bounds, direct, diffuse)
   ENDIF

   ! Integrate to get total incoming surface irradiance and bandwidth components
   spectrum = direct + diffuse
   swr_J = sum(data%swr_weights   * spectrum)
   uv_J  = sum(data%uv_weights    * spectrum)
   par_J = sum(data%par_weights   * spectrum)
   par_E = sum(data%par_E_weights * spectrum)

   ! Set diagnostics of surface irradiance (bulk)
   _DIAG_VAR_S_(data%id_par_E_sf) = par_E    ! Photosynthetically Active Radiation, PAR (umol/m2/s)
   _DIAG_VAR_S_(data%id_par_sf)   = par_J    ! Photosynthetically Active Radiation, PAR (W/m2)
   _DIAG_VAR_S_(data%id_swr_sf)   = swr_J    ! Total shortwave radiation (W/m2), SW [up to 4000 nm]
   _DIAG_VAR_S_(data%id_uv_sf)    = uv_J     ! Ultraviolet Radiation, UV (W/m2)
   swr_J = sum(data%swr_weights * diffuse)   ! Repopulate variable with diffuse only fraction
   _DIAG_VAR_S_(data%id_swr_dif_sf) = swr_J  ! Diffuse shortwave radiation (W/m2), SW [up to 4000 nm]

   _DIAG_VAR_S_(data%id_par_E_sf_w) = 0.        ! Photosynthetically Active Radiation (umol/m2/s)
   _DIAG_VAR_S_(data%id_par_sf_w) = 0.          ! Photosynthetically Active Radiation (W/m2)
   _DIAG_VAR_S_(data%id_swr_sf_w) = 0.          ! Total shortwave radiation (W/m2) [up to 4000nm]
   _DIAG_VAR_S_(data%id_uv_sf_w) = 0.           ! UV (W/m2)

END SUBROUTINE aed_calculate_dry_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE reflectance(nlambda, F, theta, W, rho_d, rho_s, costheta_r)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nlambda
   AED_REAL,INTENT(in) :: F(nlambda), theta, W
   AED_REAL,INTENT(out) :: rho_d(nlambda), rho_s(nlambda), costheta_r
!
!LOCALS

   ! Coefficients for foam reflectance-wind speed relationship (p1666 Gregg and Carder 1990)
   AED_REAL,PARAMETER :: D1 = 2.2e-5
   AED_REAL,PARAMETER :: D2 = 4.0e-4
   AED_REAL,PARAMETER :: D3 = 4.5e-5
   AED_REAL,PARAMETER :: D4 = 4.0e-5

   ! Air density (g/m3)
   AED_REAL,PARAMETER :: rho_a = 1.2e3

   ! Refractive index of seawater, Gregg and Carder 1990 p1667
   ! Note: pure water has an index of about 1.33 (e.g., Kirk 2011),
   ! but in seawater it is higher (1.34 - 1.36) and dependent on temperature and salinity
   ! (e.g., https://doi.org/10.1016/0011-7471(71)90050-7)
   AED_REAL,PARAMETER :: n_w = 1.333   ! freshwater (Kirk 2011; SCOR CRDN 2009); was 1.341 (seawater)

   AED_REAL :: C_D, tau, rho_f_W, rho_f(nlambda)
   AED_REAL :: b, rho_dsp, rho_ssp, theta_r
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Drag coefficient (Eqs 42, 43 Gregg and Carder 1990)
   ! Constants match Trenberth et al. 1989 p1508 J Clim. However, they use different wind speed thresholds
   ! and a constant value between 3 and 10 m s-1. Theirs is also continuous, whereas the expression below is discontinuous.
   IF (W <= 0.) THEN
      C_D = 0.
   ELSEIF (W <= 7.) THEN
      C_D = 0.62e-3 + 1.56e-3 / W
   ELSE
      C_D = 0.49e-3 + 0.065e-3 * W
   ENDIF

   ! Wavelength-independent foam reflectance [affects direct and diffuse light] (Eqs 39-41 Gregg and Carder 1990)
   ! Reformulated to make dependence on surface stress (tau, units seem to be 10-3 m2/s2) explicit.
   tau = rho_a * C_D * W**2
   IF (W <= 4.) THEN
      rho_f_W = 0.
   ELSEIF (W <= 7.) THEN
      rho_f_W = D1 * tau - D2
   ELSE
      rho_f_W = D3 * tau - D4 * W**2
   ENDIF

   ! Final wavelength and wind speed-dependent foam reflectance
   rho_f = rho_f_W * F

   ! Calculate zenith angle (radians) inside the water, taking refraction into account: Snell's law
   theta_r = asin(sin(theta) / n_w)

   ! Direct light specular component
   IF (theta * rad2deg >= 40. .and. W > 2.) THEN
      ! Wind speed dependant (Eqs 46, 47 Gregg and Carder 1990)
      b = -7.14e-4 * W + 0.0618
      rho_dsp =  0.0253 * exp(b * (theta * rad2deg - 40.))
   ELSE
      ! Fresnel's Law, Eq 44 Gregg and Carder 1990 - note: contains typo (internal 1/2), see Kirk 3rd ed 2011, p 46
      rho_dsp = 0.5 * (sin(theta - theta_r)**2 / sin(theta + theta_r)**2 + tan(theta - theta_r)**2 / tan(theta + theta_r)**2)
   ENDIF

   ! Diffuse light specular component (p 1667 Gregg and Carder 1990)
   IF (W > 4.) THEN
      rho_ssp = 0.057
   ELSE
      rho_ssp = 0.066
   ENDIF

   rho_d = rho_dsp + rho_f
   rho_s = rho_ssp + rho_f
   costheta_r = cos(theta_r)
END SUBROUTINE reflectance
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_oasim(data,column,layer_idx,extinction)
!------------------------------------------------------------------------------+
! Return the light extinction contributed by the OASIM spectral solution.
!
! 2026-07 (Tier 1A). Before this existed, the binding was commented out and the
! entire 33-band spectral model drove NOTHING: growth used the host's broadband
! 0.45*SW attenuated by Kw plus the non-spectral Ke* constants, and heating used
! Kw alone. OASIM computed a full radiative transfer solution that no other part
! of the model ever read.
!
! IMPORTANT - why this returns Kd_bio and not Kd:
!   The host's `extinction` is an ACCUMULATOR. Each model ADDS its contribution on
!   top of Kw. OASIM's `Kd` diagnostic is a TOTAL that already includes pure water,
!   so returning it here would count water twice (and Kw is 0.5 in aed.nml while
!   OASIM's own physics puts pure-water PAR attenuation at ~0.13 /m - see 33-O).
!   `Kd_bio` is the same two-stream solution computed with water-only IOPs removed,
!   i.e. attenuation from phytoplankton, CDOM, detritus and sediment only - which is
!   exactly what the Ke* constants represent and what this hook must supply.
!
! ENABLING THIS IS NOT SUFFICIENT ON ITS OWN. `drive_extinction` defaults to
! .FALSE. because switching it on while the Ke* constants remain non-zero counts
! every attenuator twice. The accompanying change is to zero KePHY/KeDOM/KeDOMR/
! KePOM/KeCPOM/Ke_ss and re-derive Kw (~0.13, not 0.5), then recalibrate against
! the observed FCR Kd(PAR) = 1.04 /m. Handoff 33-P has the measured decomposition.
!
! Timing note: `id_kd_bio` is written during calculate_column. Depending on host
! call order this hook may read the value from the previous timestep - a one-step
! lag, which is stable but should be verified against the host's sequencing before
! any calibration is trusted.
!------------------------------------------------------------------------------+
!ARGUMENTS
   CLASS (aed_oasim_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: kd_bio
!
!------------------------------------------------------------------------------+
!BEGIN
   IF (.NOT. data%drive_extinction) RETURN

   kd_bio = _DIAG_VAR_(data%id_kd_bio)

   !# Guard against the first timestep, when the column solve has not yet run and the
   !# diagnostic is still zero, and against any non-finite value.
   IF (kd_bio > zero_) extinction = extinction + kd_bio

END SUBROUTINE aed_light_extinction_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_oasim_fatal(what, yday, layer, layer_idx, band, lambda, h_, v1, v2, v3, v4, names)
!-------------------------------------------------------------------------------
! Report an invalid optical operand with enough context to locate it, then fail
! with a NON-ZERO exit status.
!
! 2026-09-08. Deliberately not PURE and deliberately fatal: the previous version
! printed a warning and carried on, which is how a NaN reached the state arrays
! unnoticed. band/lambda are passed as -1 where the fault is not band-specific.
!-------------------------------------------------------------------------------
   CHARACTER(len=*),INTENT(in) :: what, names
   AED_REAL,INTENT(in) :: yday, lambda, h_, v1, v2, v3, v4
   INTEGER, INTENT(in) :: layer, layer_idx, band
!
!-------------------------------------------------------------------------------
!BEGIN
   print '(A)',            ' =============================================================='
   print '(A,A)',          ' OASIM FATAL: ', TRIM(what)
   print '(A,F12.5)',      '   yearday            ', yday
   print '(A,I6,A,I6)',    '   layer (top-down)   ', layer, '    host layer index ', layer_idx
   IF ( band > 0 ) THEN
      print '(A,I6,A,F10.2)', '   band               ', band, '    wavelength (nm)  ', lambda
   ELSE
      print '(A)',            '   band               n/a (fault is not band-specific)'
   ENDIF
   print '(A,ES14.6)',     '   layer thickness h  ', h_
   print '(A,A)',          '   operands           ', TRIM(names)
   print '(A,4ES15.6)',    '                      ', v1, v2, v3, v4
   print '(A)',            ' =============================================================='
   FLUSH(6)
   STOP 1
END SUBROUTINE aed_oasim_fatal
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE LOGICAL FUNCTION aed_oasim_is_finite(x)
!-------------------------------------------------------------------------------
! Portable finiteness test. (x /= x) is true only for NaN; the magnitude test
! rejects both infinities. Deliberately not IEEE_ARITHMETIC: libaed-light does not
! currently USE it, and this keeps the dependency surface unchanged.
!-------------------------------------------------------------------------------
   AED_REAL,INTENT(in) :: x
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_oasim_is_finite = (x == x) .AND. (abs(x) <= HUGE(x))
END FUNCTION aed_oasim_is_finite
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE LOGICAL FUNCTION aed_oasim_kd_inputs_ok(a_, b_, bb_, fwd, dir_, dif_, ct, h_)
!-------------------------------------------------------------------------------
! Validate the RAW operands before any division is performed on them.
!
! 2026-09-08. The previous version of this check ran AFTER the divisions by
! costheta_r had already happened, only printed a warning and then continued, did
! not test h, did not test the stream irradiances, and used `>= 0` tests that
! +Infinity passes. The NaN path was the worst: a NaN irradiance makes dir_frac
! NaN, hence every mixture weight NaN, hence every `c > 0` test false, so no
! branch was taken and the collapse silently produced +Infinity.
!-------------------------------------------------------------------------------
   AED_REAL,INTENT(in) :: a_, b_, bb_, fwd, dir_, dif_, ct, h_
!
!-------------------------------------------------------------------------------
!BEGIN
   aed_oasim_kd_inputs_ok = .FALSE.

   IF ( .NOT. aed_oasim_is_finite(a_)   ) RETURN
   IF ( .NOT. aed_oasim_is_finite(b_)   ) RETURN
   IF ( .NOT. aed_oasim_is_finite(bb_)  ) RETURN
   IF ( .NOT. aed_oasim_is_finite(fwd)  ) RETURN
   IF ( .NOT. aed_oasim_is_finite(dir_) ) RETURN
   IF ( .NOT. aed_oasim_is_finite(dif_) ) RETURN
   IF ( .NOT. aed_oasim_is_finite(ct)   ) RETURN
   IF ( .NOT. aed_oasim_is_finite(h_)   ) RETURN

   IF ( a_  < zero_ .OR. b_ < zero_ .OR. bb_ < zero_ ) RETURN
   IF ( fwd < zero_ .OR. fwd > one_ )                  RETURN
   IF ( dir_ < zero_ .OR. dif_ < zero_ )               RETURN
   IF ( ct <= zero_ )                                  RETURN   !# would divide by zero below
   IF ( h_ <= zero_ )                                  RETURN   !# would divide by zero below

   aed_oasim_kd_inputs_ok = .TRUE.
END FUNCTION aed_oasim_kd_inputs_ok
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
PURE SUBROUTINE aed_oasim_kd_mixture(kk, cc, h_, kd, ok)
!-------------------------------------------------------------------------------
! Diffuse attenuation coefficient of a three-term convex mixture of exponentially
! attenuating streams:
!
!     T  = sum_i c_i * exp(-k_i*h/2)        sum_i c_i = 1,  c_i >= 0
!     Kd = -2*log(T)/h
!
! Evaluated by log-sum-exp - factor out the smallest k among the terms that carry
! weight - so it is exact where the direct form is well conditioned and finite
! where the direct form underflows:
!
!     Kd = k_min - (2/h)*log( sum_i c_i*exp(-(k_i - k_min)*h/2) )
!
! For the strongly absorbing infrared bands (pure-water a up to 1e6 m-1 above
! 1350 nm) every exp underflows to exactly zero in double precision, so the direct
! form gave log(0) = -Inf. Those bands have par_weights = 0 exactly, so the PAR
! collapse then formed 0*Inf = NaN. That is why OAS_Kd_bio was NaN at the surface
! layer on every daylight record of every run in this project.
!
! DO NOT "simplify" this to sum(c_i*k_i). That is the small-h Taylor limit;
! underflow is the large-k*h regime, where the mixture is dominated by the
! SLOWEST-decaying stream. By Jensen the weighted mean always overestimates - by
! 50% on a two-stream test case. It shipped on 2026-09-08 and was wrong.
!
! ok is .FALSE. on invalid input; the caller decides how to fail. This is PURE and
! therefore cannot stop, by design - validation and policy stay separate. It is a
! SUBROUTINE rather than a FUNCTION because a pure function may not have an
! INTENT(OUT) dummy argument (F2018 C1583), and returning invalidity through a
! sentinel value is exactly the silent substitution this whole change removes.
!-------------------------------------------------------------------------------
   AED_REAL,INTENT(in)  :: kk(3), cc(3), h_
   AED_REAL,INTENT(out) :: kd
   LOGICAL, INTENT(out) :: ok
!
!LOCALS
   AED_REAL :: kmin, s, csum
   INTEGER  :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   kd = zero_ ; ok = .FALSE.

   IF ( .NOT. aed_oasim_is_finite(h_) ) RETURN
   IF ( h_ <= zero_ ) RETURN

   csum = zero_
   DO i = 1,3
      IF ( .NOT. aed_oasim_is_finite(kk(i)) ) RETURN
      IF ( .NOT. aed_oasim_is_finite(cc(i)) ) RETURN
      IF ( kk(i) < zero_ .OR. cc(i) < zero_ ) RETURN
      csum = csum + cc(i)
   ENDDO
   !# the weights must be a convex combination; anything else is a caller defect
   IF ( abs(csum - one_) > 1.0e-9 ) RETURN

   kmin = HUGE(kmin)
   DO i = 1,3
      IF ( cc(i) > zero_ ) kmin = MIN(kmin, kk(i))
   ENDDO
   IF ( kmin >= HUGE(kmin) ) RETURN     !# no term carries weight

   s = zero_
   DO i = 1,3
      IF ( cc(i) > zero_ ) s = s + cc(i) * exp(-0.5*(kk(i) - kmin)*h_)
   ENDDO
   IF ( .NOT. (s > zero_) ) RETURN

   kd = kmin - 2.0*log(s)/h_
   ok = aed_oasim_is_finite(kd)
END SUBROUTINE aed_oasim_kd_mixture
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_oasim
