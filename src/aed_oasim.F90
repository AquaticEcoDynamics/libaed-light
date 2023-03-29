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
!#  Copyright 2022 -  The University of Western Australia                      #
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
!### This module is not used yet and is under development. Everything may change
!#------------------------------##################------------------------------
!#                                                                             #
!# An implementation of the 'Ocean-Atmosphere Spectral Irradiance Model (OASIM)'
!#                                                                             #
!# Lazzari, P., Salon, S., Terzić, E., Gregg, W. W., D'Ortenzio, F.,         #
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

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_oasim_data_t
!
   TYPE type_iop
      INTEGER                   :: id_c       ! concentration (could be chl, carbon, or
                                              ! something else, but its product with
                                              ! a or b below should return units m-1)
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a  ! specific absorption (m-1 concentration-1)
      AED_REAL,DIMENSION(:),ALLOCATABLE :: b  ! specific total scattering (m-1 concentration-1)
      AED_REAL :: b_b                         ! ratio of backscattering to total scattering (dimensionless)
   END TYPE
!
   TYPE,extends(aed_model_data_t) :: aed_oasim_data_t
      !# Variable identifiers
      INTEGER :: id_swr, id_uv, id_par, id_par_E, id_par_E_scalar, id_par_J_scalar, id_par_E_dif, id_swr_abs, id_secchi
      INTEGER :: id_swr_sf, id_par_sf, id_uv_sf, id_par_E_sf, id_swr_dif_sf, id_mean_wind_out, id_wind_out, id_zen
      INTEGER :: id_swr_sf_w, id_par_sf_w, id_uv_sf_w, id_par_E_sf_w
      INTEGER :: id_alpha_a, id_beta_a, id_omega_a, id_F_a
      INTEGER :: id_lon, id_lat, id_cloud, id_wind_speed, id_airpres, id_relhum
      INTEGER :: id_lwp, id_O3, id_WV, id_mean_wind_speed, id_visibility, id_air_mass_type
      INTEGER :: id_yearday
      INTEGER :: id_h
      INTEGER :: nlambda
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_surface_band_dir, id_surface_band_dif
      INTEGER,DIMENSION(:),ALLOCATABLE :: id_band_dir, id_band_dif, id_a_band, id_b_band, id_Kd
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: id_a_iop
      AED_REAL,DIMENSION(:),ALLOCATABLE :: lambda, lambda_bounds, par_weights, &
                           par_E_weights, swr_weights, uv_weights, F, lambda_out
      AED_REAL,DIMENSION(:),ALLOCATABLE :: exter
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a_o, a_u, a_v, tau_r
      AED_REAL,DIMENSION(:),ALLOCATABLE :: a_w, b_w, a_w_out, b_w_out
      TYPE (type_iop), ALLOCATABLE :: iops(:)
      INTEGER :: l490_l
      INTEGER :: spectral_output
      LOGICAL :: save_Kd, ozone
      AED_REAL :: cloud, airpres, lwp, O3, WV, visibility, air_mass_type

    CONTAINS
         PROCEDURE :: define            => aed_define_oasim
       !  PROCEDURE :: calculate         => aed_calculate_oasim
         PROCEDURE :: calculate_surface => aed_calculate_surface_oasim
         PROCEDURE :: calculate_column  => aed_calculate_column_oasim
!        PROCEDURE :: light_extinction  => aed_light_extinction_oasim
!        PROCEDURE :: delete            => aed_delete_oasim

   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!   AED_REAL,PARAMETER :: pi = 3.14159265358979323846
!   AED_REAL,PARAMETER :: deg2rad = pi / 180.
!   AED_REAL,PARAMETER :: rad2deg = 180. / pi
!
!

!   INTEGER :: nlambda_astm
!   INTEGER :: nlambda_w

!  AED_REAL,DIMENSION(:),ALLOCATABLE :: a_w, b_w, et_astm, lambda_astm, lambda_w

#include "oasim.inc"


!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed_define_oasim(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
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
   INTEGER :: lambda_ref_iop = 0, a_star_iop = 0, S_iop = 0, b_star_iop = 0, eta_iop = 0, b_b_iop = 0

   INTEGER,PARAMETER :: exter_source = 2

   ! Coefficients for wavelength dependence of foam reflectance (Eqs A11, A12 in Gregg & Casey 2009)
   AED_REAL,PARAMETER :: a0 =  0.9976
   AED_REAL,PARAMETER :: a1 =  0.2194
   AED_REAL,PARAMETER :: a2 =  0.0554
   AED_REAL,PARAMETER :: a3 =  0.0067
   AED_REAL,PARAMETER :: b0 =  5.026
   AED_REAL,PARAMETER :: b1 = -0.0114
   AED_REAL,PARAMETER :: b2 =  9.552e-6
   AED_REAL,PARAMETER :: b3 = -2.698e-9

   AED_REAL,PARAMETER :: Planck     = 6.62606957e-34 ! Planck constant (m2 kg/s)
   AED_REAL,PARAMETER :: lightspeed = 299792458      ! Speed of light (m/s)
   AED_REAL,PARAMETER :: Avogadro   = 6.02214129e23  ! Avogadro constant (/mol)

   AED_REAL :: log_T_w

!  %% NAMELIST   %%  /aed_oasim/
!  %% Last Checked never
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

! %% From Module Globals
!  INTEGER  :: diag_level = 10      ! 0 = no diagnostic outputs
!                                   ! 1 = basic diagnostic outputs
!                                   ! 2 = flux rates, and supporitng
!                                   ! 3 = other metrics
!                                   !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_oasim/

   NAMELIST /aed_oasim/ lambda_method, nlambda, lambda_min, lambda_max,        &
                        n_iop, iop_type, iop_link,                             &
                        spectral_output, nlambda_out, lambda_out,              &
                        lambda, save_Kd,                                       &
                        cloud, airpres, lwp, O3, WV, visibility, air_mass_type,&
                        compute_mean_wind, diag_level

!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_oasim configuration"

   ! Read the namelist
   read(namlst,nml=aed_oasim,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_oasim'

   data%cloud = cloud
   data%airpres = airpres
   data%lwp = lwp
   data%O3 = O3; data%ozone = ozone
   data%WV = WV
   data%visibility = visibility
   data%air_mass_type = air_mass_type

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
      ! OASIM default bandwidths, 33 set in oasmin.inc
      nlambda = nlambda_oasim
      data%nlambda = nlambda
      ALLOCATE(data%lambda(nlambda), data%lambda_bounds(nlambda + 1))
      data%lambda(:) = lambda_oasim
   END SELECT

   ALLOCATE(data%iops(n_iop))
   DO i_iop = 1, n_iop
      ALLOCATE(data%iops(i_iop)%a(nlambda))
      ALLOCATE(data%iops(i_iop)%b(nlambda))
      write(strindex, '(i0)') i_iop

      SELECT CASE (iop_type(i_iop))
      CASE (1) ! diatoms
         CALL interp(size(lambda_diatoms), lambda_diatoms, a_diatoms, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_diatoms), lambda_diatoms, b_diatoms, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.002 ! Gregg & Rousseau 2016 but originally Morel 1988
      CASE (2) ! chlorophytes
         CALL interp(size(lambda_chlorophytes), lambda_chlorophytes, a_chlorophytes, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_chlorophytes), lambda_chlorophytes, b_chlorophytes, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Ahn et al. 1992 as reported in Gregg & Rousseau 2016
      CASE (3) ! cyanobacteria
         CALL interp(size(lambda_cyanobacteria), lambda_cyanobacteria, a_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_cyanobacteria), lambda_cyanobacteria, b_cyanobacteria, nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0032 ! Gregg & Rousseau 2016 but originally Ahn et al. 1992
      CASE (4) ! coccolithophorids
         CALL interp(size(lambda_coccolithophores), lambda_coccolithophores, a_coccolithophores, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_coccolithophores), lambda_coccolithophores, b_coccolithophores, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.00071 * 10 ! Note: 10x Morel 1988 as reported in Gregg & Rousseau 2016
      CASE (5) ! dinoflagellates
         CALL interp(size(lambda_dinoflagellates), lambda_dinoflagellates, a_dinoflagellates, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%a)
         CALL interp(size(lambda_dinoflagellates), lambda_dinoflagellates, b_dinoflagellates, &
                                                                             nlambda, data%lambda, data%iops(i_iop)%b)
         data%iops(i_iop)%b_b = 0.0029 ! Gregg & Rousseau 2016 but originally Morel 1988
      CASE (6) ! detritus (small, as described in Gregg & Rousseau 2016)
         ! Parameters below match the small organic detritus parametrization of Gallegos et al. 2011 (table 2)
         !  - the latter also offers a parametrizaton for large detritus.
         ! Note: Gallegos et al. 2011 constants are specific to dry weight! Did Gregg & Rousseau 2016
         !  misinterpret them as specific to carbon weight?
         ! If so: Babin et al 2003 state 2.6 g DW per g C is representative for suspended OM
         ! NB 12.0107 converts from mg-1 to mmol-1
         ! JB 14/1/2019: adding factor 2.6 (DW/C) discussed above as that seems like to make L4/WCO match much better.
         data%iops(i_iop)%a(:) = 8e-5 * exp(-0.013 * (data%lambda - 440)) * 12.0107 * 2.6
         data%iops(i_iop)%b(:) = 0.00115 * (550. / data%lambda)**0.5 * 12.0107 * 2.6
         data%iops(i_iop)%b_b = 0.005
     !CASE (7) ! PIC
     !   data%iops(i_iop)%a(:) = 0
     !   CALL interp(size(), lambda_w, a_w, nlambda, data%lambda, data%iops(i_iop)%b)
     !   data%iops(i_iop)%b_b = 0.01
      CASE (8) ! CDOC
         ! NB 12.0107 converts from mg-1 to mmol-1
         data%iops(i_iop)%a(:) = 2.98e-4 * exp(-0.014 * (data%lambda - 443)) * 12.0107
         data%iops(i_iop)%b(:) = 0
         data%iops(i_iop)%b_b = 0
      CASE (9) ! Custom carbon-specific absorption and total scattering spectra
         ! NB 12.0107 converts from mg-1 to mmol-1
!        CALL data%get_parameter(a_star_iop, 'a_star_iop'//trim(strindex), 'm2/mg C', 'carbon-mass-specific absorption coefficient for IOP '//trim(strindex)//' at reference wavelength', minimum=0., default=0.)
         IF (a_star_iop /= 0.) THEN
!           CALL data%get_parameter(lambda_ref_iop, 'lambda_a_iop'//trim(strindex), 'nm', 'reference wavelength for absorption by IOP '//trim(strindex))
!           CALL data%get_parameter(S_iop, 'S_iop'//trim(strindex), '-', 'exponent of absorption spectrum for IOP '//trim(strindex), minimum=0.)
            data%iops(i_iop)%a(:) = a_star_iop * exp(-S_iop * (data%lambda - lambda_ref_iop)) * 12.0107
         ELSE
            data%iops(i_iop)%a(:) = 0
         ENDIF

!        CALL data%get_parameter(b_star_iop, 'b_star_iop'//trim(strindex), 'm2/mg C', 'carbon-mass-specific scattering coefficient for IOP '//trim(strindex)//' at reference wavelength', minimum=0., default=0.)
         IF (b_star_iop /= 0.) THEN
!           CALL data%get_parameter(lambda_ref_iop, 'lambda_b_iop'//trim(strindex), 'nm', 'reference wavelength for scattering by IOP '//trim(strindex))
!           CALL data%get_parameter(eta_iop, 'eta_iop'//trim(strindex), '-', 'exponent of scattering spectrum for IOP '//trim(strindex), minimum=0.)
!           CALL data%get_parameter(data%iops(i_iop)%b_b, 'b_b_iop'//trim(strindex), '-', 'backscattering-to-total-scattering ratio for IOP '//trim(strindex), minimum=0.)
            data%iops(i_iop)%b(:) = b_star_iop * (lambda_ref_iop / data%lambda)**eta_iop * 12.0107
         ELSE
            data%iops(i_iop)%b(:) = 0
            data%iops(i_iop)%b_b = 0
         ENDIF
      END SELECT

      ! Protect against negative coefficients caused by extrapolation beyond source spectrum boundaries.
      data%iops(i_iop)%a(:) = max(data%iops(i_iop)%a, 0.)
      data%iops(i_iop)%b(:) = max(data%iops(i_iop)%b, 0.)

      ! Link to concentration metric to allow us to convert *specific* absorption/scattering into
      ! actual absorption and scattering (in m-1)
      IF (iop_type(i_iop) >=1 .and. iop_type(i_iop) <= 5) THEN
         ! Phytoplankton: chlorophyll-specific absorption and scattering
!         data%iops(i_iop)%id_c = aed_define_variable('iop' // trim(strindex) // '_chl', 'mg Chl m-3', &
!                                                                             'chlorophyll in IOP ' // trim(strindex))
!        CALL data%request_coupling_to_model(data%iops(i_iop)%id_c, 'iop' // trim(strindex), &
!                                                                            type_bulk_standard_variable(name='total_chlorophyll'))
         data%iops(i_iop)%id_c = aed_locate_variable(TRIM(iop_link(i_iop)))
      ELSE
         ! POM/DOM/PIC: carbon-specific absorption and scattering
!         data%iops(i_iop)%id_c = aed_define_variable('iop' // trim(strindex) // '_c', &
!                                                                             'mmol C m-3', 'carbon in IOP ' // trim(strindex))
!        CALL data%request_coupling_to_model(data%iops(i_iop)%id_c, 'iop' // trim(strindex), standard_variables%total_carbon)
         data%iops(i_iop)%id_c = aed_locate_variable(TRIM(iop_link(i_iop)))
      ENDIF
   ENDDO

   ! Find wavelength bounds of photosynthetically active radiation
   ALLOCATE(data%par_weights(nlambda))
   ALLOCATE(data%swr_weights(nlambda))
   ALLOCATE(data%uv_weights(nlambda))
   ALLOCATE(data%par_E_weights(nlambda))
   CALL calculate_integral_weights(400., 700., nlambda, data%lambda, data%par_weights)
   CALL calculate_integral_weights(300., 4000., nlambda, data%lambda, data%swr_weights)
   CALL calculate_integral_weights(300., 400., nlambda, data%lambda, data%uv_weights)
   data%par_E_weights(:) = data%par_weights * data%lambda /(Planck*lightspeed)/Avogadro*1e-3
                                        ! divide by 1e9 to go from nm to m, multiply by 1e6 to go from mol to umol

   data%id_lon = aed_locate_global('longitude')
   data%id_lat = aed_locate_global('latitude')
   data%id_yearday = aed_locate_global('yearday') 

   data%id_h = aed_locate_global('layer_ht') 
   data%id_wind_speed = aed_locate_sheet_global('wind_speed')
   data%id_relhum = aed_locate_sheet_global('humidity')
  
  ! data%id_cloud = aed_locate_sheet_global('cloud')
  ! data%id_airpres = aed_locate_sheet_global('air_press')
  ! data%id_lwp = aed_locate_sheet_global('lwp')
  ! data%id_WV = aed_locate_sheet_global('WV ')
  ! data%id_visibility = aed_locate_sheet_global('visibility')
  ! data%id_air_mass_type = aed_locate_sheet_global('air_mass_type')

  ! Constant OR Van Heuklon (1979) function to account for seasonal / geographical variation in O3
   data%id_O3 = 0  !  data%id_O3 = aed_locate_sheet_global('atmosphere_mass_content_of_ozone')
   IF(ozone) data%id_O3 = aed_define_sheet_diag_variable('ozone', 'kg/m2', 'atmosphere_mass_content_of_ozone')


!! MAKE PARAM CAB ?
!!  data%id_cloud = aed_locate_sheet_global('cloud_area_fraction')
!   data%id_cloud = aed_define_sheet_diag_variable('cloud_area_fraction', '-', '-')
!! MAKE PARAM CAB ?
!!  data%id_airpres = aed_locate_sheet_global('surface_air_pressure')
!   data%id_airpres = aed_define_sheet_diag_variable('surface_air_pressure', '-', '-')
!! MAKE PARAM CAB ?
!!  data%id_yearday = aed_locate_sheet_global('number_of_days_since_start_of_the_year')
!   data%id_yearday = aed_define_sheet_diag_variable('number_of_days_since_start_of_the_year', '-', '-')
!
!! MAKE PARAM CAB ?
!!  data%id_lwp = aed_locate_sheet_global('atmosphere_mass_content_of_cloud_liquid_water')
!   data%id_lwp = aed_define_sheet_diag_variable('atmosphere_mass_content_of_cloud_liquid_water', '-', '-')
!! MAKE PARAM CAB ?
!! MAKE PARAM CAB ?
!!  data%id_WV = aed_locate_sheet_global('atmosphere_mass_content_of_water_vapor')
!   data%id_WV = aed_define_sheet_diag_variable('atmosphere_mass_content_of_water_vapor', '-', '-')
!! MAKE PARAM CAB ?
!!  data%id_visibility = aed_locate_sheet_global('visibility_in_air')
!   data%id_visibility = aed_define_sheet_diag_variable('visibility_in_air', '-', '-')
!! MAKE PARAM CAB ?
!!  data%id_air_mass_type = aed_locate_sheet_global('aerosol_air_mass_type')
!   data%id_air_mass_type = aed_define_sheet_diag_variable('aerosol_air_mass_type', '-', '-')

   IF (compute_mean_wind) THEN
      data%id_mean_wind_out = aed_define_sheet_diag_variable('mean_wind', 'm/s', 'daily mean wind speed')
      ! CAB need something here ....
      data%id_mean_wind_speed = aed_define_sheet_diag_variable('mean_wind', 'm/s', 'daily mean wind speed')
   ELSE
      data%id_mean_wind_speed = aed_locate_sheet_global('mean_wind')
   ENDIF

   data%id_zen = aed_define_sheet_diag_variable('zen', 'degrees', 'zenith angle')

   ! Aerosol properties
   data%id_alpha_a =  aed_define_sheet_diag_variable('alpha_a', '-', 'aerosol Angstrom exponent')
   data%id_beta_a =   aed_define_sheet_diag_variable('beta_a', '-', 'aerosol scale factor for optical thickness')
   data%id_omega_a =  aed_define_sheet_diag_variable('omega_a', '-', 'aerosol single scattering albedo')
   data%id_F_a =      aed_define_sheet_diag_variable('F_a', '-', 'aerosol forward scattering probability')
   data%id_wind_out = aed_define_sheet_diag_variable('wind', 'm/s', 'wind speed')

   ! Horizontal downwelling irradiance just above the water surface (BEFORE reflection by water surface)
   data%id_swr_sf =     aed_define_sheet_diag_variable('swr_sf', 'W/m^2', 'downwelling shortwave flux in air')
   data%id_swr_dif_sf = aed_define_sheet_diag_variable('swr_dif_sf', 'W/m^2', 'diffuse downwelling shortwave flux in air')
   data%id_uv_sf =      aed_define_sheet_diag_variable('uv_sf', 'W/m^2', 'downwelling ultraviolet radiative flux in air')
   data%id_par_sf =     aed_define_sheet_diag_variable('par_sf', 'W/m^2', 'downwelling photosynthetic radiative flux in air')
   data%id_par_E_sf =   aed_define_sheet_diag_variable('par_E_sf', 'umol/m^2/s', 'downwelling photosynthetic photon flux in air')

   ! Horizontal downwelling irradiance just below the water surface (AFTER reflection by water surface)
   data%id_swr_sf_w =   aed_define_sheet_diag_variable('swr_sf_w', 'W/m^2', 'downwelling shortwave flux in water')
   data%id_uv_sf_w =    aed_define_sheet_diag_variable('uv_sf_w', 'W/m^2', 'downwelling ultraviolet radiative flux in water')
   data%id_par_sf_w =   aed_define_sheet_diag_variable('par_sf_w', 'W/m^2', 'downwelling photosynthetic radiative flux in water')
   data%id_par_E_sf_w = aed_define_sheet_diag_variable('par_E_sf_w', 'umol/m^2/s', 'downwelling photosynthetic photon flux in water')

   ! Scalar downwelling irradiance within the water column
   data%id_swr =          aed_define_diag_variable('swr', 'W/m^2', 'downwelling shortwave flux')
   data%id_uv =           aed_define_diag_variable('uv', 'W/m^2', 'downwelling ultraviolet radiative flux')
   data%id_par =          aed_define_diag_variable('par', 'W/m^2', 'downwelling photosynthetic radiative flux')
   data%id_par_E =        aed_define_diag_variable('par_E', 'umol/m^2/s', 'downwelling photosynthetic photon flux')
   data%id_par_J_scalar = aed_define_diag_variable('par_J_scalar','W/m^2', 'scalar downwelling photosynthetic radiative flux')
   data%id_par_E_scalar = aed_define_diag_variable('par_E_scalar','umol/m^2/s', 'scalar downwelling photosynthetic photon flux')
   data%id_swr_abs =      aed_define_diag_variable('swr_abs', 'W/m^2', 'absorption of shortwave energy in layer')
   data%id_par_E_dif =    aed_define_diag_variable('par_E_dif', 'W/m^2', 'diffusive downwelling photosynthetic photon flux')

   data%id_secchi =       aed_define_sheet_diag_variable('secchi', 'm', 'Secchi depth (1.7/Kd 490)')

   ! Interpolate absorption and scattering spectra to user wavelength grid
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

   !data%l490_l = nlambda - 1
   !DO l = 1, nlambda - 1
   !   IF (data%lambda(l) >= 490.) THEN
   !      data%l490_l = l
   !      exit
   !   ENDIF
   !ENDDO

   data%spectral_output = spectral_output
   SELECT CASE (data%spectral_output)
   CASE (1)
      ALLOCATE(data%lambda_out(nlambda))
      data%lambda_out(:) = data%lambda
   CASE (2)
      ALLOCATE(data%lambda_out(nlambda_out))
      DO l = 1, nlambda_out
         data%lambda_out(l) = lambda_out(l)
      ENDDO
      ALLOCATE(data%a_w_out(nlambda_out), data%b_w_out(nlambda_out))
      CALL interp(nlambda_w, lambda_w, a_w, nlambda_out, data%lambda_out, data%a_w_out)
      CALL interp(nlambda_w, lambda_w, b_w, nlambda_out, data%lambda_out, data%b_w_out)
   END SELECT
   data%save_Kd = save_Kd

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
                                                'attenuation coefficient')
   ENDIF


   print *,'config OK'
END SUBROUTINE aed_define_oasim
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_surface_oasim(data,column,layer_idx)
   !-------------------------------------------------------------------------------
   ! Atmospheric component of the aed implementation of the OASIM model
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      CLASS (aed_oasim_data_t),INTENT(in) :: data
      TYPE (aed_column_t),INTENT(inout) :: column(:)
      INTEGER,INTENT(in) :: layer_idx
   !
   !LOCALS
      ! Environment
      AED_REAL :: temp, salt, wind, depth
      AED_REAL :: vel = 0.0001
   
      ! State
      AED_REAL :: oxy
   
   !
   !-------------------------------------------------------------------------------
   !BEGIN


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
   
      AED_REAL :: longitude, latitude, yearday, cloud_cover, wind_speed, airpres, relhum, LWP, water_vapour, WV, WM, visibility, AM
      AED_REAL :: days, hour, theta, costheta, alpha_a, beta_a
      INTEGER  :: l
      AED_REAL :: M, M_prime, M_oz
      AED_REAL :: O3
      AED_REAL,DIMENSION(data%nlambda) :: direct, diffuse, spectrum, Kd  ! Spectra at top of the water
                                               ! column (with refraction and reflection accounted for)
      AED_REAL :: par_J, swr_J, uv_J, par_E, F_a, omega_a
      AED_REAL,DIMENSION(data%nlambda) :: tau_a, T_a, T_oz, T_w, T_u, T_r, T_aa, T_as
      AED_REAL,DIMENSION(data%nlambda) :: T_g, T_dclr, T_sclr, T_dcld, T_scld
      AED_REAL,DIMENSION(data%nlambda) :: rho_d, rho_s
   
      AED_REAL,DIMENSION(data%nlambda) :: a, b, b_b, a_iop
      AED_REAL,DIMENSION(data%nlambda) :: f_att_d, f_att_s, f_prod_s
      AED_REAL, allocatable :: spectrum_out(:)
      INTEGER  :: i_iop
      AED_REAL :: c_iop, h, swr_top, costheta_r, dir_frac
   
      
   !-------------------------------------------------------------------------------
   !BEGIN
      !_MOVE_TO_TOP_

      layer_idx = layer_map(1) ! Start at the top of the column

      IF (data%spectral_output == 2) ALLOCATE(spectrum_out(size(data%lambda_out)))
   
      longitude   = _STATE_VAR_S_(data%id_lon)
      latitude    = _STATE_VAR_S_(data%id_lat)
      yearday     = _STATE_VAR_S_(data%id_yearday) 
      relhum      = _STATE_VAR_S_(data%id_relhum)          ! Relative humidity (-)
      wind_speed  = _STATE_VAR_S_(data%id_wind_speed)      ! Wind speed @ 10 m above surface (m/s)
      WM          =  _DIAG_VAR_S_(data%id_mean_wind_speed) ! Daily mean wind speed @ 10 m above surface (m/s)
   
      AM          = data%air_mass_type !_STATE_VAR_S_(data%id_air_mass_type)   ! Aerosol air mass type (1: open ocean, 10: continental)
      visibility  = data%visibility ! _STATE_VAR_S_(data%id_visibility)     ! Visibility (m)
      airpres     = data%airpres ! _STATE_VAR_S_(data%id_airpres)       ! Surface air pressure (Pa)
      cloud_cover = data%cloud ! _STATE_VAR_S_(data%id_cloud)        ! Cloud cover (fraction, 0-1)
      LWP         = data%LWP ! _STATE_VAR_S_(data%id_lwp)        ! Cloud liquid water content (kg m-2)
      water_vapour= data%WV ! _STATE_VAR_S_(data%id_wv)       ! Total precipitable water vapour (kg m-2) - equivalent to mm
      
      O3          = data%O3    ! Ozone content (kg m-2)
      IF(data%ozone) THEN 
         O3 = vH_O3(longitude, latitude, yearday)
         _DIAG_VAR_S_(data%id_O3) = O3
      ENDIF
   
      ! For debugging
      _DIAG_VAR_S_(data%id_wind_out) = wind_speed
      IF (data%id_mean_wind_out > 0) _DIAG_VAR_S_(data%id_mean_wind_out) = WM
   
      IF (cloud_cover > 0) LWP = LWP / cloud_cover  ! LWP is the mean density over a grid box (Jorn: ECMWF pers comm 26/2/2019), but we want the mean density per cloud-covered area
      WV = water_vapour / 10                        ! from kg m-2 to cm
      O3 = O3 * (1000 / 48.) / 0.4462         ! from kg m-2 to mol m-2, then from mol m-2 to atm cm (Basher 1982)
      days = floor(yearday)
      hour = mod(yearday, 1.0) * 24.0
   
      ! Calculate zenith angle (in radians)
      theta = zenith_angle(days, hour, longitude, latitude)
      theta = min(theta, 0.5 * pi)  ! Restrict the input zenith angle between 0 and pi/2
      _DIAG_VAR_S_(data%id_zen) = rad2deg * theta
      costheta = cos(theta)
   
      print *,'yearday', yearday,days, hour, longitude, latitude, theta, costheta
   
      ! Atmospheric path length, a.k.a. relative air mass (Eq 3 Bird 1984, Eq 5 Bird & Riordan 1986, Eq 13 Gregg & Carder 1990, Eq A5 in Casey & Gregg 2009)
      ! Note this should always exceed 1, but as it is an approximation it does not near theta -> 0 (Tomasi et al. 1998 p14). Hence the max operator.
      M = max(1., 1. / (costheta + 0.15 * (93.885 - theta * rad2deg)**(-1.253)))
   
      ! Pressure-corrected atmospheric path length (Eq A6 Casey & Gregg 2009)
      M_prime = M * airpres / pres0
   
      ! Atmospheric path length for ozone
      ! Eq 10, Bird & Riordan 1986, Eq 14 Gregg & Carder 1990; NB 6370 is the earth's radius in km
      ! See also Tomasi et al. 1998 Eq 5
      M_oz = (1. + H_oz / 6370.) / sqrt(costheta**2 + 2 * H_oz / 6370.)
   
      ! Transmittance due to ozone absorption (Eq 8 Bird 1984, Eq 9 Bird & Riordan 1986, Eq 17 Gregg & Carder 1990)
      T_oz = exp(-data%a_o * O3 * M_oz)
   
      ! Transmittance due to water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
      ! Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
      T_w = exp((-0.2385 * data%a_v * WV * M) / (1. + 20.07 * data%a_v * WV * M)**0.45)
   
      ! Transmittance due to uniformly mixed gas absorption - SHOULD use pressure corrected airmass
      ! Eq 10 Bird 1984, Eq 11 Bird and Riordan 1986, Eq 18 Gregg & Carder 1990
      ! Bird & Riordan use 118.93 rather than 118.3, but state 118.3 should be used in the future.
      T_u = exp(-1.41 * data%a_u * M_prime / (1. + 118.3 * data%a_u * M_prime)**0.45)
   
      ! Transmittance terms that apply for both cloudy and clear skies.
      T_g = T_oz * T_w * T_u
   
      ! -------------------
      ! clear skies part
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
      ! cloudy skies part
      ! -------------------
   
      ! Transmittance due to absorption and scattering by clouds
      CALL slingo(costheta, max(0., LWP) * 1000, r_e, data%nlambda, data%lambda, T_dcld, T_scld)
   
      ! Diffuse and direct irradiance streams (Eqs 1, 2 Gregg & Casey 2009)
      ! These combine terms for clear and cloudy skies, weighted by cloud cover fraction
      direct  = data%exter * costheta * T_g * ((1. - cloud_cover) * T_dclr + cloud_cover * T_dcld)
      diffuse = data%exter * costheta * T_g * ((1. - cloud_cover) * T_sclr + cloud_cover * T_scld)
   
      spectrum = direct + diffuse
      par_J = sum(data%par_weights * spectrum)
      swr_J = sum(data%swr_weights * spectrum)
      par_E = sum(data%par_E_weights * spectrum)
      uv_J  = sum(data%uv_weights * spectrum)
      _DIAG_VAR_S_(data%id_par_E_sf) = par_E ! Photosynthetically Active Radiation, PAR (umol/m2/s)
      _DIAG_VAR_S_(data%id_par_sf)   = par_J ! Photosynthetically Active Radiation, PAR (W/m2)
      _DIAG_VAR_S_(data%id_swr_sf)   = swr_J ! Total shortwave radiation (W/m2), SW [up to 4000 nm]
      _DIAG_VAR_S_(data%id_uv_sf)    = uv_J  ! Ultraviolet Radiation, UV (W/m2)
   
      swr_J = sum(data%swr_weights * diffuse)
      _DIAG_VAR_S_(data%id_swr_dif_sf) = swr_J ! Diffuse shortwave radiation (W/m2), SW [up to 4000 nm]
   
   
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
   
      ! Sea surface reflectance
      CALL reflectance(data%nlambda, data%F, theta, wind_speed, rho_d, rho_s, costheta_r)
   
      ! Incorporate the loss due to reflectance
      direct = direct * (1. - rho_d)
      diffuse = diffuse * (1. - rho_s)
      spectrum = direct + diffuse
   
      par_E = sum(data%par_E_weights * spectrum)
      par_J = sum(data%par_weights * spectrum)
      swr_J = sum(data%swr_weights * spectrum)
      uv_J  = sum(data%uv_weights * spectrum)
      _DIAG_VAR_S_(data%id_par_sf_w) = par_J  ! Photosynthetically Active Radiation (W/m2)
      _DIAG_VAR_S_(data%id_par_E_sf_w) =par_E ! Photosynthetically Active Radiation (umol/m2/s)
      _DIAG_VAR_S_(data%id_swr_sf_w) =  swr_J ! Total shortwave radiation (W/m2) [up to 4000 nm]
      _DIAG_VAR_S_(data%id_uv_sf_w) = uv_J    ! UV (W/m2)
   
   
   !CAB   _DOWNWARD_LOOP_BEGIN_

      DO layer = 1,SIZE(layer_map)
         layer_idx = layer_map(layer)

         !print *,'layer_idx',layer,layer_idx, swr_J

         ! Save downwelling shortwave flux at top of the layer
         swr_top = swr_J
   
         ! Compute absorption, total scattering and backscattering in current layer from IOPs
         a = data%a_w
         b = data%b_w
         b_b = 0.5 * data%b_w
         DO i_iop = 1, size(data%iops)
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
            b = b + c_iop * data%iops(i_iop)%b
            b_b = b_b + c_iop * data%iops(i_iop)%b_b * data%iops(i_iop)%b
         ENDDO
   
         ! Transmissivity of direct/diffuse attentuation and conversion from direct to diffuse - for one half of the layer
         h = _STATE_VAR_(data%id_h)
         f_att_d = exp(-0.5 * (a + b) * h / costheta_r)         ! Gregg & Rousseau 2016 Eq 8
         f_att_s = exp(-0.5 * (a + r_s * b_b) * h / mcosthetas) ! Gregg & Rousseau 2016 Eq 9
         f_prod_s = exp(-0.5 * a * h / costheta_r) - f_att_d    ! Gregg & Rousseau 2016 Eq 14 but not accounting for backscattered fraction
   
         IF (data%save_Kd .and. data%spectral_output /= 0) THEN
            DO l = 1, data%nlambda
               !Kd(l) = 2 * (log(direct(l) + diffuse(l)) - log(direct(l) * (f_att_d(l) + f_prod_s(l)) + diffuse(l) * f_att_s(l))) / h
               IF (direct(l) + diffuse(l) > 0) THEN
                  ! Direct and diffuse stream
                  dir_frac = direct(l) / (direct(l) + diffuse(l))
                  Kd(l) = - 2 *log(dir_frac * (f_att_d(l) + f_prod_s(l)) + (1.0 - dir_frac) * f_att_s(l)) / h
               ELSE
                  ! Only diffuse stream
                  Kd(l) = (a(l) + r_s * b_b(l)) / mcosthetas
               ENDIF
            ENDDO
         ENDIF

         IF (data%save_Kd) &
            _DIAG_VAR_(data%id_Kd(size(data%lambda_out)+1)) = SUM(Kd(1:data%nlambda))
   
         ! From top to centre of layer
         direct = direct * f_att_d
         diffuse = diffuse * f_att_s + direct * f_prod_s
         spectrum = direct + diffuse
   
         par_E = sum(data%par_E_weights * spectrum)
         par_J = sum(data%par_weights * spectrum)
         swr_J = sum(data%swr_weights * spectrum)
         uv_J  = sum(data%uv_weights * spectrum)
         _DIAG_VAR_(data%id_par_E) = par_E ! Photosynthetically Active Radiation (umol/m2/s)
         _DIAG_VAR_(data%id_par) = par_J   ! Photosynthetically Active Radiation (W/m2)
         _DIAG_VAR_(data%id_swr) =  swr_J  ! Total shortwave radiation (W/m2) [up to 4000 nm]
         _DIAG_VAR_(data%id_uv) = uv_J     ! UV (W/m2)
         _DIAG_VAR_(data%id_par_E_dif) = sum(data%par_E_weights * diffuse) ! Diffuse Photosynthetically Active photon flux (umol/m2/s)
   
         ! Compute scalar PAR as experienced by phytoplankton
         spectrum = direct / costheta_r + diffuse / mcosthetas
         par_J = sum(data%par_weights * spectrum)
         par_E = sum(data%par_E_weights * spectrum)
         _DIAG_VAR_(data%id_par_J_scalar) = par_J ! Scalar Photosynthetically Active Radiation (W/m2)
         _DIAG_VAR_(data%id_par_E_scalar) = par_E ! Scalar Photosynthetically Active photon flux (umol/m2/s)
   
         ! From centre to bottom of layer
         direct = direct * f_att_d
         diffuse = diffuse * f_att_s + direct * f_prod_s
         spectrum = direct + diffuse
   
         ! Save spectrally resolved outputs
         SELECT CASE (data%spectral_output)
         CASE (1)
            DO l = 1, data%nlambda
               _DIAG_VAR_(data%id_a_band(l)) = a(l) - data%a_w(l)
               _DIAG_VAR_(data%id_b_band(l)) = b(l) - data%b_w(l)
                IF (data%save_Kd) Kd(l) = _DIAG_VAR_(data%id_Kd(l))
            ENDDO
         CASE (2)
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
         swr_J = sum(data%swr_weights * spectrum)
         _DIAG_VAR_(data%id_swr_abs) = swr_top - swr_J
      END DO

   !CAB   _VERTICAL_LOOP_END_
   
      ! Put remaining shortwave in bottom layer
      ! (assumes all light is absorbed by sediment and injected into water column)
   !CAB   _MOVE_TO_BOTTOM_
   !MH   _DIAG_VAR_(data%id_swr_abs) = swr_J
   
         _DIAG_VAR_S_(data%id_secchi) = 0.
      
END SUBROUTINE aed_calculate_column_oasim
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
   AED_REAL,PARAMETER :: n_w = 1.341

   AED_REAL :: C_D, tau, rho_f_W, rho_f(nlambda)
   AED_REAL :: b, rho_dsp, rho_ssp, theta_r
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


END MODULE aed_oasim
