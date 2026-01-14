!###############################################################################
!#                                                                             #
!# aed_atmosphere.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2022-2026 : The University of Western Australia                  #
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
!# Created June 2022                                                           #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# This is currently a copy of aed_dummy except for the vH_O3 routine which    #
!# is an implementation of the Ozone model from :                              #
!#                                                                             #
!# Estimating Atmospheric Ozone for Solar Radiation Models                     #
!#                                              - van Heuklon, T.K. (1979)     #
!# "Solar Energy" Vol 22, pp. 63-68                                            #
!#                                                                             #
!#                                                                             #
!# The intent is to add other routines that contribute to light from the       #
!# atmosphere - eg cloud cover etc.                                            #
!#                                                                             #
!#                                                                             #
!###############################################################################
!                                                                              !
!         .----------------.  .----------------.  .----------------.           !
!         | .--------------. || .--------------. || .--------------. |         !
!         | |      __      | || |  _________   | || | ____    ____ | |         !
!         | |     /  \     | || | |  _   _  |  | || ||_   \  /   _|| |         !
!         | |    / /\ \    | || | |_/ | | \_|  | || |  |   \/   |  | |         !
!         | |   / ____ \   | || |     | |      | || |  | |\  /| |  | |         !
!         | | _/ /    \ \_ | || |    _| |_     | || | _| |_\/_| |_ | |         !
!         | ||____|  |____|| || |   |_____|    | || ||_____||_____|| |         !
!         | |              | || |              | || |              | |         !
!         | '--------------' || '--------------' || '--------------' |         !
!         '----------------'  '----------------'  '----------------'           !
!                                                                              !
!###############################################################################

#include "aed.h"

!
MODULE aed_atmosphere
!-------------------------------------------------------------------------------
! aed_atmosphere --- atmosphere model
!
! The AED module "atmosphere" contains only variables to provide vars required in
! other modules that we didnt provide (usually for debugging purposes)
!-------------------------------------------------------------------------------
   USE aed_core

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_atmosphere_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_atmosphere_data_t
      !# Variable identifiers
      INTEGER :: num_v, num_dv, num_sv, num_dsv
      INTEGER :: id_sine, id_vsine
      INTEGER,ALLOCATABLE :: id_atmosphere_v(:), id_atmosphere_dv(:),           &
                             id_atmosphere_sv(:), id_atmosphere_dsv(:)
      AED_REAL,ALLOCATABLE :: dm_max(:), dm_min(:)
      AED_REAL,ALLOCATABLE :: dm_smax(:), dm_smin(:)

     CONTAINS
         PROCEDURE :: define            => aed_define_atmosphere
         PROCEDURE :: calculate         => aed_calculate_atmosphere
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_atmosphere
   END TYPE

! MODULE GLOBALS
   AED_REAL :: today = 1.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS



!###############################################################################
SUBROUTINE aed_define_atmosphere(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_atmosphere_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i, num_v, num_dv, num_sv, num_dsv

   CHARACTER(len=4),POINTER :: prefix => null()

!  %% NAMELIST   %%  /aed_atmosphere/
!  %% Last Checked 20/08/2021
   CHARACTER(len=40) :: dm_vars(100) = ''
   AED_REAL          :: dm_max(100) = NaN_
   AED_REAL          :: dm_min(100) = NaN_
   AED_REAL          :: dm_init(100) = 0.
   CHARACTER(len=40) :: dm_dvars(100) = ''
   CHARACTER(len=40) :: dm_svars(100) = ''
   AED_REAL          :: dm_smax(100) = NaN_
   AED_REAL          :: dm_smin(100) = NaN_
   AED_REAL          :: dm_sinit(100) = 0.
   CHARACTER(len=40) :: dm_dsvars(100) = ''
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_atmosphere/

   NAMELIST /aed_atmosphere/ dm_vars, dm_max, dm_min, dm_init,             &
                             dm_dvars,                                     &
                             dm_svars, dm_smax, dm_smin, dm_sinit,         &
                             dm_dsvars, diag_level
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_atmosphere initialization"

   num_v = 0 ; num_dv = 0 ; num_sv = 0 ; num_dsv = 0

   ! Read the namelist
   read(namlst,nml=aed_atmosphere,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_atmosphere'

   DO i=1,100 ; IF (dm_vars(i)  .EQ. '' ) THEN  ; num_v  = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dvars(i) .EQ. '' ) THEN  ; num_dv = i-1  ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_svars(i)  .EQ. '' ) THEN ; num_sv  = i-1 ; EXIT ; ENDIF ; ENDDO
   DO i=1,100 ; IF (dm_dsvars(i) .EQ. '' ) THEN ; num_dsv = i-1 ; EXIT ; ENDIF ; ENDDO

   ALLOCATE(data%id_atmosphere_v(num_v))
   ALLOCATE(data%id_atmosphere_dv(num_dv))
   ALLOCATE(data%id_atmosphere_sv(num_sv))
   ALLOCATE(data%id_atmosphere_dsv(num_dsv))

   ALLOCATE(data%dm_min(num_v))   ; ALLOCATE(data%dm_max(num_v))
   ALLOCATE(data%dm_smin(num_sv)) ; ALLOCATE(data%dm_smax(num_sv))

   data%num_v   = num_v
   data%num_dv  = num_dv
   data%num_sv  = num_sv
   data%num_dsv = num_dsv

   CALL aed_set_prefix(prefix)

   ! Register state variables
   DO i=1,data%num_v
      data%id_atmosphere_v(i) = aed_define_variable(dm_vars(i), '', '', dm_init(i), dm_min(i), dm_max(i), 0.)
      data%dm_min(i) = dm_min(i)
      data%dm_max(i) = dm_max(i)
   ENDDO

   DO i=1,data%num_sv
      data%id_atmosphere_sv(i) = aed_define_sheet_variable(dm_svars(i), '', '', dm_sinit(i), dm_smin(i), dm_smax(i), .FALSE.)
      data%dm_smin(i) = dm_smin(i)
      data%dm_smax(i) = dm_smax(i)
   ENDDO

   DO i=1,data%num_dv
      data%id_atmosphere_dv(i) = aed_define_diag_variable(dm_dvars(i), '', '')
   ENDDO

   DO i=1,data%num_dsv
      data%id_atmosphere_dsv(i) = aed_define_sheet_diag_variable(dm_dsvars(i), '', '', .FALSE.)
   ENDDO

   data%id_vsine = aed_define_diag_variable('DUM_vol_sine', 'no units', 'DBG volume sine between 0.0 and 1.0')
   data%id_sine = aed_define_sheet_diag_variable('DUM_sine', 'no units', 'DBG sine wave between 0.0 and 1.0', .FALSE.)
END SUBROUTINE aed_define_atmosphere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_atmosphere(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_atmosphere_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: scale, offs

!-------------------------------------------------------------------------------
!BEGIN
   _DIAG_VAR_(data%id_vsine) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * 0.5) + 0.5

   DO i=1,data%num_v
      scale = (data%dm_max(i) - data%dm_min(i)) / 2.
      offs = data%dm_min(i) + scale
      _STATE_VAR_(data%id_atmosphere_v(i)) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * scale) + offs
   ENDDO
END SUBROUTINE aed_calculate_atmosphere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_atmosphere(data,column,layer_idx)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_atmosphere_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER :: i
   AED_REAL :: scale, offs

!-------------------------------------------------------------------------------
!BEGIN
   IF (layer_idx .EQ. 1) today = today + 1.0/36.5

   _DIAG_VAR_S_(data%id_sine) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * 0.5) + 0.5

   DO i=1,data%num_sv
      scale = (data%dm_smax(i) - data%dm_smin(i)) / 2.
      offs = data%dm_smin(i) + scale
      _STATE_VAR_S_(data%id_atmosphere_sv(i)) = &
        (sin(MOD((today+(layer_idx-1)*10.),365.)/365. * 2 * 3.1415) * scale) + offs
   ENDDO
END SUBROUTINE aed_calculate_benthic_atmosphere
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION vH_O3(long, lat, day)
!-------------------------------------------------------------------------------
! See :
! Estimating Atmospheric Ozone for Solar Radiation Models - van Heuklon, T.K. (1979)
! "Solar Energy" Vol 22, pp. 63-68
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: long, lat, day
!
!LOCALS
   AED_REAL :: O3

   ! Eq 4
   !  O3 = J + {A + C sin[D(E + F)] + G sin[H(lambda + I)]}[sin**2(beta * phi)]
   !
   ! Where Table 1 provides :
   !
   !------------------------------------------------------------
   ! Parm.     N.hemisphere       Both h/sphere    S.hemisphere
   !------------------------------------------------------------
   !  A          150.0                 --             100.0
   ! beta          1.28                --               1.5
   !  C           40.0                 --              30.0
   !  D                              0.9865
   !  E                       Jan 1 = 1.0; Jan 2 = 2.0 ...
   !  F          -30.0                 --             152.625
   !  G                               20.0
   !  H            3.0                 --               2.0
   !  I    if lambda > 0 then 20.0     --             -75.0
   !       else                0.0
   !  J                               235.0
   ! phi        N = +                                 S = -
   ! lambda                        E = +, W = -
   !------------------------------------------------------------

   AED_REAL,PARAMETER :: D = 0.9865, G = 20., J = 235.

   AED_REAL :: A, beta, C, E, F, H, I, phi, lambda

   AED_REAL,PARAMETER :: deg2rad = 3.14159265358979323846 / 180.
!
!-------------------------------------------------------------------------------
!BEGIN
    E = day
    phi = lat
    lambda = -180. + MOD(long + 180., 360.)

    IF (lat > 0.) THEN ! North
       A = 150.
       beta = 1.28
       C = 40.
       F = -30.
       H = 3.
       IF (lambda > 0.) THEN ! E
          I = 20.
       ELSE ! W
          I = 0.
       ENDIF
    ELSE ! South
       A = 100.
       beta = 1.5
       C = 30.
       F = 152.625
       H = 2.
       I = -75.
    ENDIF

    vH_O3 = J + (A + C * sin(D * (E + F) * deg2rad) +                          &
            G * sin(H * (lambda + I) * deg2rad)) * sin(beta * phi * deg2rad)**2
END FUNCTION vH_O3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_atmosphere
