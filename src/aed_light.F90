!###############################################################################
!#                                                                             #
!# aed_light.F90                                                               #
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
!# Created April 2022                                                          #
!#                                                                             #
!###############################################################################

!#------------------------------##################------------------------------
!### This module is not used yet and is under development. Everything may change
!#------------------------------##################------------------------------

#include "aed.h"

!
MODULE aed_light

   USE aed_core
!  USE aed_common

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC aed_light_data_t
!
   TYPE,extends(aed_model_data_t) :: aed_light_data_t
      !# Variable identifiers
      INTEGER :: id_light, id_light_par, id_light_uv, id_light_extc
      INTEGER :: id_par, id_extc
      LOGICAL :: outputLight

     CONTAINS
         PROCEDURE :: define            => aed_define_light
         PROCEDURE :: calculate         => aed_calculate_light
!        PROCEDURE :: mobility          => aed_mobility_light
!        PROCEDURE :: light_extinction  => aed_light_extinction_light
!        PROCEDURE :: delete            => aed_delete_light

   END TYPE

! MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

   LOGICAL :: link_water_clarity = .FALSE.

!===============================================================================
CONTAINS

!###############################################################################
SUBROUTINE aed_define_light(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the AED model
!
!  Here, the aed namelist is read and te variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: namlst
   CLASS (aed_light_data_t),INTENT(inout) :: data

!
!LOCALS
   INTEGER :: status
   INTEGER :: i, num_tn,num_tkn,num_tp,num_toc,num_tss,num_turb,num_tfe,num_tal

!  %% NAMELIST   %%  /aed_light/
!  %% Last Checked 20/08/2021
   AED_REAL :: initial_light = 1.0
   AED_REAL :: initial_par   = 1.0
   AED_REAL :: initial_uv    = 1.0
   AED_REAL :: initial_extc  = 1.0
   LOGICAL  :: outputLight   = .FALSE.
! %% From Module Globals
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
!                                             ! 1 = basic diagnostic outputs
!                                             ! 2 = flux rates, and supporitng
!                                             ! 3 = other metrics
!                                             !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_light/

   NAMELIST /aed_light/ initial_light, initial_par, initial_uv, initial_extc,  &
                        outputLight
!
!-------------------------------------------------------------------------------
!BEGIN
   print *,"        aed_light configuration"

   ! Read the namelist
   read(namlst,nml=aed_light,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_light'

   data%outputLight = outputLight
   IF (data%outputLight) THEN
     data%id_light = aed_define_diag_variable('light',         &
                     'W/m2', 'Shortwave light flux')

     data%id_light_par = aed_define_diag_variable('par',            &
                     'W/m2', 'Photosynthetically active light flux')

     data%id_light_uv = aed_define_diag_variable('uv',               &
                     'W/m2', 'Ultraviolet light flux')

     data%id_light_extc = aed_define_diag_variable('extc',           &
                     'W/m2', 'Light extinction coefficient')

     ! Register environmental dependencies
     data%id_par = aed_locate_global('par')
     data%id_extc = aed_locate_global('extc_coef')
   END IF
END SUBROUTINE aed_define_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_light(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of aed_light model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_light_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   AED_REAL :: val

!-------------------------------------------------------------------------------
!BEGIN
   ! Light and Extinction Coefficient
   IF (data%outputLight) THEN
     val = _STATE_VAR_(data%id_par)
     _DIAG_VAR_(data%id_light) =  val/0.43
     _DIAG_VAR_(data%id_light_par)   =  val
     _DIAG_VAR_(data%id_light_uv)    = (val/0.43)*0.05
     val = _STATE_VAR_(data%id_extc)
     _DIAG_VAR_(data%id_light_extc)  =  val
   END IF
END SUBROUTINE aed_calculate_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE BioExtinction(column, count, extc)
!-------------------------------------------------------------------------------
!
! Calculate the specific light attenuation additions due to AED modules
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(inout) :: extc(:)
!
#if 0
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: localext
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_

   CALL aed_light_extinction(column, 1, localext)
   IF (link_water_clarity) THEN
      extc(1) = localext
   ELSE
      extc(1) = localext + Kw
   END IF

   IF (count <= 1) RETURN

   DO i = 2, count
      CALL aed_light_extinction(column, i, localext)
      IF (link_water_clarity) THEN
         extc(i) = localext
      ELSE
         extc(i) = localext + Kw
      END IF
   ENDDO
#endif
END SUBROUTINE BioExtinction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE Light(column, count, Io, extc, par_, h_)
!-------------------------------------------------------------------------------
!
! Calculate photosynthetically active radiation over entire column
! based on surface radiation, and background and biotic extinction.
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,  INTENT(in)    :: count
   AED_REAL, INTENT(in)    :: Io
   AED_REAL, INTENT(inout) :: extc(:)
   AED_REAL, INTENT(inout) :: par_(:)
   AED_REAL, INTENT(inout) :: h_(:)
!
!LOCAL VARIABLES:
   INTEGER :: i
   AED_REAL :: zz, localext, localshade
!
!-------------------------------------------------------------------------------
!BEGIN
   zz = zero_
   localext = zero_

   CALL BioExtinction(column,count,extc)

   localext = extc(1)
   zz = 0.001 !0.5*h_(1)    !MH: assume top of layer
   par_(1) = 0.45 * Io * EXP( -(localext) * zz )

   IF (count <= 1) RETURN

   DO i = 2, count
      localext = extc(i)

      !zz = zz + 0.5*h_(i)
      zz = h_(i)
      par_(i) = par_(i-1) * EXP( -(localext) * zz )
   ENDDO
END SUBROUTINE Light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed_column_t), INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: nlev
!
!LOCALS
#if 0
   INTEGER :: i
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN
   localext = zero_; localext_up = zero_

   ! Surface Kd
   CALL aed_light_extinction(column, nlev, localext)

   ! Surface PAR
   par(nlev) = par_fraction * rad(nlev) * EXP( -(lKw+localext)*1e-6*dz(nlev) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = (nlev-1),1,-1
      localext_up = localext
      CALL aed_light_extinction(column, i, localext)

      par(i) = par(i+1) * EXP( -(lKw + localext_up) * dz(i+1) )

      IF (bioshade_feedback) extc_coef(i) = lKw + localext
   ENDDO
#endif
END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_light
