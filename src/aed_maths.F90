!###############################################################################
!#                                                                             #
!# aed_maths.F90                                                               #
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

!#------------------------------##################------------------------------
!### This module is not used yet and is under development. Everything may change
!#------------------------------##################------------------------------
!#

#include "aed.h"

!
MODULE aed_maths

   IMPLICIT NONE

   PRIVATE
!
!  PUBLIC interp_0d, interp_0d_scalar, interp_1d, interp_1d_scalar
   PUBLIC interp
!

   INTERFACE interp
      MODULE PROCEDURE interp_0d
      MODULE PROCEDURE interp_1d
      MODULE PROCEDURE interp_0d_scalar
      MODULE PROCEDURE interp_1d_scalar
   END INTERFACE

CONTAINS

!###############################################################################
SUBROUTINE interp_0d(nsource,x,y,ntarget,targetx,targety)
!-------------------------------------------------------------------------------
! 0D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                      :: nsource, ntarget
   AED_REAL,DIMENSION(nsource),INTENT(in)  :: x, y
   AED_REAL,DIMENSION(ntarget),INTENT(in)  :: targetx
   AED_REAL,DIMENSION(ntarget),INTENT(out) :: targety
!
!LOCALS
   INTEGER  :: i, j
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO j = 1,ntarget
      DO while (i+1<nsource)
         IF (x(i+1)>=targetx(j)) exit
         i = i+1
      ENDDO
      frac = (targetx(j)-x(i))/(x(i+1)-x(i))
      targety(j) = y(i) + frac*(y(i+1)-y(i))
   ENDDO
END SUBROUTINE interp_0d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_0d_scalar(nsource,x,y,targetx,targety)
!-------------------------------------------------------------------------------
! 0D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                     :: nsource
   AED_REAL,DIMENSION(nsource),INTENT(in) :: x, y
   AED_REAL,INTENT(in)                    :: targetx
   AED_REAL,INTENT(out)                   :: targety
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO while (i+1<nsource)
      IF (x(i+1)>=targetx) exit
      i = i+1
   ENDDO
   frac = (targetx-x(i))/(x(i+1)-x(i))
   targety = y(i) + frac*(y(i+1)-y(i))
END SUBROUTINE interp_0d_scalar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_1d(m,nsource,x,y,ntarget,targetx,targety)
!-------------------------------------------------------------------------------
! 1D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                        :: nsource, ntarget, m
   AED_REAL,DIMENSION(nsource),INTENT(in)    :: x
   AED_REAL,DIMENSION(m,nsource),INTENT(in)  :: y
   AED_REAL,DIMENSION(ntarget),INTENT(in)    :: targetx
   AED_REAL,DIMENSION(m,ntarget),INTENT(out) :: targety
!
!LOCALS
   INTEGER  :: i, j
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO j = 1,ntarget
      DO while (i+1<nsource)
         IF (x(i+1)>=targetx(j)) exit
         i = i+1
      ENDDO
      frac = (targetx(j)-x(i))/(x(i+1)-x(i))
      targety(:,j) = y(:,i) + frac*(y(:,i+1)-y(:,i))
   ENDDO
END SUBROUTINE interp_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_1d_scalar(m,nsource,x,y,targetx,targety)
!-------------------------------------------------------------------------------
! 1D interpolation, extrapolates beyond boundaries
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)                       :: nsource, m
   AED_REAL,DIMENSION(nsource),INTENT(in)   :: x
   AED_REAL,DIMENSION(m,nsource),INTENT(in) :: y
   AED_REAL,INTENT(in)                      :: targetx
   AED_REAL,DIMENSION(m),INTENT(out)        :: targety
!
!LOCALS
   INTEGER  :: i
   AED_REAL :: frac

!-------------------------------------------------------------------------------
!BEGIN
   i = 1
   DO while (i+1<nsource)
      IF (x(i+1)>=targetx) exit
      i = i+1
   ENDDO
   frac = (targetx-x(i))/(x(i+1)-x(i))
   targety(:) = y(:,i) + frac*(y(:,i+1)-y(:,i))
END SUBROUTINE interp_1d_scalar
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_maths
