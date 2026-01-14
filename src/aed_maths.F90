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
!#

#include "aed.h"

!
MODULE aed_maths

   IMPLICIT NONE

   PRIVATE
!
!  PUBLIC interp_0d, interp_0d_scalar, interp_1d, interp_1d_scalar
   PUBLIC interp, calculate_integral_weights, pi, deg2rad, rad2deg
!

   INTERFACE interp
      MODULE PROCEDURE interp_0d
      MODULE PROCEDURE interp_1d
      MODULE PROCEDURE interp_0d_scalar
      MODULE PROCEDURE interp_1d_scalar
   END INTERFACE

AED_REAL,PARAMETER :: pi = 3.14159265358979323846
AED_REAL,PARAMETER :: deg2rad = pi / 180.
AED_REAL,PARAMETER :: rad2deg = 180. / pi

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


!###############################################################################
SUBROUTINE calculate_integral_weights(xl, xr, n, x, w)
   !-------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER, INTENT(in)  :: n
      AED_REAL,INTENT(in)  :: x(n), xl, xr
      AED_REAL,INTENT(out) :: w(n)
   !
   !LOCALS
      INTEGER  :: i
      AED_REAL :: deltax, f
   !-------------------------------------------------------------------------------
   !BEGIN
      w = 0.

      IF (x(1) > xr) THEN
         ! All points to right of desired range.
         w(1) = xr - xl
         return
      ELSEIF (x(n) < xl) THEN
         ! All points to left of desired range.
         w(n) =  xr - xl
         return
      ENDIF

      DO i = 2, n
         deltax = x(i) - x(i-1)
         IF (x(i) >= xl .and. x(i) <= xr) THEN
            ! Right-hand point in desired range.
            IF (x(i-1) >= xl) THEN
               ! Whole segment in range
               ! Integral: 0.5*(y1+y2)*(x2-x1)
               w(i-1:i) = w(i-1:i) + 0.5 * deltax
            ELSE
               ! Segment crosses left boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xl-x1) = y1*(1-(xl-x1)/(x2-x1)) + y2*(xl-x1)/(x2-x1)
               ! integral = 0.5*(yb+y2)*(x2-xl)
               ! yb+y2 = y1*(1-(xl-x1)/(x2-x1)) + y2*(1+(xl-x1)/(x2-x1))
               f = (xl - x(i-1)) / deltax
               w(i-1) = w(i-1) + (1.0 - f) * (x(i) - xl) * 0.5
               w(i  ) = w(i  ) + (1.0 + f) * (x(i) - xl) * 0.5
            ENDIF
         ELSEIF (x(i) > xr .and. x(i-1) < xr) THEN
            ! Right-hand point beyond desired range, left-hand point before right range boundary.
            IF (x(i-1) >= xl) THEN
               ! Segment crosses right boundary
               ! y at boundary: yb = y1 + (y2-y1)/(x2-x1)*(xr-x1) = y1*(1-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               ! integral = 0.5*(y1+yb)*(xr-x1)
               ! y1+yb = y1*(2-(xr-x1)/(x2-x1)) + y2*(xr-x1)/(x2-x1)
               f = (xr - x(i-1)) / deltax
               w(i-1) = w(i-1) + (2.0 - f) * (xr - x(i-1)) * 0.5
               w(i  ) = w(i  ) + (         f) * (xr - x(i-1)) * 0.5
            ELSE
               ! Segment crosses both boundaries
               ! y at centre: yc = y1 + (y2-y1)/(x2-x1)*(0.5*(xl+xr)-x1) = y1*(1-(0.5*(xl+xr)-x1)/(x2-x1)) + y2*(0.5*(xl+xr)-x1)/(x2-x1)
               ! integral: (xr-xl)*yc
               f = (0.5*(xl+xr)-x(i-1))/deltax
               w(i-1) = w(i-1) + (1.-f)*(xr-xl)
               w(i  ) = w(i  ) + (      f)*(xr-xl)
            ENDIF
         ENDIF
      ENDDO
      IF (x(1)>xl) w(1) = w(1) + (x(1)-xl)
      IF (x(n)<xr) w(n) = w(n) + (xr-x(n))
   END SUBROUTINE calculate_integral_weights
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_maths
