#include "aed.h"

! I refer to the contents of file Direct_Diffuse_Notes.pdf whose title is: "Technical note on direct and diffuse surface spectral
! vector irradiance above air-sea interface for Cockburn Sound"
!
! NOTE:  The subroutine direct_diffuse_curtin is not written in a way that can be overloaded.

MODULE curtin_light

   USE aed_maths

   IMPLICIT NONE

   PRIVATE cloud_proportion

   PUBLIC direct_diffuse_curtin, maxswflux_statistical

!  MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


#include "oasim.inc"

!===============================================================================
CONTAINS



SUBROUTINE direct_diffuse_curtin(SWFlux, MaxSWFlux, theta, month, met, direct, diffuse, uv, par)
!-----------------------------------------------------------------------------------------------
! INPUT:
!   SWFlux [W/m2]
!   MaxSWFlux [W/m2], if negative, use an estimate provided by FUNCTION maxswflux_statistical(theta)
!   theta [degrees], solar zenith angle in range 5-95 (domain never sees theta<5, above 90 non-zero diffuse)
!   month [1-12], month of year - not range checked
!   met ['B' or 'W'], meteorological model, B=BARRA, W=WRF
! OUTPUT:
!   direct [W/m2/nm]
!   diffuse [W/m2/nm]
!   uv [W/m2]
!   par [W/m2]
!-----------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in)  :: SWFlux, MaxSWFlux, theta
   INTEGER,INTENT(in)          :: month
   CHARACTER(len=1),INTENT(in) :: met
   AED_REAL,INTENT(out) :: direct(nlambda_5nm_astm), diffuse(nlambda_5nm_astm), uv, par

!LOCALS
   AED_REAL, parameter :: MINSWFLUX = 0.1 ! [W/m2]
   ! Coefficients for solar zenith angle dependence from EcoLight (RADTRAN-X) simulations
   AED_REAL, DIMENSION(5, nlambda_curtin) :: ce, ci, cex, cix ! e is for direct, i is for diffuse
   AED_REAL, DIMENSION(nlambda_curtin)    :: se, si, sex, six ! spectra 300nm-1um, x is for cloudy 
   AED_REAL, DIMENSION(nlambda_5nm_astm)  :: se2, si2         ! spectra 300nm-4um
   AED_REAL, DIMENSION(nlambda_5nm_astm)  :: es, taug
   INTEGER                                :: i, i1, i2
   AED_REAL,DIMENSION(ntheta_curtin)      :: muv
   TYPE(AED_REAL)                         :: maxflux, mu, swflux_ratio, cloud_index, prop, flux_curtin, ratio, const, tranc

!------------------------------------------------------------------------------------
!BEGIN

   ! MaxSWFlux < 0 means use a function built from forecast/reanalysis timeseries of SWFlux over the Cockburn Sound domain
   if (MaxSWFlux < 0) then
      maxflux = maxswflux_statistical(theta, month, met)
   else
      maxflux = MaxSWFlux
   end if

   ! Problem if maxflux or SWFlux tiny, set all outputs to 0 and return
   if ((maxflux < MINSWFLUX) .or. (SWFlux < MINSWFLUX)) then
     direct = 0.0
     diffuse = 0.0
     uv = 0.0
     par = 0.0
     return
   end if

   ! Test for cloud-free using ratio SWFlux / MaxSWFlux
   swflux_ratio = SWFlux / maxflux

   ! mu is the cosine of solar zenith angle, muv a vector of powers of mu
   mu = cos(min(85.0,theta)*deg2rad) ! the regression coefficients not applicable above 85 degrees
   muv(1) = 1.0
   do i = 2,ntheta_curtin
     muv(i) = muv(i-1) * mu
   end do 
   
   ! se and si are direct and diffuse that, at mesh nodes, should closely match EcoLight simulations

   ! Need clear-sky values even for cloudy skies to use in estimating 1um-4um tail
   ce = reshape(coeff_curtin(:,1,:,1), (/ ntheta_curtin, nlambda_curtin /) )
   ci = reshape(coeff_curtin(:,2,:,1), (/ ntheta_curtin, nlambda_curtin /) )
   se = max(0.0, matmul(muv, ce))
   si = max(0.0, matmul(muv, ci))
   ! se=direct & si=diffuse are pretty accurate estimates of the clear-sky EcoLight simulated values 300nm-1um

   ! Only need to do this for cloudy skies
   if (swflux_ratio < 0.98) then
      cloud_index = cloud_proportion(swflux_ratio) * 10.0 - 1.0
      i1 = floor(cloud_index)
      i2 = ceiling(cloud_index)
      prop = modulo(cloud_index, 1.0)
      ce = reshape(coeff_curtin(:,1,:,i1), (/ 5, nlambda_curtin /) )
      ci = reshape(coeff_curtin(:,2,:,i1), (/ 5, nlambda_curtin /) )
      cex = reshape(coeff_curtin(:,1,:,i2), (/ 5, nlambda_curtin /) )
      cix = reshape(coeff_curtin(:,2,:,i2), (/ 5, nlambda_curtin /) )
      sex = max(0.0, (1.0-prop) * matmul(muv, ce) + prop * matmul(muv, cex))
      six = max(0.0, (1.0-prop) * matmul(muv, ci) + prop * matmul(muv, cix))
   end if
   ! sex=direct & six=diffuse are pretty accurate estimates of the cloudy-sky EcoLight simulated values 300nm-1um

   ! Add 1um-4um tail to clear-sky direct
   es = max(1e-37, es_5nm_astm)              ! to prevent a log(0)
   taug = log(et_5nm_astm/es)                ! above 700nm, taug is effectively a gaseous (absorption) optical depth
   se2 = mu * et_5nm_astm * exp(-taug/mu)    ! taug and a slant path can give a good estimate of direct irradiance
   ratio = se(140) / se2(140)                ! a-priori knowledge, index 140 is 997.5nm
   se2 = max(0.0, se2 * ratio)               ! ratio should be close to 1, i.e. a small adjustment to 1um-4um tail
   se2(1:140) = se(1:140)                    ! swap in the EcoLight-sourced clear-direct for <1um

   ! Add 1um-4um tail to clear-sky diffuse
   si2 = ed_5nm_astm - es_5nm_astm           ! above 700nm, total less (direct + circumsolar) is good approx to diffuse
   ratio = si(140) / si2(140)                ! a-priori knowledge, index 140 is 997.5nm
   si2 = max(0.0, si2 * ratio)               ! ratio can be <<1 for large solzen BUT clear diffuse is tiny at >1um
   si2(1:140) = si(1:140)                    ! swap in the EcoLight-sourced clear-diffuse for <1um

   ! If flux ratio so indicates, the 1um-4um tail needs to be for cloudy skies
   if (swflux_ratio < 0.98) then
   
     ! Add 1um-4um cloudy-sky direct
     tranc = sex(140) / se(140)              ! slant path cloud transmittance tends to constant at 1um
     se2 = tranc * se2                       ! modify clear direct to get cloudy direct
     se2(1:140) = sex(1:140)                 ! swap in the EcoLight-sourced cloudy-direct for <1um

     ! Add in 1um-4um cloudy-sky diffuse
     const = sum(sex) / sum(six)             ! Ecolight says direct/diffuse a constant for given cloud fraction
     si2 = se2 / const

   end if

   ! If sun lower than horizon, nix the direct component and attribute all to diffuse
   if (theta >= 90) se2 = 0.0

   ! use some a-priori knowledge; lambda spacing is 5nm
   flux_curtin = 5.0 * sum(se2 + si2)
   ratio = SWFlux / flux_curtin

   ! scale direct and diffuse to integrate to input SWFlux
   direct = ratio * se2
   diffuse = ratio * si2

   ! more a-priori knowledge; indices 1, 20, 21, 80 centered at 302.5, 397.5, 402.5, 697.5
   uv = 5.0 * sum(direct(1:20) + diffuse(1:20))
   par = 5.0 * sum(direct(21:80) + diffuse(21:80))

END SUBROUTINE direct_diffuse_curtin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



AED_REAL FUNCTION cloud_proportion(swflux_ratio)
!```````````````````````````````````````````````
! This function is not executed if swflux_ratio > 0.98
! swflux_ratio could be in range 0.0 to 0.98, but swflux_ratio<0.25 is effectively forced to 0.25
! NOTE:  What the above is saying is that under full cloud cover, the flux ratio is at minimum 0.25 in EcoLight.
!        Look to see if that is borne out in the reanalysis product - frankly it does not look like it. Figure C.
! cloud_proportion will return in range 0.338 to 0.999
! REF: Figure F.
!ARGUMENTS
   AED_REAL,INTENT(in) :: swflux_ratio
!LOCALS
   TYPE(AED_REAL)      :: r
!BEGIN
   r = swflux_ratio
   if (r<0.25) r=0.25
   cloud_proportion = -4.2754 + 30.4163*r + 5.1592*exp(r) + 0.3864*exp(4*r) -39.1796*tan(r)
END FUNCTION cloud_proportion



AED_REAL FUNCTION maxswflux_statistical(theta, month, met)
!`````````````````````````````````````````````````````````!
AED_REAL, parameter :: doy(12) = (/ 15., 45., 74., 105., 135., 166., 196., 227., 258., 288., 319., 349. /)
!ARGUMENTS
   AED_REAL,INTENT(in) :: theta            ! solar zenith angle in degrees
   INTEGER,INTENT(in)  :: month            ! needs to be in range 1-12, not checked
   CHARACTER(len=1),INTENT(in) :: met      ! 'B' for BARRA, else 'W' for WRF is assumed
!LOCALS
   TYPE(AED_REAL)              :: t, mu, f
!BEGIN
   if (theta<5) then
     t = 5 ! you should throw a warning, for the Cockburn Sound domain this ought not occur
   else if (theta>100) then
     t = 100 ! no warning needed, just night-time, at worst will evaluate -ve and get set to zero
   else
     t = theta
   end if
   mu = cos(t*deg2rad)
   if (met=='B') then
     maxswflux_statistical = 44.0329 + (493.04213 + (1226.9471 - 632.4993*mu)*mu)*mu  ! BARRA
   else
     maxswflux_statistical = 0.868905 + (902.6874 + (492.08104 - 239.7974*mu)*mu)*mu  ! WRF
   end if
   f = (1 + 0.034 * cos(2.*pi*doy(month)/365.)) / 1.03287 ! scaling factor for month of year
   maxswflux_statistical = f * maxswflux_statistical
   if (maxswflux_statistical < 0) maxswflux_statistical = 0
END FUNCTION maxswflux_statistical



END MODULE curtin_light
