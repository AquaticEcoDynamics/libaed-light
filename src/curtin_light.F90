#include "aed.h"

! I refer to the contents of file Direct_Diffuse_Notes.pdf whose title is: 
!"Technical note on direct and diffuse surface spectral
!   vector irradiance above air-sea interface for Cockburn Sound"
!
! NOTE:  The subroutine direct_diffuse_curtin is not written in a way that can be overloaded.

MODULE curtin_light

   USE aed_maths

   IMPLICIT NONE

   PRIVATE cloud_proportion, maxswflux_statistical

   PUBLIC direct_diffuse_bands, direct_diffuse_curtin, calculate_lambda_bounds

!  MODULE GLOBALS
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs


#include "oasim.inc"

!===============================================================================
CONTAINS

SUBROUTINE direct_diffuse_bands(SWFlux, thetarad, yearday, met, nlambda, lambda_bounds, direct, diffuse)
!-------------------------------------------------------------------------------------------------------
! IMPORTANT:  See NOTE appended to this source file for limitations and distortions.
! INPUT:
!   SWFlux [W/m2]
!   thetarad [radians], solar zenith angle
!   yearday [1-366], day of year
!   met ['B' or 'W'], meteorological model, B=BARRA, W=WRF
!   nlambda
!   lambda_bounds
! OUTPUT:
!   direct [W/m2/nm]
!   diffuse [W/m2/nm]
!-----------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in)  :: SWFlux, thetarad, yearday
   CHARACTER(len=1),INTENT(in) :: met
   INTEGER,INTENT(in)          :: nlambda
   AED_REAL,INTENT(in)         :: lambda_bounds(nlambda+1)
   AED_REAL,INTENT(out) :: direct(nlambda), diffuse(nlambda)

!LOCALS
   INTEGER                                :: i, n, i1, i2
   INTEGER, PARAMETER, DIMENSION(*)       :: key = (/(i,i=1,nlambda_5nm_astm)/)
   AED_REAL, DIMENSION(nlambda_5nm_astm)  :: direct_5nm, diffuse_5nm, rindex, value
   TYPE(AED_REAL)                         :: thetadeg, uv, par, r1, r2, w
   AED_REAL, DIMENSION(nlambda)           :: rindex_beg, rindex_end
 
!------------------------------------------------------------------------------------
!BEGIN
   ! theta is thetarad in degrees, because that is what I coded for in SUBROUTINE direct_diffuse_curtin()
   thetadeg = thetarad / deg2rad
   ! get the 5nm spectra for this theta, yearday and for the cloudiness predicted by SWFlux
   call direct_diffuse_curtin(SWFlux, -1.0, thetadeg, yearday, met, direct_5nm, diffuse_5nm, uv, par)
   ! aggregate to the nlambda output bands defined by the nlambda+1 lambda_bounds
   ! lambda_5nm_astm are the MID-POINTS of the (contiguous) 5nm bands, lambda_bounds are the EDGES of
   !   the OASIM or custom multi-spectral (also contiguous) bands.  How are you going to effieciently
   !   do the aggregation?  value and width dimensioned to 40+2 and zeroed in each loop
   !   of iband = 1, nlambda.  Why this dimension - the max size of any band aggregate 200nm
   rindex = float(key)
   call interp(nlambda_5nm_astm, lambda_5nm_astm - 2.5, rindex, nlambda, lambda_bounds(1:nlambda), rindex_beg)
   call interp(nlambda_5nm_astm, lambda_5nm_astm + 2.5, rindex, nlambda, lambda_bounds(2:nlambda+1), rindex_end)
   ! Initialise the multi-spectral direct & diffuse output vectors
   direct = 0.0
   diffuse = 0.0
   do i = 1, nlambda
      r1 = ceiling(rindex_beg(i))
      r2 = floor(rindex_end(i))
      i1 = min(nlambda_5nm_astm, int(r1))
      i2 = max(1, int(r2))
      ! The 5nm channels fully within the current multi-spectral band
      if (r2.ge.r1) then
         direct(i) = direct(i) + sum(direct_5nm(i1:i2))*5.0
         diffuse(i) = diffuse(i) + sum(diffuse_5nm(i1:i2))*5.0
      end if
      ! Ought not be necessary but better to avoid out-of-bounds possibility than a small error in irradiance
      i1 = max(i1-1,1)
      i2 = min(i2+1,nlambda_5nm_astm)
      ! The partial 5nm channel on the shortwave side of the current multi-spectral band
      if (rindex_beg(i).lt.r1) then
         w = (r1 - rindex_beg(i))*5.0
         direct(i) = direct(i) + direct_5nm(i1)*w
         diffuse(i) = diffuse(i) + diffuse_5nm(i1)*w
      end if
      ! The partial 5nm channel on the longwave side of the current multi-spectral band
      if (rindex_end(i).gt.r2) then
         w = (rindex_end(i) - r2)*5.0
         direct(i) = direct(i) + direct_5nm(i2)*w
         diffuse(i) = diffuse(i) + diffuse_5nm(i2)*w
      end if
      w = lambda_bounds(i+1)-lambda_bounds(i)
      direct(i) = direct(i)/w
      diffuse(i) = diffuse(i)/w

      if(lambda_bounds(i+1)<283.0) then
        direct(i) = 0.0
        diffuse(i) = 0.0
      endif 
   end do
END SUBROUTINE direct_diffuse_bands
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


SUBROUTINE direct_diffuse_curtin(SWFlux, MaxSWFlux, theta, day_of_year, met, direct, diffuse, uv, par)
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
   AED_REAL,INTENT(in)  :: SWFlux, MaxSWFlux, theta, day_of_year
  !INTEGER,INTENT(in)          :: month
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
   AED_REAL                               :: maxflux, mu, swflux_ratio, cloud_index, prop, flux_curtin, ratio, const, tranc

!------------------------------------------------------------------------------------
!BEGIN

   ! MaxSWFlux < 0 means use a function built from forecast/reanalysis timeseries of SWFlux over the Cockburn Sound domain
   if (MaxSWFlux < 0) then
      maxflux = maxswflux_statistical(theta, day_of_year, met)
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
   taug =cos(37.*deg2rad)*log(et_5nm_astm/es)! above 700nm, taug is effectively a gaseous (absorption) optical depth
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



SUBROUTINE calculate_lambda_bounds(nlambda, lambda, lambda_bounds)
!-----------------------------------------------------------------
! Assumption is that bands are contiguous, so we have n+1 band boundaries to service nlambda bands
! The SUBROUTINE calculate_integral_weights() in aed_maths.F90 does something a bit funky at the
!   UV end of spectrum.  I may have to do something funky too so that the numerical integration
!   over UV matches my expectation.
! N.B. This actually gets called ONLY for lambda_method=1 i.e. OASIM wavelength grid in oasim.inc
!   which means that I can safely hard-code end bounds to concur with expectations in aed_oasim.F90
!   See NOTE appended to this source file for limitations and distortions.
! INPUT:
!   nlambda is number of wavelengths
!   lambda are the wavelengths
! OUTPUT:
!   lambda_bounds are the bounds
!-------------------------------------------------------------------------------
AED_REAL, parameter :: xl = 300.0      ! This are the limits of the Curtin light
AED_REAL, parameter :: xr = 4000.0     !   estimates, so just hard-code them in.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER, INTENT(in)  :: nlambda
   AED_REAL,INTENT(in)  :: lambda(nlambda)
   AED_REAL,INTENT(out) :: lambda_bounds(nlambda+1)
!
!LOCALS
   INTEGER  :: i
   AED_REAL,DIMENSION(:),ALLOCATABLE :: x, tx
!-------------------------------------------------------------------------------
!BEGIN
   lambda_bounds = -999.9
   ALLOCATE(x(nlambda), tx(nlambda))
   do i = 1,nlambda
     x(i) = float(i)
     tx(i) = x(i) - 0.5
   end do
   call interp(nlambda, x, lambda, nlambda+1, tx, lambda_bounds(1:nlambda))
   ! finagling ...
   lambda_bounds(1) = xl
   lambda_bounds(2) = 304.166667
   lambda_bounds(nlambda+1) = xr
   !x = lambda_bounds(2:nlambda+1) - lambda_bounds(1:nlambda)
   !do i = 1, nlambda + 1
   !  if (i.le.nlambda) then
   !     write (0,*) i, lambda_bounds(i), x(i)
   !   else
   !     write (0,*) i, lambda_bounds(i)
   !   end if
   !end do
END SUBROUTINE calculate_lambda_bounds
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
   AED_REAL      :: r
!BEGIN
   r = swflux_ratio
   if (r<0.25) r=0.25
   cloud_proportion = -4.2754 + 30.4163*r + 5.1592*exp(r) + 0.3864*exp(4*r) -39.1796*tan(r)
END FUNCTION cloud_proportion



AED_REAL FUNCTION maxswflux_statistical(theta, doy, met)
!`````````````````````````````````````````````````````````!
!AED_REAL, parameter :: doy(12) = (/ 15., 45., 74., 105., 135., 166., 196., 227., 258., 288., 319., 349. /)
!ARGUMENTS
   AED_REAL,INTENT(in) :: theta, doy       ! solar zenith angle in degrees
   !INTEGER,INTENT(in)  :: month              ! needs to be in range 1-365, not checked
   CHARACTER(len=1),INTENT(in) :: met      ! 'B' for BARRA, else 'W' for WRF is assumed
!LOCALS
   AED_REAL              :: t, mu, f
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
    !maxswflux_statistical = 44.0329 + (493.04213 + (1226.9471 - 632.4993*mu)*mu)*mu  ! BARRA
     maxswflux_statistical = 42.81 + (487.54 + (1240.40 - 676.85*mu)*mu)*mu  ! BARRA
   else
    !maxswflux_statistical = 0.868905 + (902.6874 + (492.08104 - 239.7974*mu)*mu)*mu  ! WRF
     maxswflux_statistical = 3.145 + (873.26 + (571.63 - 329.36*mu)*mu)*mu   ! WRF
   end if
   !f = (1 + 0.034 * cos(2.*pi*doy(month)/365.)) / 1.03287 ! scaling factor for month of year
   !f = (1 + 0.034 * cos(2.*pi*doy/365.)) / 1.03287 ! scaling factor for month of year
   f = 1 + 0.034 * cos(2.*pi*doy/365.)   ! scaling factor for day of year
   maxswflux_statistical = f * maxswflux_statistical
   if (maxswflux_statistical < 0) maxswflux_statistical = 0
END FUNCTION maxswflux_statistical


END MODULE curtin_light


! NOTE
! ====
! These are the weights you get from running calculate_integral_weights()
! over bands PAR, SWR and UV for the lambda=lambda_oasim read from oasim.inc.
! The weights are basically the bandwidths in nm, because the spectra are in
! units W/m2/nm.  But what is going on at the UV end of spectrum?  I need to
! ensure that I finagle the lambda_bounds that support the band weights I know
! I am going to get so that PAR, SWR and UV W/m2 come out about right in the
! computations in existing codebase.  JED 18-Oct-2023
!
! Addenda
! -------
! I have a fundamental problem arising from the boundaries for UV (300-400) and PAR (400-700)
! splitting 25nm aggregated bands. e.g band oasim_lambda(5)=400nm.  The value in this band,
! direct or diffuse, is the mean value of 5nm channels 19,20,21,22 plus half of channel 18
! half of channel 23.  So what should get apportioned to PAR is the aggregate of channels 21,22
! and half of 23.  Instead it gets half (see its weight under par_weights) of the aggregate
! of 19-22 + half 18 + half 23 BECAUSE band oasim_lambda(5) has ALREADY been aggregated
! BEFORE the method of computing PAR from oasim bands and weights is applied. I cannot
! un-mix this problem without choosing different multi-spectral band boundaries so that
! UV and PAR boundaries do not split multi-spectral bands.   JED 19-Oct-2023 
! 
!
!   i    lambda  par_weights   swr_weights    uv_weights    lambda_bounds
!   -   -------  -----------   -----------    ----------    -------------
!   1   250.000         0.00        4.1666        4.1666   300.00  304.16 
!   2   325.000         0.00        33.333        33.333   304.16  337.50 
!   3   350.000         0.00        25.000        25.000   337.50  362.50 
!   4   375.000         0.00        25.000        25.000   362.50  387.50 
!   5   400.000         12.5        25.000        12.500   387.50  412.50 
!   6   425.000         25.0        25.000        0.0000   412.50  437.50 
!   7   450.000         25.0        25.000        0.0000   437.50  462.50 
!   8   475.000         25.0        25.000        0.0000   462.50  487.50 
!   9   500.000         25.0        25.000        0.0000   487.50  512.50 
!  10   525.000         25.0        25.000        0.0000   512.50  537.50 
!  11   550.000         25.0        25.000        0.0000   537.50  562.50 
!  12   575.000         25.0        25.000        0.0000   562.50  587.50 
!  13   600.000         25.0        25.000        0.0000   587.50  612.50 
!  14   625.000         25.0        25.000        0.0000   612.50  637.50 
!  15   650.000         25.0        25.000        0.0000   637.50  662.50 
!  16   675.000         25.0        25.000        0.0000   662.50  687.50 
!  17   700.000         12.5        25.000        0.0000   687.50  712.50 
!  18   725.000         0.00        37.500        0.0000   712.50  750.00 
!  19   775.000         0.00        62.500        0.0000   750.00  812.50 
!  20   850.000         0.00        87.500        0.0000   812.50  900.00 
!  21   950.000         0.00        100.00        0.0000   900.00  1000.0 
!  22   1050.00         0.00        100.00        0.0000   1000.0  1100.0 
!  23   1150.00         0.00        100.00        0.0000   1100.0  1200.0 
!  24   1250.00         0.00        100.00        0.0000   1200.0  1300.0 
!  25   1350.00         0.00        100.00        0.0000   1300.0  1400.0 
!  26   1450.00         0.00        100.00        0.0000   1400.0  1500.0 
!  27   1550.00         0.00        100.00        0.0000   1500.0  1600.0 
!  28   1650.00         0.00        100.00        0.0000   1600.0  1700.0 
!  29   1750.00         0.00        125.00        0.0000   1700.0  1825.0 
!  30   1900.00         0.00        225.00        0.0000   1825.0  2050.0 
!  31   2200.00         0.00        500.00        0.0000   2050.0  2550.0 
!  32   2900.00         0.00        750.00        0.0000   2550.0  3300.0 
!  33   3700.00         0.00        700.00        0.0000   3300.0  4000.0 
