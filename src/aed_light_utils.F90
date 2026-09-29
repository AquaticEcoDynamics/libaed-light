!###############################################################################
!#                                                                             #
!# aed_light_utils.F90                                                         #
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
!# Created Sept 2026 (from aed_slingo, curtin_light & aed_oasim)               #
!#                                                                             #
!###############################################################################

#include "aed.h"

#define FLOAT(x) (x)

!
MODULE aed_light_utils
!-------------------------------------------------------------------------------
! Stateless routines shared by the AED light modules :
!   - solar geometry (zenith_angle) and physical constants
!   - atmosphere : ozone (vH_O3), aerosols (navy_aerosol_model),
!                  clouds (slingo), and the combined sky_spectrum
!   - Curtin direct/diffuse spectra (direct_diffuse_bands, direct_diffuse_curtin)
!   - air-water interface reflectance
!-------------------------------------------------------------------------------
   USE aed_spectra, ONLY : interp
   USE aed_light_tables

   IMPLICIT NONE

   PRIVATE
!
   PUBLIC pi, deg2rad, rad2deg, Planck, lightspeed, Avogadro
   PUBLIC zenith_angle, vH_O3, navy_aerosol_model, slingo, sky_spectrum
   PUBLIC direct_diffuse_bands, direct_diffuse_curtin, calculate_lambda_bounds
   PUBLIC reflectance
!

   AED_REAL,PARAMETER :: pi = 3.14159265358979323846
   AED_REAL,PARAMETER :: deg2rad = pi / 180.
   AED_REAL,PARAMETER :: rad2deg = 180. / pi

   AED_REAL,PARAMETER :: Planck     = 6.62606957e-34  ! Planck constant (m2 kg/s)
   AED_REAL,PARAMETER :: lightspeed = 299792458       ! Speed of light (m/s)
   AED_REAL,PARAMETER :: Avogadro   = 6.02214129e23   ! Avogadro constant (/mol)

!===============================================================================
CONTAINS

! https://journals.ametsoc.org/downloadpdf/journals/atsc/46/10/1520-0469_1989_046_1419_agpfts_2_0_co_2.pdf

!###############################################################################
SUBROUTINE slingo(mu0, LWP, r_e, nlambda, lambda, T_dcld, T_scld)
!-------------------------------------------------------------------------------
! mu0: cosine of zenith angle
! LWP: liquid water path (g m-2)
! r_e: equivalent radius of drop size distribution (um)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: nlambda
   AED_REAL,INTENT(in) :: mu0, LWP, r_e, lambda(nlambda)
   AED_REAL,INTENT(out),DIMENSION(nlambda) :: T_dcld, T_scld
!
!LOCALS
!
   AED_REAL,DIMENSION(nlambda_slingo) :: tau, omega, g
   AED_REAL,DIMENSION(nlambda_slingo) :: beta0, beta_mu0, f, U2
   AED_REAL,DIMENSION(nlambda_slingo) :: alpha1, alpha2, alpha3, alpha4, epsilon, M, E, gamma1, gamma2
   AED_REAL,DIMENSION(nlambda_slingo) :: one_minus_omega_f, gamma_denom, DIF_denom
   AED_REAL,DIMENSION(nlambda_slingo) :: T_DB, R_DIF, T_DIF, R_DIR, T_DIR
   AED_REAL,PARAMETER :: U1 = 7. / 4.

!-------------------------------------------------------------------------------
!BEGIN
   ! Slingo 1989 Eqs 1-3
   ! cloud optical depth, single scatter albedo, asymmetry parameter
   ! derived for r_e = [4.2, 16.6]
   ! NB OASIM uses a default r_e that equals the mean of 10 um (Kiehl et al., 1998 J. Clim.) and 11.8 um (Han et al., 1994 J. Clim.)
   ! NB Stephens 1978 J Atmos Sciences Eq 10 links tau directly to LWP, which allows estimation of r_e as
   !                                                                                   e.g. in Slingo 1989 section #4
   tau = LWP * (a_slingo + b_slingo / r_e)
   omega = 1. - (c_slingo + d_slingo * r_e)
   g = e_slingo + f_slingo * r_e

   ! fraction of scattered diffuse radation which is scattered into the backward hemisphere (Slingo 1989 Eq 6)
   beta0 = 3. / 7. * (1. - g)

   ! fraction of scattered direct radation which is scattered into the backward hemisphere (Slingo 1989 Eq 7)
   beta_mu0 = 0.5 - 0.75 * mu0 * g / (1. + g)

   ! fraction of scattered direct flux which emerges at zenith angles close to that of the incident beam (Slingo 1989 Eq 8)
   f = g * g

   ! reciprocals of the effective cosine for diffuse upward and downward fluxes (Slingo 1989 Eqs 9 and 10)
   ! JB 2018-09-12: Setting lower limit of U2 to 1, given that cosines cannot exceed 1 (thus its reciprocal cannot be lower than 1)
   U2 = max(1., U1 * (1. - (1. - omega) / 7. / omega / beta0))

   alpha1 = U1 * (1. - omega * (1. - beta0))
   alpha2 = U2 * omega * beta0
   alpha3 = (1 - f) * omega * beta_mu0
   alpha4 = (1 - f) * omega * (1 - beta_mu0)
   epsilon = sqrt(max(0., alpha1 * alpha1 - alpha2 * alpha2))
   M = alpha2 / (alpha1 + epsilon)
   E = exp(-epsilon * tau)
   one_minus_omega_f = 1. - omega * f
   gamma_denom = one_minus_omega_f**2 - epsilon**2 * mu0**2
   gamma1 = ( one_minus_omega_f * alpha3 - mu0 * (alpha1 * alpha3 + alpha2 * alpha4)) / gamma_denom
   gamma2 = (-one_minus_omega_f * alpha4 - mu0 * (alpha1 * alpha4 + alpha2 * alpha3)) / gamma_denom

   ! Transmissivity for the direct solar beam (Slingo 1989 Eq 20)
   T_DB = exp(-one_minus_omega_f * tau / mu0)

   ! Diffuse reflectivity for diffuse incident radiation (Slingo 1989 Eq 21)
   DIF_denom = 1. - E**2 * M**2
   R_DIF = M * (1. - E**2) / DIF_denom

   ! Diffuse transmissivity for diffuse incident radiation (Slingo 1989 Eq 22)
   T_DIF = E * (1. - M**2) / DIF_denom

   ! Diffuse reflectivity for direct incident radiation (Slingo 1989 Eq 21)
   R_DIR = max(0., -gamma2 * R_DIF - gamma1 * T_DB * T_DIF + gamma1)

   ! Diffuse transmissivity for direct incident radiation (Slingo 1989 Eq 21)
   T_DIR = min(1., -gamma2 * T_DIF - gamma1 * T_DB * R_DIF + gamma2 * T_DB)

   CALL interp(nlambda_slingo, lambda_slingo, T_DB, nlambda, lambda, T_dcld)
   CALL interp(nlambda_slingo, lambda_slingo, T_DIR, nlambda, lambda, T_scld)
END SUBROUTINE slingo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
AED_REAL FUNCTION vH_O3(long, lat, day)
!-------------------------------------------------------------------------------
! Van Heuklon (1979) function to account for seasonal and geographical variation
! in O3. See :
! Estimating Atmospheric Ozone for Solar Radiation Models - van Heuklon, T.K. (1979)
! "Solar Energy" Vol 22, pp. 63-68
! Cf Bird 1984 p 466 reports O3=0.344 atm-cm
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: long, lat, day
!
!LOCALS
   AED_REAL :: O3 ! kg m-2

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

!  AED_REAL,PARAMETER :: deg2rad = 3.14159265358979323846 / 180.
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

    ! From matm-cm to mol m-2 to kg m-2
    vH_O3 = vH_O3 * 0.4462e-3 * 48. / 1000


END FUNCTION vH_O3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
AED_REAL FUNCTION zenith_angle(days, hour, dlon, dlat)
!-------------------------------------------------------------------------------
! Spencer, J.W. (1971)
!  Fourier Series Representation of the Position of the Sun.
!  Search, 2, 162-172.
!-------------------------------------------------------------------------------
!ARGUMENTS
    AED_REAL,INTENT(in)  :: days, hour, dlon, dlat
!
!LOCALS
    AED_REAL,PARAMETER :: yrdays = 365.24
    AED_REAL           :: th0, th02, th03, sundec
    AED_REAL           :: thsun, coszen, et
    AED_REAL           :: rlon, rlat
!-------------------------------------------------------------------------------
!BEGIN
    ! from now on everything in radians
    rlon = deg2rad * dlon
    rlat = deg2rad * dlat

    ! Sun declination from Fourier expansion of Spencer (1971, Search 2:172).
    th0 = 2.*pi*days/yrdays
    th02 = 2.*th0
    th03 = 3.*th0
    sundec =  0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) &
            - 0.006758*cos(th02) + 0.000907*sin(th02)                &
            - 0.002697*cos(th03) + 0.001480*sin(th03)
    et = (0.000075+0.001868*cos(th0)-0.032077*sin(th0)-0.014615*cos(th02)-0.040849*sin(th02))*12./pi

    ! sun hour angle :
    thsun = (hour+et-12.)*15.*deg2rad + rlon

    ! cosine of the solar zenith angle :
    coszen = sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)

    zenith_angle = acos(coszen)
END FUNCTION zenith_angle


!###############################################################################
SUBROUTINE navy_aerosol_model(AM, WM, W, RH, V, costheta, alpha, beta, F_a, omega_a)
!-------------------------------------------------------------------------------
! AM: air-mass type (1 = marine aerosol-dominated, 10 continental aerosol-dominated)
! WM: wind speed averaged over past 24 h (m s-1)
! W: instantaneous wind speed (m s-1)
! RH: relative humidity (-)
! V: visibility (m)
! costheta: cosine of zenith angle
! alpha: Angstrom exponent. NB optical thickness tau = beta * lambda**(-alpha)
! beta: scale factor for optical thickness
! F_a: forward scattering probability
! omega_a: single scattering albedo
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: AM, WM, W, RH, V, costheta
   AED_REAL,INTENT(out) :: alpha, beta, F_a, omega_a
!
!LOCALS
   AED_REAL,PARAMETER :: R = 0.05
   INTEGER,PARAMETER :: nsize = 3, gridsize = 3
   AED_REAL,PARAMETER :: r_o(nsize) = (/0.03, 0.24, 2.0/)
   AED_REAL,PARAMETER :: H_a = 1000. ! Aerosol scale height (m) Gregg & Carder 1990 p1665
   AED_REAL,PARAMETER :: r_grid(gridsize) = (/0.1, 1., 10./)

   AED_REAL :: relhum, A(nsize), f, gamma, tau_a550
   AED_REAL :: dNdr(gridsize), y(gridsize), x(gridsize)
   INTEGER :: i
   AED_REAL :: c_a550
   AED_REAL :: B1, B2, B3, cos_theta_bar

!-------------------------------------------------------------------------------
!BEGIN
   ! Impose upper limit on relative humidity (RH = 1 causes division by 0 in expression for f below)
   relhum = min(0.999, RH)

   ! Amplitude functions for aerosol components (Eqs 21-23 Gregg & Carder 1990)
   ! Units are number of particles per cubic centimeter per micrometer (Gathman 1983)
   A(1) = 2000 * AM * AM
   A(2) = max(0.5, 5.866 * (WM - 2.2))
   A(3) = max(1.4e-5, 0.01527 * (W - 2.2) * R)

   ! function relating particle size to relative humidity (Eq 24 Gregg & Carder 1990)
   f = ((2. - relhum) / (6 * (1. - relhum)))**(1. / 3.)

   ! Particle density at different radii (Gregg & Carder 1990 p 1665)
   DO i = 1, gridsize
       dNdr(i) = sum(A * exp(-log(r_grid(i) / f / r_o)**2) / f)
   ENDDO

   ! Assume Junge distribution (i.e., particle density is power law of radius): dN/dr = C r^gamma
   ! Estimate gamma with least squares.
   ! Note: least-squares estimate of slope is x-y covariance / variance of x
   ! Due to the choice of r_grid (0.1, 1, 10), the sum of x_i = log10 r_i is 0.
   ! As a result, we do not need to subtract x_mean*y_mean and x_mean^2 to compute (co)variances.
   x = log10(r_grid)
   y = log10(dNdr)
   gamma = sum(x * y) / sum(x**2)

   ! Calculate Angstrom exponent from exponent of Junge distribution (Eq 26 Gregg & Carder 1990).
   ! Junge distribution: dN/d(ln r) = r dN/dr = C r^(-v)
   ! Thus, dN/dr = C r^(-v-1)
   ! The relation to gamma defined above: gamma = -v - 1. Thus, v = -gamma - 1
   ! The relationship between Junge distribution and Angstrom exponent is alpha = v - 2 (e.g. Tomasi et al. 1983)
   ! Thus, alpha = -gamma - 3
   alpha = -(gamma + 3)

   ! Estimate concentrationPARAMETER (Eqs 28, 29 Gregg & Carder 1990)
   c_a550 = 3.91 / V
   tau_a550 = c_a550 * H_a
   beta = tau_a550 * 550.**alpha

   ! AsymmetryPARAMETER (Eq 35 Gregg & Carder 1990) - called alpha in Gregg & Casey 2009
   ! Range: 0.65 (alpha >= 1.2) to 0.82 (alpha <= 0)
   ! For comparison:
   ! - Bird 1984 (Eq 15) uses cos_theta_bar = 0.64 [implied by value of F_a]
   ! - Bird & Riordan 1986 (p 91) use cos_theta_bar = 0.65
   cos_theta_bar = -0.1417 * min(max(0., alpha), 1.2) + 0.82

   ! Forward scattering probability (Eqs 31-34 Gregg & Carder 1990)
   ! NB B1-B3 are A, B, C in Gregg & Casey 2009 Eqs 3-6
   ! NB B1-B3 are AFS, BFS, ALG in Bird & Riordan 1986 Eqs 22-26
   B3 = log(1. - cos_theta_bar)
   B1 = B3 * (1.4590 + B3 * ( 0.1595 + 0.4129 * B3))
   B2 = B3 * (0.0783 + B3 * (-0.3824 - 0.5874 * B3))
   F_a = 1. - 0.5 * exp((B1 + B2 * costheta) * costheta)

   ! Single scattering albedo (Eq 36 Gregg & Carder 1990)
   ! For comparison:
   ! - Bird 1984 (Eq 15) uses omega_a = 0.928
   ! - Bird & Riordan 1986 (p 91) use omega_a = 0.945 at 400 nm for rural aerosols (AM = 10)
   ! - Shettle and Fenn 1979 (tables 28-35) report 0.982 at RH=0% to 0.9986 at RH=99% at 550 nm for their maritime aerosol model (AM=1)
   omega_a = (-0.0032 * AM + 0.972) * exp(3.06e-2 * relhum)
END SUBROUTINE navy_aerosol_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE sky_spectrum(nlambda, lambda, exter, a_o, a_v, a_u, tau_r,          &
                        theta, costheta, AM, WM, wind_speed, relhum,           &
                        visibility, airpres, cloud_cover, LWP_in, water_vapour,&
                        O3_in, direct, diffuse, alpha_a, beta_a, F_a, omega_a)
!-------------------------------------------------------------------------------
! Direct and diffuse spectral irradiance just above the water surface, from
! Bird (1984) / Bird & Riordan (1986) / Gregg & Carder (1990) clear-sky
! transmittances (incl. the Navy aerosol model), combined with Slingo (1989)
! cloud transmittances weighted by cloud cover fraction (Gregg & Casey 2009)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: nlambda
   AED_REAL,DIMENSION(nlambda),INTENT(in) :: lambda  ! wavelength (nm)
   AED_REAL,DIMENSION(nlambda),INTENT(in) :: exter   ! extraterrestrial irradiance (W/m2/nm)
   AED_REAL,DIMENSION(nlambda),INTENT(in) :: a_o, a_v, a_u ! ozone, water vapour and mixed gas absorption
   AED_REAL,DIMENSION(nlambda),INTENT(in) :: tau_r   ! Rayleigh optical thickness
   AED_REAL,INTENT(in)  :: theta, costheta           ! zenith angle (radians) and its cosine
   AED_REAL,INTENT(in)  :: AM                        ! aerosol air mass type (1: open ocean, 10: continental)
   AED_REAL,INTENT(in)  :: WM, wind_speed            ! daily mean & current wind speed @10m (m/s)
   AED_REAL,INTENT(in)  :: relhum                    ! relative humidity (-)
   AED_REAL,INTENT(in)  :: visibility                ! visibility (m)
   AED_REAL,INTENT(in)  :: airpres                   ! surface air pressure (Pa)
   AED_REAL,INTENT(in)  :: cloud_cover               ! cloud cover (fraction, 0-1)
   AED_REAL,INTENT(in)  :: LWP_in                    ! cloud liquid water content (kg m-2)
   AED_REAL,INTENT(in)  :: water_vapour              ! total precipitable water vapour (kg m-2)
   AED_REAL,INTENT(in)  :: O3_in                     ! ozone content (kg m-2)
   AED_REAL,DIMENSION(nlambda),INTENT(out) :: direct, diffuse ! (W/m2/nm)
   AED_REAL,INTENT(out) :: alpha_a, beta_a, F_a, omega_a      ! aerosol properties
!
!LOCALS
   AED_REAL,PARAMETER :: pres0 = 101300.    ! Reference air pressure in Pa (Bird 1984 p461, Bird & Riordan p89)
  !AED_REAL,PARAMETER :: ga = 0.05          ! Ground albedo (used by Bird & Riordan 1986, but discarded by Gregg & Carder 1990)
   AED_REAL,PARAMETER :: H_oz = 22.         ! Height of maximum ozone concentration (km)
   AED_REAL,PARAMETER :: r_e = (10. + 11.8) / 2 ! Equivalent radius of cloud drop size distribution (um)
                                            ! based on mean of Kiehl et al. & Han et al. (cf OASIM)
   AED_REAL :: LWP, WV, O3, M, M_prime, M_oz
   AED_REAL,DIMENSION(nlambda) :: tau_a, T_a, T_oz, T_w, T_u, T_r, T_aa, T_as
   AED_REAL,DIMENSION(nlambda) :: T_g, T_dclr, T_sclr, T_dcld, T_scld
!
!-------------------------------------------------------------------------------
!BEGIN
   LWP = LWP_in
   O3  = O3_in

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
   T_oz = exp(-a_o * O3 * M_oz)

   ! Transmittance due to water vapour absorption - should NOT use pressure-corrected airmass (see Bird and Riordan 1986)
   !  Eq 8, Bird and Riordan 1986 (Eq 7 Bird 1984 is wrong, as mentioning in B&R, p 89), Eq 19 Gregg & Carder 1990
   T_w = exp((-0.2385 * a_v * WV * M) / (1. + 20.07 * a_v * WV * M)**0.45)

   ! Transmittance due to uniformly mixed gas absorption - SHOULD use pressure corrected airmass
   !  Eq 10 Bird 1984, Eq 11 Bird and Riordan 1986, Eq 18 Gregg & Carder 1990
   !  Bird & Riordan use 118.93 rather than 118.3, but state 118.3 should be used in the future.
   T_u = exp(-1.41 * a_u * M_prime / (1. + 118.3 * a_u * M_prime)**0.45)

   ! Transmittance terms that apply for both cloudy and clear skies.
   T_g = T_oz * T_w * T_u

   ! -------------------
   ! Clear skies part
   ! -------------------

   ! Transmittance due to Rayleigh scattering (use precomputed optical thickness)
   T_r = exp(-M_prime * tau_r)

   ! Transmittance due to aerosol absorption (Eq 26 Gregg & Carder 1990)
   CALL navy_aerosol_model(AM, WM, wind_speed, relhum, visibility, costheta, alpha_a, beta_a, F_a, omega_a)
   tau_a = beta_a * lambda**(-alpha_a)
   T_a = exp(-tau_a * M)

   ! Direct transmittance
   T_dclr = T_r * T_a

   ! Separate absorption and scattering components of aerosol transmittance
   T_aa = exp(-(1. - omega_a) * tau_a * M)
   T_as = exp(-omega_a * tau_a * M)
   T_sclr = T_aa * 0.5 * (1. - T_r**0.95) + T_r**1.5 * T_aa * F_a * (1 - T_as)

   ! -------------------
   ! Cloudy skies part
   ! -------------------

   ! Transmittance due to absorption and scattering by clouds
   CALL slingo(costheta, max(0., LWP) * 1000, r_e, nlambda, lambda, T_dcld, T_scld)

   ! Diffuse and direct irradiance streams (Eqs 1, 2 Gregg & Casey 2009)
   !  These combine terms for clear and cloudy skies, weighted by cloud cover fraction
   direct  = exter * costheta * T_g * ((1. - cloud_cover) * T_dclr + cloud_cover * T_dcld)
   diffuse = exter * costheta * T_g * ((1. - cloud_cover) * T_sclr + cloud_cover * T_scld)
END SUBROUTINE sky_spectrum
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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
   rindex = FLOAT(key)
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
     x(i) = FLOAT(i)
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


END MODULE aed_light_utils



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
