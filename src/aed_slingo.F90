!

#include "aed.h"

MODULE aed_slingo

   USE aed_maths

PUBLIC

AED_REAL,PARAMETER :: a_slingo(24) = (/ &
      3.094000e-02, 2.944000e-02, 3.308000e-02,  2.801000e-02, 2.668000e-02,   &
      2.698000e-02, 2.672000e-02, 2.838000e-02,  2.831000e-02, 2.895000e-02,   &
      3.115000e-02, 2.650000e-02, 2.622000e-02,  2.497000e-02, 2.632000e-02,   &
      2.589000e-02, 2.551000e-02, 2.463000e-02,  2.237000e-02, 1.970000e-02,   &
      1.850000e-02, 1.579000e-02, 1.950000e-02, -1.023000e-02                  &
   /)

AED_REAL,PARAMETER :: b_slingo(24) = (/ &
      1.252000e+00, 1.270000e+00, 1.246000e+00, 1.293000e+00, 1.307000e+00,    &
      1.315000e+00, 1.320000e+00, 1.300000e+00, 1.317000e+00, 1.315000e+00,    &
      1.244000e+00, 1.349000e+00, 1.362000e+00, 1.376000e+00, 1.365000e+00,    &
      1.385000e+00, 1.401000e+00, 1.420000e+00, 1.452000e+00, 1.501000e+00,    &
      1.556000e+00, 1.611000e+00, 1.540000e+00, 1.933000e+00 &
   /)

AED_REAL,PARAMETER :: c_slingo(24) = (/ &
     7.900000e-07, -6.500000e-07, -3.000000e-07,  1.000000e-06,  0.000000e+00, &
     1.000000e-06,  0.000000e+00,  0.000000e+00, -1.200000e-06, -1.200000e-07, &
    -2.700000e-07,  2.300000e-06,  3.300000e-06,  9.800000e-06, -4.600000e-05, &
    -2.800000e-05,  6.200000e-05,  2.400000e-04,  1.200000e-04,  1.200000e-03, &
     1.900000e-04,  1.230000e-01,  4.490000e-01,  2.500000e-02 &
   /)

AED_REAL,PARAMETER :: d_slingo(24) = (/ &
      3.690000e-07, 4.330000e-07, 2.360000e-07, 0.000000e+00, 0.000000e+00,    &
      0.000000e+00, 0.000000e+00, 0.000000e+00, 4.000000e-07, 4.400000e-07,    &
      1.400000e-06, 1.700000e-06, 2.800000e-06, 2.100000e-05, 5.000000e-05,    &
      8.000000e-05, 2.600000e-04, 8.560000e-04, 6.670000e-04, 2.160000e-03,    &
      2.540000e-03, 9.350000e-03, 1.540000e-03, 1.220000e-02 &
   /)

AED_REAL,PARAMETER :: e_slingo(24) = (/ &
      8.440000e-01, 8.410000e-01, 8.390000e-01, 8.360000e-01, 8.400000e-01,    &
      8.200000e-01, 8.280000e-01, 8.250000e-01, 8.280000e-01, 8.180000e-01,    &
      8.040000e-01, 8.090000e-01, 8.060000e-01, 7.830000e-01, 7.840000e-01,    &
      7.800000e-01, 7.730000e-01, 7.540000e-01, 7.490000e-01, 7.400000e-01,    &
      7.690000e-01, 8.510000e-01, 8.310000e-01, 7.260000e-01 &
   /)

AED_REAL,PARAMETER :: f_slingo(24) = (/ &
      1.558000e-03, 1.680000e-03, 1.946000e-03, 2.153000e-03, 1.881000e-03,    &
      3.004000e-03, 2.467000e-03, 2.776000e-03, 2.492000e-03, 2.989000e-03,    &
      3.520000e-03, 3.387000e-03, 3.355000e-03, 5.035000e-03, 4.745000e-03,    &
      4.989000e-03, 5.405000e-03, 6.555000e-03, 6.931000e-03, 7.469000e-03,    &
      5.171000e-03, 2.814000e-03, 6.102000e-03, 6.652000e-03 &
   /)

AED_REAL,PARAMETER :: lambda_max_slingo(24) = (/ &
      3.000000e+02, 3.300000e+02, 3.600000e+02, 4.000000e+02, 4.400000e+02,    &
      4.800000e+02, 5.200000e+02, 5.700000e+02, 6.400000e+02, 6.900000e+02,    &
      7.500000e+02, 7.800000e+02, 8.700000e+02, 1.000000e+03, 1.100000e+03,    &
      1.190000e+03, 1.280000e+03, 1.530000e+03, 1.640000e+03, 2.130000e+03,    &
      2.380000e+03, 2.910000e+03, 3.420000e+03, 4.000000e+03 &
   /)

AED_REAL,PARAMETER :: lambda_min_slingo(24) = (/ &
      2.500000e+02, 3.000000e+02, 3.300000e+02, 3.600000e+02, 4.000000e+02,    &
      4.400000e+02, 4.800000e+02, 5.200000e+02, 5.700000e+02, 6.400000e+02,    &
      6.900000e+02, 7.500000e+02, 7.800000e+02, 8.700000e+02, 1.000000e+03,    &
      1.100000e+03, 1.190000e+03, 1.280000e+03, 1.530000e+03, 1.640000e+03,    &
      2.130000e+03, 2.380000e+03, 2.910000e+03, 3.420000e+03 &
   /)

AED_REAL,PARAMETER :: lambda_slingo(24) = (/ &
      2.750000e+02, 3.150000e+02, 3.450000e+02, 3.800000e+02, 4.200000e+02,    &
      4.600000e+02, 5.000000e+02, 5.450000e+02, 6.050000e+02, 6.650000e+02,    &
      7.200000e+02, 7.650000e+02, 8.250000e+02, 9.350000e+02, 1.050000e+03,    &
      1.145000e+03, 1.235000e+03, 1.405000e+03, 1.585000e+03, 1.885000e+03,    &
      2.255000e+03, 2.645000e+03, 3.165000e+03, 3.710000e+03 &
   /)

INTEGER,PARAMETER :: nlambda_slingo = 24


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
    AED_REAL           :: thsun, coszen
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

    ! sun hour angle :
    thsun = (hour-12.)*15.*deg2rad + rlon

    ! cosine of the solar zenith angle :
    coszen = sin(rlat)*sin(sundec)+cos(rlat)*cos(sundec)*cos(thsun)

    zenith_angle = acos(coszen)
END FUNCTION zenith_angle
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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



END MODULE aed_slingo

