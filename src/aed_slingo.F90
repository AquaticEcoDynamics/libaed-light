!

#include "aed.h"

MODULE aed_slingo

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
   ! NB Stephens 1978 J Atmos Sciences Eq 10 links tau directly to LWP, which allows estimation of r_e as e.g. in Slingo 1989 section #4
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


END MODULE aed_slingo
