!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 module constants

  implicit none

  ! General constants
  integer, parameter      :: my_form       = 2
  integer, parameter      :: imax          = 665
  real(kind=8), parameter :: dx            = 6.25d-3
  integer, parameter      :: jmax          = 121
  real(kind=8), parameter :: dy0           = 5.d-4
  real(kind=8), parameter :: stf           = 1.01d0
  real(kind=8), parameter :: x0            = 1.d0
  real(kind=8), parameter :: Re            = 3.65d5
  real(kind=8), parameter :: Pr            = 0.72d0
  integer, parameter      :: lvls          = 4

  ! Falkner-Skan constants
  real(kind=8), parameter :: tol_fs        = 1.0d-13
  real(kind=8), parameter :: gtol_fs       = 1.0d-23
  real(kind=8), parameter :: etapp_guess   = 1.1d0
  real(kind=8), parameter :: eta_zero      = .5d0
  real(kind=8), parameter :: eta_variation = 0.01d0
  real(kind=8), parameter :: deta          = 1.d-4
  real(kind=8)            :: beta_fs(imax)

  contains

 end module constants

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 module derivative_commom_variables

  use constants

  ! first derivative coefficients
  real(kind=8)            :: fp_fd_coef(7)
  real(kind=8)            :: sp_fd_coef(9)
  real(kind=8)            :: cp_fd_coef(8)
  real(kind=8)            :: pp_fd_coef(9)
  real(kind=8)            :: lp_fd_coef(7)
 
  ! second derivative coefficients
  real(kind=8)            :: fp_sd_coef(8)
  real(kind=8)            :: sp_sd_coef(10)
  real(kind=8)            :: cp_sd_coef(8)
  real(kind=8)            :: pp_sd_coef(10)
  real(kind=8)            :: lp_sd_coef(8)
 
  ! Poisson coefficients
  real(kind=8)            :: sp_poi_coef(7)
  real(kind=8)            :: cp_poi_coef(8)
  real(kind=8)            :: pp_poi_coef(4)
  real(kind=8)            :: lp_poi_coef(5)

  ! Vorticity coefficients
  real(kind=8)            :: w_at_w_coef(9)
 
  real(kind=8)            :: dwydy_coef(8)
 
  real(kind=8)            ::  a1x(imax),  b1x(imax),  c1x(imax)
  real(kind=8)            ::  a2x(imax),  b2x(imax),  c2x(imax)
  real(kind=8)            ::  a1y(jmax),  b1y(jmax),  c1y(jmax)
  real(kind=8)            ::  a2y(jmax),  b2y(jmax),  c2y(jmax)
  real(kind=8)            ::  a1fy(jmax), b1fy(jmax), c1fy(jmax)
 
  real(kind=8)            :: sp_integ_coef(7)
  real(kind=8)            :: mp_integ_coef(8)
  real(kind=8)            :: pp_integ_coef(7)
  real(kind=8)            :: lp_integ_coef(7)

 end module derivative_commom_variables

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
