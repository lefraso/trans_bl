      integer imax, ptsx, jmax, kfour, kphys, tt, i0, i1, i2, i3, i4,
     &        stpp, start, stencil, meshpx, meshdx, msh, tt_base, 
     &        my_form, np, cc
      real*8 Re, dx, dxx, dyy, dy, dt, x0, omega, alpha, beta, dt_base,
     &       dyypdxx, dz, irf, stf, alphaf, af, bf, cf, df, Pr, U_1,
     &       L_1, N_1, Raio, R, pi, d1, fac_y, C_s, C_w, delta_les
      complex*16 im

c     Immersed simulations my_form = 0
c     Gortler simulations my_form = 1
c     Gortler simulations with heat transfer my_form = 2
c     LES simulations my_form = 4
      parameter ( my_form = 4 )

c     start the program from t = 0 (start=0) or from a given time t (start=1)
      parameter ( start = 0 )

c     imaginary number
      parameter ( im = (0.d0,1.d0) )

c     initial and end point of disturbance strip
      parameter ( i0 = 10, i1 = 15, i2 = 50 )

c     constant to calculate ue      
      parameter ( U_1 = 27.935d0 )         
c     parameter ( U_1 = 5.0d0 )         

c     lenght scale parameter L
      parameter ( L_1 = 0.18d0 )
c     parameter ( L_1 = 0.1d0 )

c     dynamic viscosity (nu = 1.56d-5 )
      parameter ( N_1 = 1.56d-5 )
c     parameter ( N_1 = 1.5095d-5 )

c     Reynolds Number, Prandtl Number and initial x of the domain
      parameter ( Re = U_1*L_1/N_1, Pr = 0.72d0, x0 = 1.d0 )

c     adimensionalization parameter in the y direction (Re or 1.d0)
      parameter (fac_y = 1.d0)
c     parameter (fac_y = Re)

c     displacement thickness at the roughness center
      parameter ( d1 = 0.784048d-3 * dsqrt(fac_y) / L_1 )   ! zero gradient at x = 2.3722

c     roughness radio (dimensional, nondimensional)
      parameter ( Raio = 1.d-2, R = Raio / L_1 )

c     pi number
      parameter ( pi = 3.141592653589793d0 )

c     alpha (streamwise wavelength) and omega (frequency) omega = F(Hz)*2*pi*L/U
      parameter ( alpha = 0.45723605d0 )
      parameter ( omega = 549.316406d0*2.d0*pi*L_1/U_1 )

c     value of beta (spanwise wavelength) (2*pi/lambda_z)
c     parameter ( beta = (2.d0 * pi) / (18.d-3 / L_1) ) ! 18mm
      parameter ( beta = (2.d0 * pi) / (8.d0 * R) )     ! lambda_z = 0.4 (4 times the roughness diameter)

c     number of points in y direction and delta y(dy/sqrt(re*x0))
c     parameter ( jmax = 257, dy = 5.0d-4*dsqrt(fac_y), dyy = dy * dy)
      parameter ( jmax = 177, dy = 1.8d-4*dsqrt(fac_y), dyy = dy * dy)
c     parameter ( jmax = 137, dy = 1.8d-4*dsqrt(fac_y), dyy = dy * dy)
c     parameter ( jmax = 169, dy = 1.8d-4*dsqrt(fac_y), dyy = dy * dy)
      parameter ( stf = 1.01d0 )
      
c     number of processing elements
      parameter ( np = 8 )

c     number of points in x direction and delta x
c     parameter ( imax = 1177, ptsx = (imax + (np - 1) * 25) / np )
      parameter ( imax = 473, ptsx = (imax + (np - 1) * 25) / np )
      parameter ( dx = 1.5d-2, dxx = dx * dx )
c     parameter ( dx = 3.125d-3, dxx = dx * dx )

c     initial and end point of damping zone
      parameter ( i3 = imax-110, i4 = imax-60 )

c     steps per period, number of time steps and time step(2*pi/omega/stpp)
      parameter ( stpp = 128, tt = 35 * stpp )
c     parameter ( stpp = 256, tt = 35 * stpp )
      parameter ( dt = (2.d0 * pi) / omega / stpp )

c     for baseflow use these parameters
      parameter ( dt_base = 0.01d0 * dx )
c     parameter ( tt_base = 500 * imax  )
      parameter ( tt_base = 3 * imax / dt )

c     number of meshes used in the multigrid solver
      parameter ( msh = 4 )

c     parameter used in poisson subroutines
      parameter ( dyypdxx = dyy / dxx )

c     number of modes in fourier and physical space
      parameter ( kfour = 2, kphys = 4 )
c     parameter ( kfour = 21, kphys = 64 )
      parameter ( dz = (2.d0 * pi) / (kphys * beta) )

c     usados para solucao de poisson paralelizado
      parameter ( stencil = 5, meshpx  = ( stencil - 1 ) / 2 )
      parameter ( meshdx  = 12 )

c     Filter constants (Lele C.2.5)
      parameter ( alphaf = 0.48d0 )
      parameter ( af = (  11.d0 + 10.d0 * alphaf) / 16.d0 )
      parameter ( bf = (  15.d0 + 34.d0 * alphaf) / 64.d0 ) !/2
      parameter ( cf = (-  3.d0 +  6.d0 * alphaf) / 32.d0 ) !/2
      parameter ( df = (   1.d0 -  2.d0 * alphaf) / 64.d0 ) !/2

c     immersed boundary method constant
      parameter ( irf = - 0.1d0 * Re )

c     Smagorinsky and Wale LES constants
      parameter ( C_s = 0.032d0, C_w = 0.104d0 )

c     Delta for the LES model
      parameter ( delta_les = (dx * dy * dz)**(1.d0/3.d0) )

C     curvature case 
      parameter ( cc = 1 )
