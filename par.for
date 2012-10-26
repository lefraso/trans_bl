      integer imax, ptsx, jmax, kfour, kphys, tt, i0, i1, i2, i3, i4,
     &        stpp, start, stencil, meshpx, meshdx, msh, tt_base, 
     &        my_form, np
      real*8 Re, dx, dxx, dyy, dy, dt, x0, omega, alpha, beta, dt_base,
     &       dyypdxx, dz, irf, stf, alphaf, af, bf, cf, df, Pr
      complex*16 im

c     Immersed simulations my_form = 0, Gortler simulations my_form = 1
c     Gortler simulations with heat transfer my_form = 2
      parameter ( my_form = 0 )

c     start the program from t = 0 (start=0) or from a given time t (start=1)
      parameter ( start = 0 )

c     imaginary number
      parameter ( im = (0.d0,1.d0) )

c     alpha (streamwise wavelength) and omega (frequency) omega = F(Hz)*2*pi*L/U
c     parameter ( alpha = 0.45723605d0, omega = 16.3362818d0 ) ! 130 Hz
c     parameter ( alpha = 22.6d0, omega = 10.758878952d0 ) ! 250 Hz
      parameter ( alpha = 22.6d0, omega = 23.669533695d0 ) ! 550 Hz
c     parameter ( alpha = 22.6d0, omega = 30.124861066d0 ) ! 700 Hz

c     value of beta (spanwise wavelength)
      parameter ( beta = 5.d0 )
c     parameter ( beta = 12.566370616d0 )
c     parameter ( beta = 34.9065850399d0 )

c     initial and end point of disturbance strip
      parameter ( i0 = 30, i1 = 50, i2 = 80 )

c     Reynolds Number, Prandtl Number and initial x of the domain
c     parameter ( Re = 33124.d0, Pr = 0.72d0, x0 = 1.d0 )
      parameter ( Re = 3.65d5, Pr = 0.72d0, x0 = 1.d0 )

c     number of points in y direction and delta y(dy/sqrt(re*x0))
c     parameter ( jmax = 185, dy = 8.d-4, dyy = dy * dy )
      parameter ( jmax = 121, dy = 5.d-4, dyy = dy*dy )
      parameter ( stf = 1.01d0 )
      
c     number of processing elements
      parameter ( np = 8 )

c     number of points in x direction and delta x
c     parameter ( imax = 857, ptsx = (imax + (np - 1) * 25) / np )
      parameter ( imax = 665, ptsx = (imax + (np - 1) * 25) / np )
c     parameter ( dx = 1.5d-2, dxx = dx * dx )
      parameter ( dx = 6.25d-3, dxx = dx * dx )

c     steps per period, number of time steps and time step(2*pi/omega/stpp)
      parameter ( stpp = 128, tt = 400 * stpp )
      parameter ( dt = 6.283185307179586d0 / omega / stpp )

c     for baseflow use these parameters
      parameter ( dt_base = 0.1d0 * dx )
      parameter ( tt_base = 100 * imax )
c     parameter ( tt_base = 5000 )

c     number of meshes used in the multigrid solver
      parameter ( msh = 4 )

c     parameter used in poisson subroutines
      parameter ( dyypdxx = dyy / dxx )

c     initial and end point of damping zone
      parameter ( i3 = imax-110, i4 = imax-60 )

c     number of modes in fourier and physical space
      parameter ( kfour = 11, kphys = 32 )
      parameter ( dz = 6.283185307179586d0 / (kphys * beta) )

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
