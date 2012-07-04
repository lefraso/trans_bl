      integer imax, ptsx, jmax, kfour, kphys, tt, i0, i1, i2, i3, i4,
     &        stpp, start, stencil, meshpx, meshdx, msh, stpp_base, 
     &        tt_base 
      real*8 Re, dx, dxx, dyy, dy, dt, x0, omega, alpha, beta, dt_base,
     &       dyypdxx, dz, irf, stf, alphaf, af, bf, cf, df, beta_fs, m
      complex*16 im

c     start the program from t = 0 (start=0) or from a given time t (start=1)
      parameter ( start = 0 )

c     imaginary number
      parameter ( im = (0.d0,1.d0) )

c     alpha (streamwise wavelength) and omega (frequency)
      parameter ( alpha = 22.6d0, omega = 10.758878952d0 ) ! 250 Hz
c     parameter ( alpha = 22.6d0, omega = 23.669533695d0 ) ! 550 Hz
c     parameter ( alpha = 22.6d0, omega = 30.124861066d0 ) ! 700 Hz

c     value of beta (spanwise wavelength)
      parameter ( beta = 12.566370616d0 )

c     initial and end point of disturbance strip
      parameter ( i0 = 20, i1 = 32, i2 = 48 )

c     Reynolds number and initial x of the domain
      parameter ( Re = 3.65d5, x0 = 1.d0 )

c     number of points in y direction and delta y(dy/sqrt(re*x0))
      parameter ( jmax = 145, dy = 5.d-4, dyy = dy*dy )
      parameter ( stf = 1.01d0 )
      
c     number of points in x direction and delta x
c     parameter ( imax = 665, ptsx = 105, dx = 6.25d-3, dxx=dx*dx )
      parameter ( imax = 665, ptsx = 185, dx = 6.25d-3, dxx=dx*dx )
c     parameter ( imax = 665, ptsx = 665, dx = 6.25d-3, dxx=dx*dx )

c     steps per period, number of time steps and time step(2*pi/omega/stpp)
      parameter ( stpp = 96, tt = 15*stpp )
c     parameter ( stpp = 64, tt = 30*stpp )
c     parameter ( stpp = 64, tt = 40*stpp )
      parameter ( dt = 6.283185307179586d0/omega/stpp )

c     for baseflow use these parameters
      parameter ( stpp_base = 2500, tt_base = 100*stpp_base )
      parameter ( dt_base = 6.283185307179586d0/stpp_base )

c     number of meshes used in the multigrid solver
      parameter ( msh = 4 )

c     parameter used in poisson subroutines
      parameter ( dyypdxx = dyy / dxx )

c     initial and end point of damping zone
      parameter ( i3 = imax-100, i4 = imax-40 )
c     parameter ( i3 = imax-50, i4 = imax-20 )

c     number of modes in fourier and physical space
      parameter ( kfour = 2, kphys = 4,
     &            dz = 6.283185307179586d0/ (kphys*beta) )

c     usados para solucao de poisson paralelizado      
      parameter ( stencil = 5, meshpx  = ( stencil - 1 ) / 2 )
      parameter ( meshdx  = 12 )

c     Filter constants (Lele C.2.5)
      parameter ( alphaf = 0.40d0 )
      parameter ( af = (11.d0+10.d0*alphaf)/16.d0 )
      parameter ( bf = (15.d0+34.d0*alphaf)/64.d0 ) !/2
      parameter ( cf = (-3.d0+ 6.d0*alphaf)/32.d0 ) !/2
      parameter ( df = ( 1.d0- 2.d0*alphaf)/64.d0 ) !/2

c     parametros do falkner-skan
      parameter ( beta_fs = -0.06d0 )
      parameter ( m = beta_fs / (2.d0-beta_fs) )

c     immersed boundary method constant
      parameter ( irf = - 0.1d0*Re )
