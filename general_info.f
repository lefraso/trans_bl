
      program general_info
      
       implicit none
       include 'par.for'
       integer ny
       real*8 fac, delta, x_ad, x_dim, yl, outboundheigh

       ! number of boundary layer thickness at outflow boundary
       fac   = 2.5d0

       ! verification of number of points adopted for the boundary 
       ! layer thickness at inflow boundary
       x_ad  = 1.d0                              ! non dimensional x at inflow
       x_dim = x_ad * L_1                        ! dimensional x at inflow
       delta = 5.d0 * dsqrt( x_dim * N_1 / U_1 ) ! delta calculation based in blasius BL
       yl    = delta                             ! dimensional BL height
       yl    = yl * dsqrt(fac_y) / L_1           ! non dimensional BL height
       if ( stf .eq. 1.d0 ) then
        ny = yl / dy
       else
        ny = log((yl * (stf - 1.d0)/ dy + 1.d0)) / log(stf) + 1
       endif
       write(*,*) 'Number of points used in the BL at inflow:', ny
       ! verification of number of points adopted for the boundary 
       ! layer thickness at inflow boundary

       ! verification of number of points the should be adopted
       ! at outflow boundary
       x_ad  = 1.d0 + dble(imax-1) * dx          ! non dimensional x at outflow
       x_dim = x_ad * L_1                        ! dimensional x at outflow
       delta = 5.d0 * dsqrt( x_dim * N_1 / U_1 ) ! delta calculation based in blasius BL
       yl    = fac * delta                       ! fac x dimensional BL height
       yl    = yl * dsqrt(fac_y) / L_1           ! fac x non dimensional BL height
       if ( stf .eq. 1.d0 ) then
        ny = yl / dy
       else
        ny = log((yl * (stf - 1.d0)/ dy + 1.d0)) / log(stf) + 1
       endif
       if ( mod(ny-1, 2**(msh-1)) .ne. 0.d0 ) then ! calculation to garantee the correct
        ny = ny + ( np - mod(ny-1, 2**(msh-1)) )   ! use of multigrid method
       endif
       write(*,*)
       write(*,*) 'You should use at least ', ny, 'points in the normal
     &             direction, to have ', fac,'BL heights at outflow'
       write(*,*) 'The first dy is:', dy
       write(*,*) 'The last dy with ',ny,' points is:', dy*stf**(ny-1) 
       ! verification of number of points the should be adopted
       ! at outflow boundary

       ! domain height compared to boundary layer thickness at outflow boundary
       if ( stf.eq.1.d0 ) then
        outboundheigh = dble(jmax-1) * dy
       else
c       outboundheigh = dy * (stf**(jmax-1)-1.d0)/(stf-1.d0)
       endif
       write(*,*) 'You are adopting ', outboundheigh / (yl / fac), 
     &            'BL at outflow'
       ! domain height compared to boundary layer thickness at outflow boundary

       ! CFL and dt calculation
       write(*,*)
       write(*,*) 'dtmin based in viscous terms is ', Re/2.d0*
     &            (dx**2*dy**2)/(dx**2+dy**2)
       write(*,*) 'dtmin based in convective terms is ', dt / dx
       write(*,*) 'You are using a delta t = ', dt
       ! CFL and dt calculation
      
      end program general_info
