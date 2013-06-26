      program integ 10072002

      ! this program reads the output of fs_bound.f and calculates the
      ! values of displacement thickness (alfa1), momentum thickness (alfa2)
      ! and calculates the H12=alfa1/alfa2
      implicit none
      include '../par.for'
      integer i, j, jmaxh
      real*8 delta1(imax), delta2(imax), H12(imax), fc1, fc2, umax, 
     &       dya, ux(imax,jmax), uy(imax,jmax), wz(imax,jmax)
      
      dya = dy * dsqrt(Re)

      open(1,file='baseflow.bin',form='unformatted')
      read(1) ux, uy, wz
      close(unit=1)
      
      fc1 = 0.d0
      fc2 = 0.d0
      jmaxh = jmax

      do i = 1, imax
        umax = ux(i,jmax)
        fc1  = 3.d0 * (umax-ux(i,1)) / 8.d0
        fc1  = fc1 + 7.d0 * (umax-ux(i,2)) / 6.d0
        fc1  = fc1 + 23.d0 * (umax-ux(i,3)) / 24.d0
        fc2  = 3.d0 * ux(i,1) * (umax-ux(i,1)) / 8.d0
        fc2  = fc2 + 7.d0 * ux(i,2) * (umax-ux(i,2)) / 6.d0
        fc2  = fc2 + 23.d0 * ux(i,3) * (umax-ux(i,3)) / 24.d0
        do j = 4, jmaxh - 3
          fc1 = fc1 + umax - ux(i,j)
          fc2 = fc2 + ux(i,j) * ( umax-ux(i,j) )
        end do
        fc1       = fc1 + 23.d0 * (umax - ux(i,jmaxh-2)) / 24.d0
        fc1       = fc1 +  7.d0 * (umax - ux(i,jmaxh-1)) / 6.d0
        fc1       = fc1 +  3.d0 * (umax - ux(i,jmaxh)) / 8.d0
        fc2       = fc2 + 23.d0 * 
     &              ux(i,jmaxh-2) * (umax-ux(i,jmaxh-2)) / 24.d0
        fc2       = fc2 +  7.d0 * 
     &              ux(i,jmaxh-1) * (umax-ux(i,jmaxh-1)) / 6.d0
        fc2       = fc2 +  3.d0 * 
     &              ux(i,jmaxh) * (umax-ux(i,jmaxh)) / 8.d0
        delta1(i) = dya * fc1
        delta2(i) = dya * fc2 / umax
        H12(i)    = delta1(i) / delta2(i)
      end do

      open (1, file = 'H12.dat',status = 'unknown')
      write(1,*) 'VARIABLES="i","delta1","delta2","H12"'
      write(1,*) 'ZONE I=',imax,',  F=POINT'
      write(*,*) ' The results are stored in the file H12.dat' 
      do i = 1, imax
        write(1,5) dble(i), delta1(i), delta2(i), H12(i)
      end do
      close (unit=1)
 
    5 format(1x,1f8.2,3d25.17)
      return
      end
