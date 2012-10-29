      program geradatalst

      implicit none
      include '../par.for'
      integer i

      write(*,*) 'TS 2D (0) or TS 3D (1) analysis ?'
      read(*,*) i

      open(2,file='datalst.dat',form='formatted')
        write(2,*) dsqrt(Re), omega/dsqrt(Re), beta*dble(i)/dsqrt(Re),
     &             imax, dx
      close(unit=2)

      stop
      end
