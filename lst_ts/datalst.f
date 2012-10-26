      program geradatalst

      implicit none
      include '../par.for'

      open(2,file='datalst.dat',form='formatted')
        write(2,*) dsqrt(Re), omega/dsqrt(Re), beta/dsqrt(Re), 
     &             imax, dx
      close(unit=2)

      stop
      end
