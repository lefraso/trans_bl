ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c              subroutines transforms binary data in asc              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 12062012

      implicit none
      include '../par.for'
      character*20 nm
      integer i, j, my_rank, inter, shift
      real*8 x, y,
     &       uxt(imax,jmax), uyt(imax,jmax), wzt(imax,jmax),
     &       ux(ptsx,jmax),  uy(ptsx,jmax),  wz(ptsx,jmax)

      ! Disturbances variables data
      inter = 2**( msh - 1 ) * ( stencil - 2 )
      do my_rank = 0, np - 1
        write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
        open(1,file=nm,form='unformatted')
        read(1)
        read(1) ux, uy, wz
        close (unit=1)
        shift = my_rank * (ptsx - inter - 1)
        do j = 1, jmax
          do i = 1, ptsx
            uxt(i+shift,j) = ux(i,j)
            uyt(i+shift,j) = uy(i,j)
            wzt(i+shift,j) = wz(i,j)
          end do
        end do
      end do

      ! writes data to spacial space to be open by tecplot
      open (3, file = 'basens.dat',status = 'unknown')
      write(3,*) 'VARIABLES="x","y","velu","vely","vortz"'
      write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

      do j = 1, jmax
        if(stf.eq.1.d0) then
         y = dble(j-1) * dy
        else
         y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
        endif
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(3,5)x,y,uxt(i,j), uyt(i,j), wzt(i,j)
        end do
      end do
    5 format(1x,2d14.6,3d17.9)

      open(1,file='basens.bin',form='unformatted')
      write(1) stf
      write(1) uxt, uyt, wzt
      close (unit=1)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c            end subroutines transforms binary data in asc            c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
