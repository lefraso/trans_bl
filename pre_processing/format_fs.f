ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c              subroutines transforms binary data in asc              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 04092012

      implicit none
      include '../par.for'
      character*20 nm
      integer i, j, my_rank, inter, shift
      real*8 x, y,
     &       ux(imax,jmax), uy(imax,jmax), 
     &       wz(imax,jmax), th(imax,jmax)

      if (my_form.eq.2) then
        open(1,file='base_fs.bin',form='unformatted')
        read(1) ux, uy, wz, th
        close (unit=1)
        
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'base_fs.dat',status = 'unknown')
        write(3,*) 'VARIABLES="x","y","velu","vely","vortz","theta"'
        write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
          endif
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(3,6) x, y, ux(i,j), uy(i,j), wz(i,j), th(i,j)
          end do
        end do
    6   format(1x,2d14.6,4d17.9)
        
       else
        open(1,file='base_fs.bin',form='unformatted')
        read(1) ux, uy, wz
        close (unit=1)
        
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'base_fs.dat',status = 'unknown')
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
            write(3,5) x, y, ux(i,j), uy(i,j), wz(i,j)
          end do
        end do
    5   format(1x,2d14.6,3d17.9)
        
      end if

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c            end subroutines transforms binary data in asc            c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
