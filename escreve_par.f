ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c                 subroutines writes the results        c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine escreve(t)

      ! write the results in binary form
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character c1
      character*2 c2
      character*11 nome
      integer t

      if (my_rank.le.9) then
        write (c1,'(I1)'),my_rank
        nome='data_0'//c1//'.bin'
       else
        write (c2,'(I2)'),my_rank
        nome = 'data_'//c2//'.bin'
      end if

      write(*,*) ' The results are stored in the file data_xx.bin' 
      open(1,file=nome,form='unformatted')
      write(1) t
      write(1) ux,uy,uz,wx,wy,wz
      close (unit=1)

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine escreve2(t,fanal)

      ! write the results for fourier analysis in time
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      character*20 nm
      integer i, j, k, t, fanal, var
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax), dya
      complex*16 uxp(ptsx,jmax,kfour)
      common/blas/ uxb,uyb,wzb

      if (dble(t-tt)/dble(fanal).gt.(t-tt)/fanal) return
      var = (t - tt) / fanal + 1

      write(nm,'(a,i0.2,a,i0.2)')'pert_',my_rank,'_',var
      write(*,*) ' The results are stored in the files pert_cc_XX' 
      do j = 1, jmax
        do i = 1, ptsx
          uxp(i,j,1) = ux(i,j,1) - uxb(i,j)
        end do
      end do
      do k = 2, kfour
        do j = 1, jmax
          do i = 1, ptsx
            uxp(i,j,k) = ux(i,j,k)
          end do
        end do
      end do

      open(1,file=nm,form='unformatted')
      write(1) uxp
      close (unit=1)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c             end subroutines writes the results        c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
