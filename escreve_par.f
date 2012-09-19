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
      character*15 nome
      integer t

      if (my_form.eq.2) then
          write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file=nome,form='unformatted')
          write(1) t
          write(1) ux,uy,uz,wx,wy,wz,th
          close (unit=1)
        else
          write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
          write(*,*) ' The results are stored in the file data_xx.bin' 
          open(1,file=nome,form='unformatted')
          write(1) t
          write(1) ux,uy,uz,wx,wy,wz
          close (unit=1)
      end if

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
      real*8 uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax)
      complex*16 uxp(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb

      if (dble(t-tt)/dble(fanal).gt.(t-tt)/fanal) return
      var = (t - tt) / fanal + 1

      write(nm,'(a,i0.2,a,i0.2)')'pert_',my_rank,'_',var
      write(*,*) ' The results are stored in the files pert_cc_XX' 
      select case(my_form)

       case(0)
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

       case default
        open(1,file=nm,form='unformatted')
        write(1) ux
        close (unit=1)

      end select

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c             end subroutines writes the results        c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
