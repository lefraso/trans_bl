      subroutine init

      complex*16 uxo(ptsx,jmax,11), uyo(ptsx,jmax,11),
     &           uzo(ptsx,jmax,11), wxo(ptsx,jmax,11),
     &           wyo(ptsx,jmax,11), wzo(ptsx,jmax,11)

      if (start.eq.1) then 
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(3,file=nome,form='unformatted')
        read(3) t0
        read(3) uxo,uyo,uzo,wxo,wyo,wzo
        close(3)
        t0 = t0 + 1
      end if

      do k = 1, 11
        do i = 1, ptsx
          do j = 1, jmax
            ux(i,j,k)   = uxo(i,j,k)
            uy(i,j,k)   = uyo(i,j,k)
            uz(i,j,k)   = uzo(i,j,k)
            wx(i,j,k)   = wxo(i,j,k)
            wy(i,j,k)   = wyo(i,j,k)
            wz(i,j,k)   = wzo(i,j,k)
          end do
        end do
      end do

      return
      end
