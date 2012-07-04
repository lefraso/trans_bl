ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c         immersed boundary method                      c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine delta
      
      ! calculate the delta for the immerser boundary layer
      implicit none
      include 'par.for'
      include 'comm.par'
      integer i, j, k
      real*8 sigx, sigy, x, y, z, xc, zc, xs, zs, dist, R, h
      real*8 d(imax,jmax,kphys), delta_no(ptsx,jmax,kphys)
      common/del/ delta_no
      
      sigx = 1.d0                        ! sigma_x (gaussiana)
      sigy = 1.d0                        ! sigma_y (gaussiana)
      xc   = 1.7625                      ! x centro
      zc   = dble(kphys+1) * dz * 0.5d0  ! z centro
      R    = 0.045d0                     ! raio
      h    = 0.825d-3                    ! altura

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, imax
            d(i,j,k) = 0.d0
          end do
        end do
      end do

      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          y = dble(j-1) * dy
          do i = 2, imax
            x = x0 + dble(i-1) * dx
            if (y.le.h) then
              dist     = dsqrt( (x-xc)**2 + (z-zc)**2 )
              xs       = xc - ( R * (xc - x) ) / dist
              zs       = zc - ( R * (zc - z) ) / dist
              d(i,j,k) = dexp( - ( (x-xs) / sigx / dz)**2 -
     &                   ( (z-zs) / sigy / dz)**2 ) / (sigx * sigy)
              ! inside the circle delta = 1
              if (dist.lt.R)         d(i,j,k) = 1.d0
              if (d(i,j,k).lt.1d-14) d(i,j,k) = 0.d0
            end if
          end do
        end do
      end do

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            delta_no(i,j,k) = d(i+shift,j,k)
          end do
        end do
      end do

c     if (my_rank.eq.0) then 
c       open (1, file = 'delta.dat',status = 'unknown')
c       write(1,*) 'VARIABLES="x","y","z", "delta",'
c       write(1,*) 'ZONE I=',imax,', J=',jmax,', K=,',kphys,', F=POINT'
c       do k = 1, kphys
c         z = dble(k-1) * dz
c         do j = 1, jmax
c           y = dble(j-1)*dy
c           do i = 1, imax
c             x = x0 + dble(i-1)*dx
c             write(1,3)x,y,z,d(i,j,k)
c           end do
c         end do
c       end do
c       close (unit=1)
c   3   format(1x,3f8.4,1d17.9)
c     end if

      return
      end

********************************************
      subroutine cvirt

      ! calculates the forcing terms values
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      integer i, j, k
      real*8 fxp(ptsx,jmax,kphys),      fyp(ptsx,jmax,kphys),
     &       fzp(ptsx,jmax,kphys), delta_no(ptsx,jmax,kphys),
     &       uxp(ptsx,jmax,kphys),      uyp(ptsx,jmax,kphys),
     &       uzp(ptsx,jmax,kphys),      erro
      complex*16 fx(ptsx,jmax,kfour),    fy(ptsx,jmax,kfour),
     &           fz(ptsx,jmax,kfour)
      common/del/ delta_no
      common/force/ fx,fy,fz

      call f_to_p(uxp,ux)
      call f_to_p(uyp,uy)
      call f_to_p(uzp,uz)

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            fxp(i,j,k) =  irf * delta_no(i,j,k) * uxp(i,j,k)
            fyp(i,j,k) =  irf * delta_no(i,j,k) * uyp(i,j,k)
            fzp(i,j,k) =  irf * delta_no(i,j,k) * uzp(i,j,k)
          end do
        end do
      end do

      call p_to_f(fxp,fx)
      call p_to_f(fyp,fy)
      call p_to_f(fzp,fz)

      erro = 0.d0
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
           if (delta_no(i,j,k).eq.1.d0) erro = max(abs(uxp(i,j,k)),
     &                 abs(uyp(i,j,k)), abs(uzp(i,j,k)),erro)
          end do
        end do
      end do
      write(*,*)' Max velocity in immersed body -> ', my_rank, erro

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c         end of immersed boundary method               c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
