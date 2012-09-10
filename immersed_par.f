ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c         immersed boundary method                      c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine delta_falesia
      
c calculate the delta for the immerser boundary layer
      implicit none
      include 'par.for'
      include 'comm.par'
      integer i, j, k
      real*8 ys, xs, sigx, sigy, x, y, z, yex, posy
      real*8 pi, dist, beta1, distx, a, gama, x_ini, x_end
      real*8 d(imax,jmax,kphys), delta_no(ptsx,jmax,kphys)
      common/del/ delta_no

      pi     = 4.d0*datan(1.d0)

      beta1  = 90.d0              ! inclinação da falesia
      beta1  = pi*beta1/180.d0    ! passa para radianos
      gama   = beta1 - 0.5d0 * pi ! angulo complentar

      sigx   = 1.d0*dy ! influencia o número de pontos da gaussiana em x
      sigy   = 1.d0*dy ! influencia o número de pontos da gaussiana em y

      x_ini  = 1.2d0          ! ponto x onde começa a falesia
      x_end  = 2.15d0         ! ponto x onde termina a falesia
      yex    = 7.73222d-2     ! altura da falesia em y

!     zerando o valor inicial de delta
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, imax
            d(i,j,k) = 0.d0
          end do
        end do
      end do

!     loop para dar valores para delta diferentes de zero
      do k = 1, kphys ! loop em z
        posy = 0.d0
        do j = 1, jmax ! loop em y
          if(stf.eq.1.d0) then ! posição atual em y
            y = dble(j-1) * dy
           else
            y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          do i = 1, imax ! loop em x
            x = dble(i-1) * dx + x0 ! posição atual em x
            distx = x_ini - x
            a     = ( yex - y ) * dtan(gama)
            dist  = ( distx - a ) * dcos(gama) ! mede a distância do ponto atual à fronteira imersa (<1 dentro da fronteira)
            if (dist.ge.0.and.y.le.dist*dsin(gama)+yex) then
              xs = dist * dcos(gama) ! distancia x do ponto (i,j) a reta
              ys = dist * dsin(gama) ! distancia y do ponto (i,j) a reta
             else 
              ys = y - posy
              if (posy.le.yex) posy = y
              xs = 0.d0
              if (x.lt.x_ini) xs = x - x_ini
              if (x.ge.x_end) xs = x - x_end
            end if
            d(i,j,k) = dexp( -(xs**2/(2.d0*sigx**2))
     &                       -(ys**2/(2.d0*sigy**2)) ) ! calculo do delta
            if (dist.le.0.d0.and.x.le.x_end.and.y.le.yex) d(i,j,k)=1.d0
            if (d(i,j,k).lt.1d-14) d(i,j,k) = 0.d0 ! para garantir que é zero depois de alguns pontos
          end do
        end do
      end do

!     distribuição do delta nos nós
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            delta_no(i,j,k) = d(i+shift,j,k)
          end do
        end do
      end do

!     para escrever o valor de delta em todos domínio
      open (1, file = 'delta.dat',status = 'unknown')
      write(1,*) 'VARIABLES="x","y","z","delta",'
      write(1,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,', F=POINT'
      z = -1.d0 * dz
      do j = 1, jmax
        if(stf.eq.1.d0) then
         y = dble(j-1) * dy
        else
         y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
        endif
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(1,5)x,y,z,d(i,j,kphys)
        end do
      end do
      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
          endif
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(1,5)x,y,z,d(i,j,k)
          end do
        end do
      end do
      close (unit=1)
    5 format(1x,4d17.9)
c     stop

      return
      end

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
      
      sigx = 1.d0
      sigy = 1.d0
      xc   = 1.7625                  ! x centro
      zc   = dble(kphys+1)*dz*0.5d0  ! z centro
      R    = 0.045d0                 ! raio
      h    = 0.825d-3                ! ALTURA

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
              xs       = xc-(R*(xc-x))/dist
              zs       = zc-(R*(zc-z))/dist
              d(i,j,k) = dexp(-((x-xs)/sigx/dz)**2-
     &                   ((z-zs)/sigy/dz)**2)/(sigx*sigy)
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
     &                 abs(uyp(i,j,k)),abs(uzp(i,j,k)),erro)
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
