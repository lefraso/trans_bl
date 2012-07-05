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
      integer i, j, k, ist,iex,pontox,pontoy
      real*8 yy,ys,xx,xs,xe,ye,sigx,sigy,x,y, z, h, yex, posy
      real*8 pi,dist,beta1,distx,a,hipot,gama
      real*8 d(imax,jmax,kphys), delta_no(ptsx,jmax,kphys)
      common/del/ delta_no

      pi     = 4.d0*datan(1.d0)
      beta1  = 130.d0 ! inclinação da falesia
      beta1  = pi*beta1/(180.d0)
      sigx   = 1.d0 ! influencia o número de pontos da gaussiana em x
      sigy   = 1.d0 ! influencia o número de pontos da gaussiana em y

      ist    = 160 ! ponto i onde começa a falesia em x
      iex    = 460 ! ponto i onde termina a falesia em x
      yex    = 2d-2 ! ponto y onde termina a falesia em y

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
        pontoy = 1
        posy = 0.d0
        do j = 1, jmax ! loop em y
          y      = dble(j-1) * dy ! posição atual em y !!! VER FORMULA DO Y COM JOSUEL
          pontox = ist
          do i = 1, imax ! loop em x
            x = dble(i-1) * dx ! posição atual em x
            gama  = beta1 - 0.5d0 * pi
            distx = dble(ist-i) * dx
            a     = (yex-y) * dtan(gama)
            dist    = ( distx - a )/dcos(gama) ! mede a distância do ponto atual à fronteira imersa (<1 dentro da fronteira)
            if (dist.gt.0.and.y.lt.dist*dsin(gama)+yex) then
              xs = dist * dcos(gama) ! distancia x do ponto (i,j) a reta
              ys = dist * dsin(gama) ! distancia y do ponto (i,j) a reta
             else 
              xs = x - (pontox-1)*dx
c             ys = y - (pontoy-1)*dy
              ys = y - posy
              if (i.gt.pontox.and.pontox.lt.iex) pontox = pontox+1
c             if (j.gt.pontoy.and.pontoy.lt.jex) pontoy = pontoy+1
              if (y.gt.posy.and.posy.lt.yex) posy = dble(j-1)*dy
            end if
            d(i,j,k) = dexp( -( xs/sigx/dx )**2 
     &                     -( ys/sigy/dy )**2 )/(sigx*sigy) ! calculo do delta
c           if (dist.lt.0.d0.and.i.lt.iex.and.j.lt.jex)  d(i,j,k) = 1.d0
            if (dist.lt.0.d0.and.i.lt.iex.and.y.lt.yex)  d(i,j,k) = 1.d0
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
        y = dble(j-1) * dy
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(1,5)x,y,z,d(i,j,kphys)
        end do
      end do
      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          y = dble(j-1) * dy
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(1,5)x,y,z,d(i,j,k)
          end do
        end do
      end do
      close (unit=1)
    5 format(1x,4d17.9)
      stop

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
