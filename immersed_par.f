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
      real*8 dist, beta1, distx, a, gama, x_ini, x_end
      real*8 d(imax,jmax,kphys), delta_no(ptsx,jmax,kphys)
      common/del/ delta_no

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
          do i = 1, ptsx ! loop em x
            x = x0 + dble(i+shift-1) * dx ! posição atual em x
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
            delta_no(i,j,k) = dexp( -(xs**2/(2.d0*sigx**2))
     &                       -(ys**2/(2.d0*sigy**2)) ) ! calculo do delta
            if (dist.le.0.d0.and.x.le.x_end.and.y.le.yex) 
     &          delta_no(i,j,k)=1.d0
            if (delta_no(i,j,k).lt.1d-14) delta_no(i,j,k) = 0.d0 ! para garantir que é zero depois de alguns pontos
          end do
        end do
      end do

!     para escrever o valor de delta em todos domínio
!     write(nome,'(a,i0.2,a)')'delta_',my_rank,'.dat'
!     open (1, file = nome,status = 'unknown')
!     write(1,*) 'VARIABLES="x","y","z", "delta",'
!     write(1,*) 'ZONE I=',ptsx,', J=',jmax,', K=,',kphys,', F=POINT'
!     do k = 1, kphys
!       z = dble(k-1) * dz
!       do j = 1, jmax
!         if(stf.eq.1.d0) then
!          y = dble(j-1) * dy
!         else
!          y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
!         endif
!         do i = 1, ptsx
!           x = x0 + dble(i+shift-1)*dx
!           write(1,3)x,y,z,delta_no(i,j,k)
!         end do
!       end do
!     end do
!     close (unit=1)
!   3 format(1x,3f8.4,1d17.9)


      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine delta
      
      ! calculate the delta for the immerser boundary layer
      implicit none
      include 'par.for'
      include 'comm.par'
      integer i, j, k
      character*15 nome
      real*8 sigx, sigy, x, y, z, xc, zc, xs, zs, dist, h
      real*8 delta_no(ptsx,jmax,kphys)
      complex*16 delta_nof(ptsx,jmax,kfour)
      common/del/ delta_no
      
      sigx = 1.d0
      sigy = 1.d0
c     xc   = 2.3722d0                ! x centro zero gradient
      xc   = 1.7d0                   ! x centro zero gradient TESTE
      zc   = dble(kphys+1)*dz*0.5d0  ! z centro
c     h    = 0.1d0  * d1             ! ALTURA 10%
c     h    = 0.2d0  * d1             ! ALTURA 20%
c     h    = 0.3d0  * d1             ! ALTURA 30%
c     h    = 0.4d0  * d1             ! ALTURA 40%
      h    = 0.5d0  * d1             ! ALTURA 50%

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            delta_no(i,j,k) = 0.d0
          end do
        end do
      end do

      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          if (stf.eq.1.d0) then ! posição atual em y
            y = dble(j-1) * dy
           else
            y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          do i = 1, ptsx
            x = x0 + dble(i+shift-1) * dx
            if (y.le.h) then
              dist     = dsqrt( (x-xc)**2 + (z-zc)**2 )
              xs       = xc-(R*(xc-x))/dist
              zs       = zc-(R*(zc-z))/dist
              delta_no(i,j,k) = dexp(-((x-xs)/sigx/dz)**2-
     &                   ((z-zs)/sigy/dz)**2)/(sigx*sigy)
              if (dist.lt.R)                delta_no(i,j,k) = 1.d0
              if (delta_no(i,j,k).lt.1d-14) delta_no(i,j,k) = 0.d0
            end if
          end do
        end do
      end do

      call p_to_f(delta_no,delta_nof)
      call gibbs_filter(delta_nof)
      call f_to_p(delta_no,delta_nof)

      if (my_rank.eq.2) then
      write(nome,'(a,i0.2,a)')'delta_',my_rank,'.dat'
      open (1, file = nome,status = 'unknown')
      write(1,*) 'VARIABLES="x","y","z", "delta",'
      write(1,*) 'ZONE I=',ptsx,', J=',jmax,', K=,',kphys,', F=POINT'
      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
          endif
          do i = 1, ptsx
            x = x0 + dble(i+shift-1)*dx
            write(1,3)x,y,z,delta_no(i,j,k)
          end do
        end do
      end do
      close (unit=1)
      end if
    3 format(1x,3f8.4,1d17.9)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gibbs_filter(fc)
      
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      integer i, j, k
      complex*16 fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            fc(i,j,k) = gibbs(k) * fc(i,j,k)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gibbs_function

      implicit none
      include 'par.for'
      include 'comm.var'
      integer k
      real*8 theta, gibbs0, ag, bg, cg, dg, alphag, betag

      alphag = -0.50d0
      betag  = 0.d0
      ag = (  11.d0 + 10.d0 * alphag - 10.d0 * betag) / 16.d0
      bg = (  15.d0 + 34.d0 * alphag + 30.d0 * betag) / 32.d0
      cg = (-  3.d0 +  6.d0 * alphag + 26.d0 * betag) / 16.d0
      dg = (   1.d0 -  2.d0 * alphag +  2.d0 * betag) / 32.d0

      do k = 1, kfour
!       Leandro e João Henrique
        theta    = pi * dble(k-1) / dble(kfour-1)
        gibbs(k) = ( ag + bg * dcos(theta) + 
     &               cg * dcos(2.d0*theta) + 
     &               dg * dcos(3.d0*theta) )/
     &             ( 1.d0 + 2.d0 * alphag * dcos(theta) 
     &                    + 2.d0 * betag  * dcos(2.d0*theta)) !lele

        gibbs(1) = 1.d0
        gibbs(k) = gibbs(k) * (0.5d0 * (1.d0 + dcos(theta)))**2.0
!       Leandro e João Henrique

!       Leandro e Larissa
c       theta    = dble(k-1) / dble(kfour-1)
c       gibbs0   = 1.d0 + (( - 6.d0 * theta + 15.d0) * theta - 10.d0)
c    &           * theta**3
c       theta    = pi * dble(k) / dble(kfour)
c       gibbs(k) = gibbs0 * (0.5d0 * (1.d0 + dcos(theta)))**0.6
!       Leandro e Larissa

!       Lanczos
c       theta    = pi * dble(k-1) / dble(kfour-1)
c       gibbs(k) = dsin(theta) / theta
c       if (theta .eq. 0.d0) gibbs(k) = 1.d0
!       Lanczos

!       Exponencial
        theta    = pi * dble(k-1) / dble(kfour-1)
        gibbs(k) = dexp(-theta**2)
!       Exponencial

!       Raised cosine
c       theta    = pi * dble(k-1) / dble(kfour-1)
c       gibbs(k) = 0.5d0 * (1.d0 + dcos(theta))
!       Raised cosine
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

      call gibbs_filter(fx)
      call gibbs_filter(fy)
      call gibbs_filter(fz)

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
