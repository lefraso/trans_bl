ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c       program that calculates Q for visualization     c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program isoq 19012009

      implicit none
      include 'par.for'
      include 'comm.fourier'
      integer i, j, k, var
      real*8    duxdxp(imax,jmax,kphys), duxdyp(imax,jmax,kphys)
      real*8    duydxp(imax,jmax,kphys), duydyp(imax,jmax,kphys)
      real*8    duzdxp(imax,jmax,kphys), duzdyp(imax,jmax,kphys)
      real*8    duzdzp(imax,jmax,kphys), duxdzp(imax,jmax,kphys)
      real*8    duydzp(imax,jmax,kphys),      Q(imax,jmax,kphys)
      complex*16 duxdx(imax,jmax,kfour),  duxdy(imax,jmax,kfour)
      complex*16 duydx(imax,jmax,kfour),  duydy(imax,jmax,kfour)
      complex*16 duzdx(imax,jmax,kfour),  duzdy(imax,jmax,kfour)
      complex*16 duzdz(imax,jmax,kfour),  duxdz(imax,jmax,kfour)
      complex*16 duydz(imax,jmax,kfour)
      complex*16   uxt(imax,jmax,kfour),    uyt(imax,jmax,kfour)
      complex*16   uzt(imax,jmax,kfour)
      common/iso/ uxt,uyt,uzt

      do var = 1, 16
        call initval(var)
        call derivs_kt
       
        call derxt(duxdx,uxt)
        call derxt(duydx,uyt)
        call derxt(duzdx,uzt)
        call deryt(duxdy,uxt)
        call deryfvt(duydy,uyt)
        call deryt(duzdy,uzt)
       
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, imax
              duxdz(i,j,k) = v_kb(k) * uxt(i,j,k)
              duydz(i,j,k) = v_kb(k) * uyt(i,j,k)
              duzdz(i,j,k) = v_kb(k) * uzt(i,j,k)
            end do
          end do
        end do
       
        call f_to_pt(duxdxp,duxdx)
        call f_to_pt(duydxp,duydx)
        call f_to_pt(duzdxp,duzdx)
        call f_to_pt(duxdyp,duxdy)
        call f_to_pt(duydyp,duydy)
        call f_to_pt(duzdyp,duzdy)
        call f_to_pt(duxdzp,duxdz)
        call f_to_pt(duydzp,duydz)
        call f_to_pt(duzdzp,duzdz)
       
        do k = 1, kphys
          do j = 1, jmax
            do i = 1, imax
                Q(i,j,k) = - 0.5d0 * (  duxdxp(i,j,k) * duxdxp(i,j,k) +
     &                       duydyp(i,j,k) * duydyp(i,j,k) +
     &                       duzdzp(i,j,k) * duzdzp(i,j,k) +
     &                       2.d0 * ( duxdyp(i,j,k) * duydxp(i,j,k) +
     &                       duxdzp(i,j,k) * duzdxp(i,j,k) +
     &                       duydzp(i,j,k) * duzdyp(i,j,k) )  )
            end do
          end do
        end do
       
        call escreveq(Q, var)
      end do

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initval(var)

      implicit none
      include 'par.for'
      character*15 nome
      integer i, j, t, my_rank, k, inter, shift, var
      complex*16  ux(ptsx,jmax,kfour),  uy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour), uxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour), uzt(imax,jmax,kfour)
      common/iso/ uxt,uyt,uzt
      
      ! variables data
      inter = 2**( msh - 1 ) * ( stencil - 2 )
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a,i0.2)')'pert_',my_rank,'_',var
        open(1,file=nome,form='unformatted')
        read(1) ux,uy,uz
        close (unit=1)
        shift = my_rank * (ptsx - inter - 1)
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, ptsx
              uxt(i+shift,j,k) = ux(i,j,k)
              uyt(i+shift,j,k) = uy(i,j,k)
              uzt(i+shift,j,k) = uz(i,j,k)
            end do
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine escreveq(Q, var)

      ! write the results in fourier modes
      implicit none
      include 'par.for'
      character*15 nome
      integer i, j, k, var
      real*8 x, y, z
      real*8 Q(imax,jmax,kphys)

      ! writes data to spacial space to be open by tecplot
      write(nome,'(a,i0.2,a)')'isoq_',var,'.dat'
      open (3, file = nome,status = 'unknown')
      write(3,*) 'VARIABLES="x","y","z","Q"'
      write(3,*) 'ZONE I=',imax/2+1,',J=',jmax,',K=',kphys+1,',F=POINT'

      z = - dz
      do j = 1, jmax
        y = dble(j-1) * dy
        do i = 1, imax, 2
          x = x0 + dble(i-1) * dx
          if (abs(Q(i,j,kphys)).lt.1e-15) Q(i,j,kphys) = 0.d0
          write(3,5)x,y,z,Q(i,j,kphys)
        end do
      end do

      do k = 1, kphys
        z = dble(k-1) * dz
        do j = 1, jmax
          y = dble(j-1) * dy
          do i = 1, imax, 2
            x = x0 + dble(i-1) * dx
            if (abs(Q(i,j,k)).lt.1e-15) Q(i,j,k) = 0.d0
            write(3,5)x,y,z,Q(i,j,k)
          end do
        end do
      end do
      close (unit=3)

    5 format(1x,3d14.6,1d17.9)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs_kt

      ! calculate the tri-diagonal LHS for the derivative calculation
      implicit none
      include 'par.for'
      real*8  a1x(imax),  b1x(imax),  c1x(imax)
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/ a1x,b1x,c1x
      common/der1y/ a1y,b1y,c1y
      common/der1fv/ a1fv,b1fv,c1fv

      call coeft(a1x,b1x,c1x,imax)
      call coeft(a1y,b1y,c1y,jmax)
      call coeffvt(a1fv,b1fv,c1fv,jmax)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derxt(ddx,fc)

      ! first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1x(imax), b1x(imax), c1x(imax)
      complex*16 rhs(imax), u(imax)
      complex*16 fc(imax,jmax,kfour), ddx(imax,jmax,kfour)
      common/der1x/ a1x,b1x,c1x

      do k = 1, kfour
        do j = 1, jmax
          call rhsxt(fc,rhs,j,k)
          call tridxt(rhs,a1x,b1x,c1x,u)
          do i = 1, imax
            ddx(i,j,k) = u(i)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryt(ddy,fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 rhs(jmax), u(jmax)
      complex*16 fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1y/ a1y,b1y,c1y

      do k = 1, kfour
        do i = 1, imax
          call rhsyt(fc,rhs,i,k)
          call tridyt(rhs,a1y,b1y,c1y,u)
          do j = 1, jmax
            ddy(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryfvt(ddy,fc)

      ! first derivatives calculation in y direction for uy
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      complex*16 rhs(jmax), u(jmax)
      complex*16 fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1fv/ a1fv,b1fv,c1fv

      do k = 1, kfour
        do i = 1, imax
          call rhsyfvt(fc,rhs,i,k)
          call tridyt(rhs,a1fv,b1fv,c1fv,u)
          do j = 1, jmax
            ddy(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsxt(fc,rhs,j,k)

      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(imax), fc(imax,jmax,kfour)

      rhs(1)      = ( -74.d0*fc(1,j,k) + 16.d0*fc(2,j,k) +
     &                 72.d0*fc(3,j,k) - 16.d0*fc(4,j,k) +
     &                  2.d0*fc(5,j,k))/(24.d0*dx)

      rhs(2)      = ( -406.d0*fc(1,j,k) - 300.d0*fc(2,j,k) +
     &                 760.d0*fc(3,j,k) - 80.d0*fc(4,j,k) +
     &                  30.d0*fc(5,j,k) - 4.d0*fc(6,j,k))/(120.d0*dx)

      do i = 3, imax - 2
        rhs(i)    = ( fc(i+2,j,k) + 28.d0 * (fc(i+1,j,k)-fc(i-1,j,k))-
     &                fc(i-2,j,k) )/(12.d0*dx)
      end do

      rhs(imax-1) = ( - 406.d0*fc(imax,j,k) - 300.d0*fc(imax-1,j,k)+
     &                  760.d0*fc(imax-2,j,k) - 80.d0*fc(imax-3,j,k)+
     &                   30.d0*fc(imax-4,j,k) - 4.d0*fc(imax-5,j,k) )/
     &              (-120.d0*dx)

      rhs(imax)   = ( -74.d0*fc(imax,j,k)+16.d0*fc(imax-1,j,k)+
     &                 72.d0*fc(imax-2,j,k)-16.d0*fc(imax-3,j,k)+
     &                  2.d0*fc(imax-4,j,k) )/(-24.d0*dx)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsyt(fc,rhs,i,k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = ( -74.d0*fc(i,1,k) + 16.d0*fc(i,2,k) +
     &            72.d0*fc(i,3,k) - 16.d0*fc(i,4,k) +
     &             2.d0*fc(i,5,k) )/(24.d0*dy)

      rhs(2) = ( -406.d0*fc(i,1,k) - 300.d0*fc(i,2,k) +
     &            760.d0*fc(i,3,k) -  80.d0*fc(i,4,k) +
     &             30.d0*fc(i,5,k) -   4.d0*fc(i,6,k) )/(120.d0*dy)

      do j = 3, jmax - 2
        rhs(j)  = ( fc(i,j+2,k) + 28.d0*( fc(i,j+1,k) - fc(i,j-1,k) )-
     &              fc(i,j-2,k) )/(12.d0*dy)
      end do

      rhs(jmax-1) = ( -406.d0*fc(i,jmax,k)   - 300.d0*fc(i,jmax-1,k) +
     &                 760.d0*fc(i,jmax-2,k) -  80.d0*fc(i,jmax-3,k) +
     &                  30.d0*fc(i,jmax-4,k) -   4.d0*fc(i,jmax-5,k) )/
     &              (-120.d0*dy)

      rhs(jmax)   = ( -74.d0*fc(i,jmax,k)   + 16.d0*fc(i,jmax-1,k) +
     &                 72.d0*fc(i,jmax-2,k) - 16.d0*fc(i,jmax-3,k) +
     &                  2.d0*fc(i,jmax-4,k) )/(-24.d0*dy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsyfvt(fc,rhs,i,k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = 0.d0

      rhs(2) = ( -406.d0*fc(i,1,k) - 300.d0*fc(i,2,k)+
     &            760.d0*fc(i,3,k) -  80.d0*fc(i,4,k)+
     &             30.d0*fc(i,5,k) -   4.d0*fc(i,6,k) )/(120.d0*dy)

      do j = 3, jmax - 2
        rhs(j) = ( fc(i,j+2,k) + 28.d0 * ( fc(i,j+1,k)-fc(i,j-1,k) )-
     &             fc(i,j-2,k) )/(12.d0*dy)
      end do

      rhs(jmax-1) = ( -406.d0*fc(i,jmax,k)   - 300.d0*fc(i,jmax-1,k)+
     &                 760.d0*fc(i,jmax-2,k) -  80.d0*fc(i,jmax-3,k)+
     &                 30.d0*fc(i,jmax-4,k)  -   4.d0*fc(i,jmax-5,k) )/
     &              (-120.d0*dy)

      rhs(jmax)   = ( -74.d0*fc(i,jmax,k)   + 16.d0*fc(i,jmax-1,k) +
     &                 72.d0*fc(i,jmax-2,k) - 16.d0*fc(i,jmax-3,k) +
     &                  2.d0*fc(i,jmax-4,k) )/(-24.d0*dy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridxt(r,a,b,c,u)
      
      ! solves the tridiagonal matrix for derivatives in x direction
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax), b(imax), c(imax), gam(imax), bet
      complex*16 r(imax), u(imax)
      
      bet  = b(1)
      u(1) = r(1) / bet
      do i = 2, imax
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
c       if (bet.eq.0.) pause 'tridag failed'
        u(i)   = ( r(i) - a(i) * u(i-1) ) / bet
      end do
      do i = imax - 1, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridyt(r,a,b,c,u)

      ! solves the tridiagonal matrix for derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 r(jmax), u(jmax)

      bet  = b(1)
      u(1) = r(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
c       if (bet.eq.0.) pause 'tridag failed'
        u(j)   = ( r(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeft(a,b,c,lmax)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      integer l, lmax
      real*8 a(lmax), b(lmax), c(lmax)
      
      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do l = 3, lmax - 2
        a(l)    = 1.d0
        b(l)    = 3.d0
        c(l)    = 1.d0
      end do

      a(lmax-1) = 2.d0
      b(lmax-1) = 6.d0
      c(lmax-1) = 1.d0

      a(lmax)   = 4.d0
      b(lmax)   = 1.d0
      c(lmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeffvt(a,b,c,lmax)
      
      implicit none
      integer l, lmax
      real*8 a(lmax), b(lmax), c(lmax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do l = 3, lmax - 2
        a(l)    = 1.d0
        b(l)    = 3.d0
        c(l)    = 1.d0
      end do

      a(lmax-1) = 2.d0
      b(lmax-1) = 6.d0
      c(lmax-1) = 1.d0

      a(lmax)   = 4.d0
      b(lmax)   = 1.d0
      c(lmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f_to_pt(datap,dataf)

      ! transform from Fourier space to Physical space
      ! dataf-> Fourier space (in)
      ! datap -> Physical space (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp
      real*8 datap(imax,jmax,kphys)
      complex*16 dataf(imax,jmax,kfour)

      theta = -6.28318530717959d0/dble(kphys)
      c1    = 0.5d0
      wpr   = -2.0d0*dsin(0.5d0*theta)**2
      wpi   = dsin(theta)
      
      do j = 1, jmax
        do i = 1, imax
          do k = 0, kfour-1
            datap(i,j,2*k+1) = dreal(dataf(i,j,k+1))
            datap(i,j,2*k+2) = dimag(dataf(i,j,k+1))
          end do
          do k = 2 * kfour + 1, kphys ! if 2*kfour > 2/3 kphys = alias
            datap(i,j,k) = 0.d0
          end do
          datap(i,j,1) = datap(i,j,1)*2.d0
          wr           = 1.0d0 + wpr
          wi           = wpi
          n2p3         = kphys + 3
          do p = 2, kphys/4
            p1            = 2*p-1
            p2            = p1+1
            p3            = n2p3-p2
            p4            = p3+1
            wrs           = wr
            wis           = wi
            h1r           = c1*(datap(i,j,p1)+datap(i,j,p3))
            h1i           = c1*(datap(i,j,p2)-datap(i,j,p4))
            h2r           = -c1*(datap(i,j,p2)+datap(i,j,p4))
            h2i           = c1*(datap(i,j,p1)-datap(i,j,p3))
            datap(i,j,p1) = h1r+wrs*h2r-wis*h2i
            datap(i,j,p2) = h1i+wrs*h2i+wis*h2r
            datap(i,j,p3) = h1r-wrs*h2r+wis*h2i
            datap(i,j,p4) = -h1i+wrs*h2i+wis*h2r
            wtemp         = wr
            wr            = wr*wpr-wi*wpi+wr
            wi            = wi*wpr+wtemp*wpi+wi
          end do
          h1r          = datap(i,j,1)
          datap(i,j,1) = c1*(h1r+datap(i,j,2))
          datap(i,j,2) = c1*(h1r-datap(i,j,2))
          call four1t(datap,-1,i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine four1t(dataff,isig,ii,jj)

      implicit none
      include 'par.for'
      integer isig, i, istep, z, m, mmax, ii, jj
      real*8 tempi, tempr, dataff(imax,jmax,kphys)
      real*8 theta, wi, wpi, wpr, wr, wtemp

      z = 1
      do i = 1, kphys, 2 
        if (z.gt.i) then
          tempr             = dataff(ii,jj,z)
          tempi             = dataff(ii,jj,z+1)
          dataff(ii,jj,z)   = dataff(ii,jj,i)
          dataff(ii,jj,z+1) = dataff(ii,jj,i+1)
          dataff(ii,jj,i)   = tempr
          dataff(ii,jj,i+1) = tempi
        end if
        m = kphys/2
    1   if ((m.ge.2).and.(z.gt.m)) then
          z = z - m
          m = m / 2
          goto 1
        endif
        z = z + m
      end do 
      mmax = 2
    2 if (kphys.gt.mmax) then
        istep = 2 * mmax
        theta = 6.28318530717959d0/dble(isig*mmax)
        wpr   = -2.d0*dsin(0.5d0*theta)**2
        wpi   = dsin(theta)
        wr    = 1.d0
        wi    = 0.d0
        do m = 1, mmax, 2
          do i = m, kphys, istep
            z                 = i + mmax
            tempr             = wr*dataff(ii,jj,z)-wi*dataff(ii,jj,z+1)
            tempi             = wr*dataff(ii,jj,z+1)+wi*dataff(ii,jj,z)
            dataff(ii,jj,z)   = dataff(ii,jj,i)-tempr
            dataff(ii,jj,z+1) = dataff(ii,jj,i+1)-tempi
            dataff(ii,jj,i)   = dataff(ii,jj,i)+tempr
            dataff(ii,jj,i+1) = dataff(ii,jj,i+1)+tempi
          end do 
          wtemp = wr
          wr    = wr*wpr-wi*wpi+wr
          wi    = wi*wpr+wtemp*wpi+wi
        end do 
        mmax = istep
        goto 2
      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c   end program that calculates Q for visualization     c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
