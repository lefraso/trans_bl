ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c       subroutines transforms binary data in asc       c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program analysis 29052009

      implicit none
      include 'par.for'
      character c1
      character*2 c2
      character*11 nome
      character*19 nome2
      integer i, j, k, t, my_rank, inter, shift, var
      real*8  x, y, z,
     &        uxb(imax,jmax), uyb(imax,jmax), wzb(imax,jmax),
     &        thb(imax,jmax),
     &        uxp(imax,jmax,kphys),    wxp(imax,jmax,kphys),
     &        uyp(imax,jmax,kphys),    wyp(imax,jmax,kphys),
     &        uzp(imax,jmax,kphys),    wzp(imax,jmax,kphys),
     &        thp(imax,jmax,kphys),
     &        duxbdy(imax,jmax),       duxbdz(imax,jmax),
     &        dthbdy(imax,jmax),       dthbdz(imax,jmax),
     &        dthmdy(imax),            duxmdy(imax),
     &        duxdyp(imax,jmax,kphys), dthdyp(imax,jmax,kphys)

      complex*16 uxt(imax,jmax,kfour),    wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour),    wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour),    wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour),   
     &            ux(ptsx,jmax,kfour),     wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),     wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),     wz(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour),  duxdy(imax,jmax,kfour),
     &         dthdy(imax,jmax,kfour)

      call init_thermal

      ! Boundary layer data
      open(1,file='basef_ns.bin',form='unformatted')
      read(1) uxb, uyb, wzb, thb
      close(unit=1)

      ! Disturbances variables data
      inter = 2**( msh - 1 ) * ( stencil - 2 )
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(2,file=nome,form='unformatted')
        read(2) t
        read(2) ux,uy,uz,wx,wy,wz,th
        close (unit=2)
        shift = my_rank * (ptsx - inter - 1)
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, ptsx
              uxt(i+shift,j,k) = ux(i,j,k)
              uyt(i+shift,j,k) = uy(i,j,k)
              uzt(i+shift,j,k) = uz(i,j,k)
              wxt(i+shift,j,k) = wx(i,j,k)
              wyt(i+shift,j,k) = wy(i,j,k)
              wzt(i+shift,j,k) = wz(i,j,k)
              tht(i+shift,j,k) = th(i,j,k)
            end do
          end do
        end do
      end do

      call f_to_p(uxp,uxt)
      call f_to_p(uyp,uyt)
      call f_to_p(uzp,uzt)
      call f_to_p(wxp,wxt)
      call f_to_p(wyp,wyt)
      call f_to_p(wzp,wzt)
      call f_to_p(thp,tht)


      ! writes data to spacial space to be open by tecplot
      do i = 201, 738, 67
        var = 100*((i-1)*dx+x0)
        write(nome2,'(a,i0.4,a)')'crosscutxeq',var,'.dat'
        open (3, file = nome2,status = 'unknown')
        write(3,*) 'VARIABLES=,"z","y","uxb","uyb","wzb","thb",
     &    "uxp","uyp","uzp","thp","qq"'
        write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
        
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          z = -1.d0 * 6.28318530717959d0 / ( dble(kphys) * beta )
          write(3,5)z,y,
     &     uxb(i,j),uyb(i,j),wzb(i,j),thb(i,j),
     &     uxp(i,j,kphys),uyp(i,j,kphys),uzp(i,j,kphys),thp(i,j,kphys),
     &     0.d0
        
          do k = 1, kphys
            z = dble(k-1) * 6.28318530717959d0 / ( dble(kphys) * beta )
            write(3,5)z,y,
     &                uxb(i,j),uyb(i,j),wzb(i,j),thb(i,j),
     &                uxp(i,j,k),uyp(i,j,k),
     &                uzp(i,j,k),thp(i,j,k),0.d0
          end do
        end do
        close (unit=3)
      end do

    5 format(1x,2d14.6,9d17.9)

      call deryt(duxdy,uxt)
      call deryb(duxbdy,uxb)
      call deryt(dthdy,tht)
      call deryb(dthbdy,thb)

      call f_to_p(duxdyp,duxdy)
      call f_to_p(dthdyp,dthdy)

      do i = 1, imax
        duxmdy(i) = 0.d0
        dthmdy(i) = 0.d0
      end do
      do i = 1, imax
        do k = 1, kphys
          duxmdy(i) = duxmdy(i) + duxdyp(i,1,k)
          dthmdy(i) = dthmdy(i) + dthdyp(i,1,k)
        end do
        duxmdy(i) = duxmdy(i) / dble(kphys) + duxbdy(i,1)
        dthmdy(i) = dthmdy(i) / dble(kphys) + dthbdy(i,1)
      end do

      open (4, file = 'heat_coeffs.dat',status = 'unknown')
      write(4,*) 'VARIABLES=,"x","z","dudy","dthdy","dubdy","dthbdy",
     &"dumdy","dthmdy"'
      write(4,*) 'ZONE I=',imax,', K=',kphys+1,', F=POINT'

      z = -1.d0 * 6.28318530717959d0 / ( dble(kphys) * beta )
      do i = 1, imax
        x = x0 + dble(i-1) * dx
        write(4,6)x, z, duxdyp(i,1,kphys)+duxbdy(i,1), 
     &            dthdyp(i,1,kphys)+dthbdy(i,1),
     &            duxbdy(i,1), dthbdy(i,1),duxmdy(i), dthmdy(i)
      end do

      do k = 1, kphys
        z = dble(k-1) * 6.28318530717959d0 / ( dble(kphys) * beta )
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(4,6)x, z, duxdyp(i,1,k)+duxbdy(i,1), 
     &              dthdyp(i,1,k)+dthbdy(i,1),
     &              duxbdy(i,1), dthbdy(i,1),duxmdy(i), dthmdy(i)
        end do
      end do

    6 format(1x,2d14.6,6d17.9)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_thermal

      ! initialize the program main variables
      implicit none
      include 'par.for'
      include 'comm.coef'

      open(1,file='pre_processing/coefs.bin',form='unformatted')
      read(1) fp_fd_coef
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      close(unit=1)

      call create_ctes

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine create_ctes

      implicit none
      include 'par.for'
      include 'comm.multi'   
 
      integer lvl, j
 
      ! dy0 at each multigrid level
      do lvl = 1 , msh
       if(stf.ne.1.d0) then
        v_dy0(lvl) = dy * ( ( stf**(2**(lvl-1))-1.d0) / (stf-1.d0) )
       else
        v_dy0(lvl) = dy * 2.d0**(lvl-1)
       endif 
      enddo

      do lvl = 1 , msh
       v_stf(lvl) = stf ** ( 2**(lvl-1) )
      enddo

      ! dy at each space
       if (stf.ne.1.d0) then
        do j = 1 , jmax - 1
         v_qdy(j)  = 1.d0 / (v_dy0(1) * v_stf(1)**(j-1))
        enddo
       else
        do j = 1 , jmax - 1
         v_qdy(j)  = 1.d0 / (v_dy0(1))
        enddo
      endif

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs_k

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      common/der1y/ a1y,b1y,c1y

      call coef(a1y,b1y,c1y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coef(a,b,c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_fd_coef(1)
      c(1)      = fp_fd_coef(2)

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax - 2
        a(j)    = cp_fd_coef(1)
        b(j)    = cp_fd_coef(2)
        c(j)    = cp_fd_coef(3)
      end do

      a(jmax-1) = pp_fd_coef(3)
      b(jmax-1) = pp_fd_coef(2)
      c(jmax-1) = pp_fd_coef(1)

      a(jmax)   = lp_fd_coef(2)
      b(jmax)   = lp_fd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryt(ddy,fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1y/ a1y,b1y,c1y

      call rhsy(fc,ddy)
      call tridy(a1y,b1y,c1y,ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsy(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(imax,jmax,kfour), fc(imax,jmax,kfour)

      do k = 1, kfour
        do i = 1, imax
          rhs(i,1,k) = v_qdy(1)* ( fp_fd_coef(3) * fc(i,1,k) +
     &                             fp_fd_coef(4) * fc(i,2,k) +
     &                             fp_fd_coef(5) * fc(i,3,k) +
     &                             fp_fd_coef(6) * fc(i,4,k) +
     &                             fp_fd_coef(7) * fc(i,5,k) )
          
          rhs(i,2,k) = v_qdy(1)* ( sp_fd_coef(4) * fc(i,1,k) +
     &                             sp_fd_coef(5) * fc(i,2,k) +
     &                             sp_fd_coef(6) * fc(i,3,k) +
     &                             sp_fd_coef(7) * fc(i,4,k) +
     &                             sp_fd_coef(8) * fc(i,5,k) +
     &                             sp_fd_coef(9) * fc(i,6,k) )

         do j = 3, jmax - 2
           rhs(i,j,k) = v_qdy(j-2)* ( cp_fd_coef(4) * fc(i,j-2,k) +
     &                                cp_fd_coef(5) * fc(i,j-1,k) +
     &                                cp_fd_coef(6) * fc(i,j,k)   +
     &                                cp_fd_coef(7) * fc(i,j+1,k) +
     &                                cp_fd_coef(8) * fc(i,j+2,k) )
         end do
         
         rhs(i,jmax-1,k) = v_qdy(jmax-5)*(pp_fd_coef(4)*fc(i,jmax,k)
     &                                   +pp_fd_coef(5)*fc(i,jmax-1,k)
     &                                   +pp_fd_coef(6)*fc(i,jmax-2,k)
     &                                   +pp_fd_coef(7)*fc(i,jmax-3,k)
     &                                   +pp_fd_coef(8)*fc(i,jmax-4,k)
     &                                   +pp_fd_coef(9)*fc(i,jmax-5,k))
         
         rhs(i,jmax,k) = v_qdy(jmax-4)*(lp_fd_coef(3)*fc(i,jmax,k)
     &                                 +lp_fd_coef(4)*fc(i,jmax-1,k)
     &                                 +lp_fd_coef(5)*fc(i,jmax-2,k)
     &                                 +lp_fd_coef(6)*fc(i,jmax-3,k)
     &                                 +lp_fd_coef(7)*fc(i,jmax-4,k) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryb(ddy,fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      real*8 fc(imax,jmax), ddy(imax,jmax)
      common/der1y/ a1y,b1y,c1y

      call rhsyb(fc,ddy)
      call tridyb(a1y,b1y,c1y,ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsyb(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      real*8 rhs(imax,jmax), fc(imax,jmax)

      do i = 1, imax
        rhs(i,1) = v_qdy(1)* ( fp_fd_coef(3) * fc(i,1) +
     &                         fp_fd_coef(4) * fc(i,2) +
     &                         fp_fd_coef(5) * fc(i,3) +
     &                         fp_fd_coef(6) * fc(i,4) +
     &                         fp_fd_coef(7) * fc(i,5) )
        
        rhs(i,2) = v_qdy(1)* ( sp_fd_coef(4) * fc(i,1) +
     &                         sp_fd_coef(5) * fc(i,2) +
     &                         sp_fd_coef(6) * fc(i,3) +
     &                         sp_fd_coef(7) * fc(i,4) +
     &                         sp_fd_coef(8) * fc(i,5) +
     &                         sp_fd_coef(9) * fc(i,6) )

       do j = 3, jmax - 2
         rhs(i,j) = v_qdy(j-2)* ( cp_fd_coef(4) * fc(i,j-2) +
     &                            cp_fd_coef(5) * fc(i,j-1) +
     &                            cp_fd_coef(6) * fc(i,j)   +
     &                            cp_fd_coef(7) * fc(i,j+1) +
     &                            cp_fd_coef(8) * fc(i,j+2) )
       end do
       
       rhs(i,jmax-1) = v_qdy(jmax-5)*(pp_fd_coef(4)*fc(i,jmax)
     &                               +pp_fd_coef(5)*fc(i,jmax-1)
     &                               +pp_fd_coef(6)*fc(i,jmax-2)
     &                               +pp_fd_coef(7)*fc(i,jmax-3)
     &                               +pp_fd_coef(8)*fc(i,jmax-4)
     &                               +pp_fd_coef(9)*fc(i,jmax-5))
       
       rhs(i,jmax) = v_qdy(jmax-4)*(lp_fd_coef(3)*fc(i,jmax)
     &                             +lp_fd_coef(4)*fc(i,jmax-1)
     &                             +lp_fd_coef(5)*fc(i,jmax-2)
     &                             +lp_fd_coef(6)*fc(i,jmax-3)
     &                             +lp_fd_coef(7)*fc(i,jmax-4) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridy(a,b,c,rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rhs(imax,jmax,kfour), u(jmax)

      do k = 1, kfour
        do i = 1, imax
          bet  = b(1)
          u(1) = rhs(i,1,k) / bet
          do j = 2, jmax
            gam(j) = c(j-1) / bet
            bet    = b(j) - a(j) * gam(j)
            u(j)   = ( rhs(i,j,k) - a(j) * u(j-1) ) / bet
          end do
          do j = jmax - 1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
          end do
          do j = 1, jmax
            rhs(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridyb(a,b,c,rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      real*8 rhs(imax,jmax), u(jmax)

      do i = 1, imax
        bet  = b(1)
        u(1) = rhs(i,1) / bet
        do j = 2, jmax
          gam(j) = c(j-1) / bet
          bet    = b(j) - a(j) * gam(j)
          u(j)   = ( rhs(i,j) - a(j) * u(j-1) ) / bet
        end do
        do j = jmax - 1, 1, -1
          u(j) = u(j) - gam(j+1) * u(j+1)
        end do
        do j = 1, jmax
          rhs(i,j) = u(j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f_to_p(datap,dataf)

      ! transform from Fourier space to Physical space
      ! this subroutine works for a cosine for md=1
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
      c1    =  0.5d0
      wpr   = -2.0d0*dsin(0.5d0*theta)**2
      wpi   = dsin(theta)
      
      do j = 1, jmax
        do i = 1, imax
          do k = 0, kfour - 1
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
          call four1(datap,-1,i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine four1(dataff,isig,ii,jj)

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
        m=kphys/2
    1   if ((m.ge.2).and.(z.gt.m)) then
          z = z-m
          m = m/2
          goto 1
        endif
        z = z+m
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
c   end subroutines transforms binary data in asc       c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
