ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c       program that calculates Q for visualization     c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program isoq 10092012

      implicit none
      include 'par.for'
      include 'comm.fourier'   
      integer i, j, k
      real*8    duxdxp(imax,jmax,kphys), duxdyp(imax,jmax,kphys),
     &          duydxp(imax,jmax,kphys), duydyp(imax,jmax,kphys),
     &          duzdxp(imax,jmax,kphys), duzdyp(imax,jmax,kphys),
     &          duzdzp(imax,jmax,kphys), duxdzp(imax,jmax,kphys),
     &          duydzp(imax,jmax,kphys),      Q(imax,jmax,kphys)
      complex*16 duxdx(imax,jmax,kfour),  duxdy(imax,jmax,kfour),
     &           duydx(imax,jmax,kfour),  duydy(imax,jmax,kfour),
     &           duzdx(imax,jmax,kfour),  duzdy(imax,jmax,kfour),
     &           duzdz(imax,jmax,kfour),  duxdz(imax,jmax,kfour),
     &           duydz(imax,jmax,kfour),  
     &             uxt(imax,jmax,kfour),    uyt(imax,jmax,kfour),
     &             uzt(imax,jmax,kfour)
      common/iso/ uxt,uyt,uzt

      call initval

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
     &                     duydyp(i,j,k) * duydyp(i,j,k) +
     &                     duzdzp(i,j,k) * duzdzp(i,j,k) +
     &                     2.d0 * ( duxdyp(i,j,k) * duydxp(i,j,k) +
     &                     duxdzp(i,j,k) * duzdxp(i,j,k) +
     &                     duydzp(i,j,k) * duzdyp(i,j,k) )  )
          end do
        end do
      end do

      call escreveq(Q)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initval

      implicit none
      include 'par.for'
      include 'comm.coef'
      character*15 nome
      integer i, j, t, my_rank, k, inter, shift
      real*8 uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax),
     &       thbt(imax,jmax)
      complex*16  ux(ptsx,jmax,kfour),  wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),  wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),  wz(ptsx,jmax,kfour),
     &           uxt(imax,jmax,kfour), uyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour), th(ptsx,jmax,kfour)
      common/iso/ uxt,uyt,uzt
      
      inter = 2**( msh - 1 ) * ( stencil - 2 )

      select case (my_form)

        case(0)
         do my_rank = 0, np - 1
           write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
           open(2,file=nome,form='unformatted')
           read(2) t
           read(2) ux,uy,uz,wx,wy,wz
           close (unit=2)
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

        case(1)
         open(1,file='baseflow2D/basens.bin',form='unformatted')
         read(1) uxbt, uybt, wzbt
         close(unit=1)
         do my_rank = 0, np - 1
           write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
           open(2,file=nome,form='unformatted')
           read(2) t
           read(2) ux,uy,uz,wx,wy,wz
           close (unit=2)
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
         do j = 1, jmax
           do i = 1, imax
             uxt(i,j,1) = uxt(i,j,1) + uxbt(i,j)
             uyt(i,j,1) = uyt(i,j,1) + uybt(i,j)
           end do
         end do

        case(2)
         open(1,file='baseflow2D/basens.bin',form='unformatted')
         read(1) uxbt, uybt, wzbt, thbt
         close(unit=1)
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
               end do
             end do
           end do
         end do
         do j = 1, jmax
           do i = 1, imax
             uxt(i,j,1) = uxt(i,j,1) + uxbt(i,j)
             uyt(i,j,1) = uyt(i,j,1) + uybt(i,j)
           end do
         end do
      end select

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

      call derivs_k

      call create_ctes

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine create_ctes

      implicit none
      include 'par.for'
      include 'comm.multi'   
      include 'comm.fourier'   
 
      integer lvl, j, k 
      real*8 aux
 
      ! Multigrid spatial calculations
 
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
       v_dx2(lvl) = (dx * dble(2**(lvl-1)))**2
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

      do lvl = 1 , msh
       if (stf.ne.1.d0) then
        do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
         aux           = (v_dy0(lvl) * v_stf(lvl)**(j-1))
         v_qdy2(j,lvl) = 1.d0 / (aux**2)
        enddo
       else
        do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
         aux           = (v_dy0(lvl))
         v_qdy2(j,lvl) = 1.d0 / (aux**2)
        enddo
       endif
      enddo

      v_ptsx(1) = ptsx
      v_ptsy(1) = jmax
      do lvl = 2 , msh
       v_ptsx(lvl) = (v_ptsx(lvl-1) + 1) / 2
       v_ptsy(lvl) = (v_ptsy(lvl-1) + 1) / 2
      enddo

      do k = 1 , kfour
       v_k2b2(k) = - dble(k-1) * dble(k-1) * beta * beta
       v_kb(k)   = - im * dble(k-1) * beta
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs_k

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include 'par.for'
      real*8  a1x(imax),  b1x(imax),  c1x(imax)
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/ a1x,b1x,c1x
      common/der1y/ a1y,b1y,c1y
      common/der1fv/ a1fv,b1fv,c1fv

      call coefx(a1x,b1x,c1x)
      call coefy(a1y,b1y,c1y)
      call coefy_fv(a1fv,b1fv,c1fv)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derxt(ddx,fc)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a1x(imax), b1x(imax), c1x(imax)
      complex*16 fc(imax,jmax,kfour), ddx(imax,jmax,kfour)
      common/der1x/ a1x,b1x,c1x

      call rhsx(fc,ddx)
      call tridseqx(a1x,b1x,c1x,ddx)
      
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
      subroutine deryfvt(ddy,fc)

      ! first derivative calculation in y direction for uy
      implicit none
      include 'par.for'
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      complex*16 fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1fv/ a1fv,b1fv,c1fv

      call rhsyfv(fc,ddy)
      call tridy(a1fv,b1fv,c1fv,ddy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsx(fc,rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(imax,jmax,kfour), fc(imax,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          rhs(1,j,k)=( - 74.d0 * fc(1,j,k) + 16.d0 * fc(2,j,k) +
     &                   72.d0 * fc(3,j,k) - 16.d0 * fc(4,j,k) +
     &                    2.d0 * fc(5,j,k) ) / ( 24.d0 * dx )
          
          rhs(2,j,k)=( - 406.d0 * fc(1,j,k) - 300.d0 * fc(2,j,k) +
     &                   760.d0 * fc(3,j,k) -  80.d0 * fc(4,j,k) +
     &                    30.d0 * fc(5,j,k) -   4.d0 * fc(6,j,k) ) /
     &               ( 120.d0 * dx )
          
          do i = 3, imax - 2
            rhs(i,j,k)=(           fc(i+2,j,k) - fc(i-2,j,k) +
     &                   28.d0 * ( fc(i+1,j,k) - fc(i-1,j,k) ) ) /
     &                 ( 12.d0 * dx )
          end do
          
          rhs(imax-1,j,k)=( - 406.d0 * fc(imax,j,k)
     &                      - 300.d0 * fc(imax-1,j,k)
     &                      + 760.d0 * fc(imax-2,j,k)
     &                      -  80.d0 * fc(imax-3,j,k)
     &                      +  30.d0 * fc(imax-4,j,k)
     &                      -   4.d0 * fc(imax-5,j,k) ) /
     &                    ( - 120.d0 * dx )
          
          rhs(imax,j,k)=( - 74.d0 * fc(imax,j,k)
     &                    + 16.d0 * fc(imax-1,j,k)
     &                    + 72.d0 * fc(imax-2,j,k)
     &                    - 16.d0 * fc(imax-3,j,k)
     &                    +  2.d0 * fc(imax-4,j,k) ) /
     &                  ( - 24.d0 * dx )
        end do
      end do

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
      subroutine rhsyfv(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(imax,jmax,kfour), fc(imax,jmax,kfour)

      do k = 1, kfour
        do i = 1, imax
          rhs(i,1,k) = dcmplx(0.d0,0.d0)
          
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
      subroutine tridseqx(a,b,c,rhs)

      ! solves tridiagonal matrix for the derivatives in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a(imax), b(imax), c(imax), gam(imax), bet
      complex*16 rhs(imax,jmax,kfour), u(imax)

      do k = 1, kfour
        do j = 1, jmax
          bet  = b(1)
          u(1) = rhs(1,j,k) / bet
          do i = 2, imax
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i)   = ( rhs(i,j,k) - a(i) * u(i-1) ) / bet
          end do
          do i = imax - 1, 1, -1
            u(i) = u(i) - gam(i+1) * u(i+1)
          end do
          do i = 1, imax
            rhs(i,j,k) = u(i)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coefy(a,b,c)

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
      subroutine coefy_fv(a,b,c)

      ! mount the LHS of the matrix for the first derivative of uy
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

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
      subroutine coefx(a,b,c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax), b(imax), c(imax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, imax - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(imax-1) = 2.d0
      b(imax-1) = 6.d0
      c(imax-1) = 1.d0

      a(imax)   = 4.d0
      b(imax)   = 1.d0
      c(imax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine escreveq(Q)

      ! write the results in fourier modes
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 x, y, z
      real*8 Q(imax,jmax,kphys)

      ! writes data to spacial space to be open by tecplot
      open (3, file = 'isoq.dat',status = 'unknown')
      write(3,*) 'VARIABLES="x","y","z","Q"'
      write(3,*) 'ZONE I=',imax/2+1,',J=',jmax,',K=',kphys+1,',F=POINT'

      z = -1.d0 * dz
      do j = 1, jmax
        if(stf.eq.1.d0) then
         y = dble(j-1) * dy
        else
         y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
        endif        
        do i = 1, imax, 2
          x = x0 + dble(i-1) * dx
          if (abs(Q(i,j,kphys)).lt.1e-15) Q(i,j,kphys) = 0.d0
          write(3,5)x,y,z,Q(i,j,kphys)
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

      theta = - 6.28318530717959d0 / dble(kphys)
      c1    = 0.5d0
      wpr   = - 2.0d0 * dsin(0.5d0 * theta)**2
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
          datap(i,j,1) = datap(i,j,1) * 2.d0
          wr           = 1.0d0 + wpr
          wi           = wpi
          n2p3         = kphys + 3
          do p = 2, kphys/4
            p1            = 2 * p - 1
            p2            = p1 + 1
            p3            = n2p3 - p2
            p4            = p3 + 1
            wrs           = wr
            wis           = wi
            h1r           =   c1 * ( datap(i,j,p1) + datap(i,j,p3) )
            h1i           =   c1 * ( datap(i,j,p2) - datap(i,j,p4) )
            h2r           = - c1 * ( datap(i,j,p2) + datap(i,j,p4) )
            h2i           =   c1 * ( datap(i,j,p1) - datap(i,j,p3) )
            datap(i,j,p1) =   h1r + wrs * h2r - wis * h2i
            datap(i,j,p2) =   h1i + wrs * h2i + wis * h2r
            datap(i,j,p3) =   h1r - wrs * h2r + wis * h2i
            datap(i,j,p4) = - h1i + wrs * h2i + wis * h2r
            wtemp         = wr
            wr            = wr * wpr - wi    * wpi + wr
            wi            = wi * wpr + wtemp * wpi + wi
          end do
          h1r          = datap(i,j,1)
          datap(i,j,1) = c1 * ( h1r + datap(i,j,2) )
          datap(i,j,2) = c1 * ( h1r - datap(i,j,2) )
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
        theta = 6.28318530717959d0 / dble(isig * mmax)
        wpr   = - 2.d0 * dsin(0.5d0 * theta)**2
        wpi   = dsin(theta)
        wr    = 1.d0
        wi    = 0.d0
        do m = 1, mmax, 2
          do i = m, kphys, istep
            z                 = i + mmax
            tempr             = wr * dataff(ii,jj,z)
     &                        - wi * dataff(ii,jj,z+1)
            tempi             = wr * dataff(ii,jj,z+1)
     &                        + wi * dataff(ii,jj,z)
            dataff(ii,jj,z)   = dataff(ii,jj,i)   - tempr
            dataff(ii,jj,z+1) = dataff(ii,jj,i+1) - tempi
            dataff(ii,jj,i)   = dataff(ii,jj,i)   + tempr
            dataff(ii,jj,i+1) = dataff(ii,jj,i+1) + tempi
          end do 
          wtemp = wr
          wr    = wr * wpr - wi    * wpi + wr
          wi    = wi * wpr + wtemp * wpi + wi
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
