ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c                 derivative calculations               c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs_k

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include 'par.for'
      real*8  a1x(ptsx),  b1x(ptsx),  c1x(ptsx)
      real*8  a2x(ptsx),  b2x(ptsx),  c2x(ptsx)
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      real*8  a2y(jmax),  b2y(jmax),  c2y(jmax)
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/ a1x,b1x,c1x
      common/der2x/ a2x,b2x,c2x
      common/der1y/ a1y,b1y,c1y
      common/der2y/ a2y,b2y,c2y
      common/der1fv/ a1fv,b1fv,c1fv

      call coefx(a1x,b1x,c1x)
      call coeffx(a2x,b2x,c2x)
      call coefy(a1y,b1y,c1y)
      call coefy_fv(a1fv,b1fv,c1fv)
      call coeffy(a2y,b2y,c2y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derparx(ddx,fc)

      ! first derivatives calculation in x direction
      implicit none
      include 'par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx)
      complex*16 fc(ptsx,jmax,kfour), ddx(ptsx,jmax,kfour)
      common/der1x/ a1x,b1x,c1x

      call rhsx(fc,ddx)
      call tridparx(a1x,b1x,c1x,ddx)
c     call tridseqx(a1x,b1x,c1x,ddx)
c     call boundary_exchange_derivs(ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derparxx(d2dx2,fc)

      ! second derivative calculations in x direction
      implicit none
      include 'par.for'
      real*8 a2x(ptsx), b2x(ptsx), c2x(ptsx)
      complex*16 fc(ptsx,jmax,kfour), d2dx2(ptsx,jmax,kfour)
      common/der2x/ a2x,b2x,c2x

      call rhsxx(fc,d2dx2)
      call tridparx(a2x,b2x,c2x,d2dx2)
c     call tridseqx(a2x,b2x,c2x,d2dx2)
c     call boundary_exchange_derivs(d2dx2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dery(ddy,fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 fc(ptsx,jmax,kfour), ddy(ptsx,jmax,kfour)
      common/der1y/ a1y,b1y,c1y

      call rhsy(fc,ddy)
      call tridy(a1y,b1y,c1y,ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryfv(ddy,fc)

      ! first derivative calculation in y direction for uy
      implicit none
      include 'par.for'
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      complex*16 fc(ptsx,jmax,kfour), ddy(ptsx,jmax,kfour)
      common/der1fv/ a1fv,b1fv,c1fv

      call rhsyfv(fc,ddy)
      call tridy(a1fv,b1fv,c1fv,ddy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine deryy(d2dy2,fc)

      ! second derivative calculation in y direction
      implicit none
      include 'par.for'
      real*8 a2y(jmax), b2y(jmax), c2y(jmax)
      complex*16 fc(ptsx,jmax,kfour), d2dy2(ptsx,jmax,kfour)
      common/der2y/ a2y,b2y,c2y

      call rhsyy(fc,d2dy2)
      call tridy(a2y,b2y,c2y,d2dy2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsx(fc,rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          rhs(1,j,k)=( - 74.d0 * fc(1,j,k) + 16.d0 * fc(2,j,k) +
     &                   72.d0 * fc(3,j,k) - 16.d0 * fc(4,j,k) +
     &                    2.d0 * fc(5,j,k) ) / ( 24.d0 * dx )
          
          rhs(2,j,k)=( - 406.d0 * fc(1,j,k) - 300.d0 * fc(2,j,k) +
     &                   760.d0 * fc(3,j,k) -  80.d0 * fc(4,j,k) +
     &                    30.d0 * fc(5,j,k) -   4.d0 * fc(6,j,k) ) /
     &               ( 120.d0 * dx )
          
          do i = 3, ptsx - 2
            rhs(i,j,k)=(           fc(i+2,j,k) - fc(i-2,j,k) +
     &                   28.d0 * ( fc(i+1,j,k) - fc(i-1,j,k) ) ) /
     &                 ( 12.d0 * dx )
          end do
          
          rhs(ptsx-1,j,k)=( - 406.d0 * fc(ptsx,j,k)
     &                      - 300.d0 * fc(ptsx-1,j,k)
     &                      + 760.d0 * fc(ptsx-2,j,k)
     &                      -  80.d0 * fc(ptsx-3,j,k)
     &                      +  30.d0 * fc(ptsx-4,j,k)
     &                      -   4.d0 * fc(ptsx-5,j,k) ) /
     &                    ( - 120.d0 * dx )
          
          rhs(ptsx,j,k)=( - 74.d0 * fc(ptsx,j,k)
     &                    + 16.d0 * fc(ptsx-1,j,k)
     &                    + 72.d0 * fc(ptsx-2,j,k)
     &                    - 16.d0 * fc(ptsx-3,j,k)
     &                    +  2.d0 * fc(ptsx-4,j,k) ) /
     &                  ( - 24.d0 * dx )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsxx(fc,rhs)
      
      ! RHS for the second derivative calculation in x direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do j = 1, jmax
          rhs(1,j,k)=(  9775.d0 * fc(1,j,k) - 20285.d0 * fc(2,j,k)
     &               + 11170.d0 * fc(3,j,k) -   550.d0 * fc(4,j,k)
     &               -   145.d0 * fc(5,j,k) +    35.d0 * fc(6,j,k) ) /
     &               ( 60.d0 * dxx )
          
          rhs(2,j,k)=( 4834.d0 * fc(1,j,k) - 8424.d0 * fc(2,j,k)
     &               + 1890.d0 * fc(3,j,k) + 2320.d0 * fc(4,j,k)
     &               -  810.d0 * fc(5,j,k) +  216.d0 * fc(6,j,k)
     &               -   26.d0 * fc(7,j,k) ) / ( 360.d0 * dxx )
          
          do i = 3, ptsx - 2
            rhs(i,j,k)=(   3.d0 * ( fc(i+2,j,k) + fc(i-2,j,k) )
     &                 +  48.d0 * ( fc(i+1,j,k) + fc(i-1,j,k) )
     &                 - 102.d0 *   fc(i,j,k) ) / ( 4.d0 * dxx )
          end do
          
         rhs(ptsx-1,j,k)=(   4834.d0 * fc(ptsx,j,k)
     &                     - 8424.d0 * fc(ptsx-1,j,k)
     &                     + 1890.d0 * fc(ptsx-2,j,k)
     &                     + 2320.d0 * fc(ptsx-3,j,k)
     &                     -  810.d0 * fc(ptsx-4,j,k)
     &                     +  216.d0 * fc(ptsx-5,j,k)
     &                     -   26.d0 * fc(ptsx-6,j,k) ) /
     &                   ( 360.d0 * dxx )
          
          rhs(ptsx,j,k)=(    9775.d0 * fc(ptsx,j,k)
     &                    - 20285.d0 * fc(ptsx-1,j,k)
     &                    + 11170.d0 * fc(ptsx-2,j,k)
     &                    -   550.d0 * fc(ptsx-3,j,k)
     &                    -   145.d0 * fc(ptsx-4,j,k)
     &                    +    35.d0 * fc(ptsx-5,j,k) ) /
     &                  ( 60.d0 * dxx )
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
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do i = 1, ptsx
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
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour)

      do k = 1, kfour
        do i = 1, ptsx
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
      subroutine rhsyy(fc,rhs)

      ! RHS for the second derivative calculation in y direction
      implicit none
      include 'par.for'
      include 'comm.coef'
      include 'comm.multi'
      integer i, j, k
      complex*16 rhs(ptsx,jmax,kfour), fc(ptsx,jmax,kfour),a

      do k = 1, kfour
        do i = 1, ptsx

          rhs(i,1,k) = ( fp_sd_coef(3) * fc(i,1,k) +
     &                   fp_sd_coef(4) * fc(i,2,k) +
     &                   fp_sd_coef(5) * fc(i,3,k) +
     &                   fp_sd_coef(6) * fc(i,4,k) +
     &                   fp_sd_coef(7) * fc(i,5,k) +
     &                   fp_sd_coef(8) * fc(i,6,k) ) *
     &                   v_qdy2(1,1)
          
          rhs(i,2,k) = ( sp_sd_coef(4)  * fc(i,1,k) +
     &                   sp_sd_coef(5)  * fc(i,2,k) +
     &                   sp_sd_coef(6)  * fc(i,3,k) +
     &                   sp_sd_coef(7)  * fc(i,4,k) +
     &                   sp_sd_coef(8)  * fc(i,5,k) +
     &                   sp_sd_coef(9)  * fc(i,6,k) +
     &                   sp_sd_coef(10) * fc(i,7,k) ) *
     &                   v_qdy2(1,1)
          
         do j = 3, jmax - 2
           rhs(i,j,k) = ( cp_sd_coef(4) * fc(i,j-2,k) +
     &                    cp_sd_coef(5) * fc(i,j-1,k) +
     &                    cp_sd_coef(6) * fc(i,j,k)   +
     &                    cp_sd_coef(7) * fc(i,j+1,k) +
     &                    cp_sd_coef(8) * fc(i,j+2,k) ) * 
     &                    v_qdy2(j-2,1)
         end do
          
        rhs(i,jmax-1,k) = ( pp_sd_coef(4)  * fc(i,jmax,k)   +
     &                      pp_sd_coef(5)  * fc(i,jmax-1,k) +
     &                      pp_sd_coef(6)  * fc(i,jmax-2,k) +
     &                      pp_sd_coef(7)  * fc(i,jmax-3,k) +
     &                      pp_sd_coef(8)  * fc(i,jmax-4,k) +
     &                      pp_sd_coef(9)  * fc(i,jmax-5,k) +
     &                      pp_sd_coef(10) * fc(i,jmax-6,k) ) *
     &                      v_qdy2(jmax-6,1)
         
        rhs(i,jmax,k) = ( lp_sd_coef(3) * fc(i,jmax,k)   +
     &                    lp_sd_coef(4) * fc(i,jmax-1,k) +
     &                    lp_sd_coef(5) * fc(i,jmax-2,k) +
     &                    lp_sd_coef(6) * fc(i,jmax-3,k) +
     &                    lp_sd_coef(7) * fc(i,jmax-4,k) +
     &                    lp_sd_coef(8) * fc(i,jmax-5,k) ) * 
     &                    v_qdy2(jmax-5,1)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridparx(a,b,c,rhs)
      
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, k, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx)
      complex*16 aux(2*inter), rhs(ptsx,jmax,kfour)
      complex*16 u(ptsx,jmax), gam(ptsx), bet
      
      do k = 1, kfour
        do j = 1, jmax
          if (my_rank.eq.0) then
            bet    = b(1)
            u(1,j) = rhs(1,j,k) / bet
            gam(2) = c(1) / bet
            bet    = b(2) - a(2) * gam(2)
            u(2,j) = ( rhs(2,j,k) - a(2) * u(1,j) ) / bet
           else
            call MPI_Recv(aux, 2*inter, MPI_complex16, my_rank-1,
     &                    57, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(i,j) = aux(i)
              gam(i) = aux(i+inter)
            end do
          end if
          i_ini = inter
          if (my_rank.eq.0) i_ini = 3
          i_fim = ptsx - 2
          if (my_rank.eq.numproc) i_fim = ptsx
          
          bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
          do i = i_ini, i_fim
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i,j) = ( rhs(i,j,k) - a(i) * u(i-1,j) ) / bet
          end do
          
          if (my_rank.lt.numproc) then
            do i = 1, inter
              aux(i)       = u(ptsx-inter+i-1,j)
              aux(i+inter) = gam(ptsx-inter+i-1)
            end do
            call MPI_Send(aux, 2*inter, MPI_complex16, my_rank+1,
     &                    57, MPI_COMM_WORLD, ierr)
          end if
        end do
        
        do j = 1, jmax
          if (my_rank.lt.numproc) then
            call MPI_Recv(aux, 2*inter, MPI_complex16, my_rank+1, 
     &                    67, MPI_COMM_WORLD, status, ierr)
            do i = 1, inter
              u(ptsx-inter+i,j) = aux(i)
              gam(ptsx-inter+i) = aux(i+inter)
            end do
          end if
          
          i_fim = ptsx - inter
          if (my_rank.eq.numproc) i_fim = ptsx - 1
          do i = i_fim, 1, -1
            u(i,j) = u(i,j) - gam(i+1) * u(i+1,j)
          end do
          
          if (my_rank.gt.0) then
            do i = 1, inter
              aux(i)       = u(i+1,j)
              aux(i+inter) = gam(i+1)
            end do
            call MPI_Send(aux, 2*inter, MPI_complex16, my_rank-1,
     &                    67, MPI_COMM_WORLD, ierr)
          end if
          do i = 1, ptsx
            rhs(i,j,k) = u(i,j)
          end do
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
      complex*16 rhs(ptsx,jmax,kfour), u(jmax)

      do k = 1, kfour
        do i = 1, ptsx
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
      real*8 a(ptsx), b(ptsx), c(ptsx), gam(ptsx), bet
      complex*16 rhs(ptsx,jmax,kfour), u(ptsx)

      do k = 1, kfour
        do j = 1, jmax
          bet  = b(1)
          u(1) = rhs(1,j,k) / bet
          do i = 2, ptsx
            gam(i) = c(i-1) / bet
            bet    = b(i) - a(i) * gam(i)
            u(i)   = ( rhs(i,j,k) - a(i) * u(i-1) ) / bet
          end do
          do i = ptsx - 1, 1, -1
            u(i) = u(i) - gam(i+1) * u(i+1)
          end do
          do i = 1, ptsx
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
      subroutine coeffy(a,b,c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include 'par.for'
      include 'comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_sd_coef(1)
      c(1)      = fp_sd_coef(2)

      a(2)      = sp_sd_coef(1)
      b(2)      = sp_sd_coef(2)
      c(2)      = sp_sd_coef(3)

      do j = 3, jmax - 2
        a(j)    = cp_sd_coef(1)
        b(j)    = cp_sd_coef(2)
        c(j)    = cp_sd_coef(3)
      end do

      a(jmax-1) = pp_sd_coef(3)
      b(jmax-1) = pp_sd_coef(2)
      c(jmax-1) = pp_sd_coef(1)

      a(jmax)   = lp_sd_coef(2)
      b(jmax)   = lp_sd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coefx(a,b,c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, ptsx - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(ptsx-1) = 2.d0
      b(ptsx-1) = 6.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 4.d0
      b(ptsx)   = 1.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coeffx(a,b,c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 13.d0
      c(1)      = 137.d0

      a(2)      = 1.d0
      b(2)      = 12.d0
      c(2)      = 3.d0

      do i = 3, ptsx - 2
        a(i)    = 2.d0
        b(i)    = 11.d0
        c(i)    = 2.d0
      end do

      a(ptsx-1) = 3.d0
      b(ptsx-1) = 12.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 137.d0
      b(ptsx)   = 13.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundary_exchange_derivs(var)

      ! exchange values of the boundaries
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, j, k, tam
      complex*16 aux(meshdx*jmax*kfour), var(ptsx,jmax,kfour)

      ! variable used to calculate the number of columns needed 
      ! to go forward or backward (interm)

      tam = jmax*meshdx 

      if (my_rank.lt.numproc) then
        ! Sending the near right boundary columns to node + 1
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax+(k-1)*tam) = var(i+ptsx-24-1,j,k)
          end do
        end do
       end do 
        call MPI_Send(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank+1,
     &                10, MPI_COMM_WORLD, ierr)
        ! Receiving the new right boundary columns from node + 1
        call MPI_Recv(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank+1,
     &                20, MPI_COMM_WORLD, status, ierr)
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            var(i+ptsx-meshdx,j,k) = aux(j+(i-1)*jmax+(k-1)*tam)
          end do
        end do
       end do 
       
       end if

      if (my_rank.gt.0) then
        ! Receiving the new left boundary columns from node - 1
        call MPI_Recv(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank-1,
     &                10, MPI_COMM_WORLD, status, ierr)
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            var(i,j,k) = aux(j+(i-1)*jmax+(k-1)*tam)
          end do
        end do
       end do 
        ! Sending the new right boundary columns to node - 1
       do k = 1, kfour 
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax+(k-1)*tam) = var(i+24-meshdx+1,j,k)
          end do
        end do
       end do 
        call MPI_Send(aux, meshdx*jmax*kfour, MPI_COMPLEX16, my_rank-1,
     &                20, MPI_COMM_WORLD, ierr)
      end if

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c            end of derivative calculations             c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
