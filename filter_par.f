ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c                 filter calculations                   c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filter

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      integer i, j, k, i_ini
      real*8 luf(imax,5)
      complex*16 rhsx(ptsx,jmax), rhsy(ptsx,jmax), rhsz(ptsx,jmax)
      common/fil/ luf

      i_ini = 1
      if (my_rank.eq.0) i_ini = 2

      do k = 1, kfour
        call rhsf(rhsx,rhsy,rhsz,k)
        call lusolver(luf,rhsx)
        call lusolver(luf,rhsy)
        call lusolver(luf,rhsz)
        do j = 2, jmax
          do i = i_ini, ptsx
            wx(i,j,k) = rhsx(i,j)
            wy(i,j,k) = rhsy(i,j)
            wz(i,j,k) = rhsz(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lhsf(a)

      ! gives the LHS of the pentadiagonal matrix for the filter
      implicit none
      include 'par.for'
      integer i
      real*8 a(imax,5)

      a(1,1)      = 0.d0
      a(1,2)      = 0.d0
      a(1,3)      = 1.d0
      a(1,4)      = 0.d0
      a(1,5)      = 0.d0

      a(2,1)      = 0.d0
      a(2,2)      = 0.d0
      a(2,3)      = 1.d0
      a(2,4)      = 0.d0
      a(2,5)      = 0.d0

      a(3,1)      = 0.d0
      a(3,2)      = 0.d0
      a(3,3)      = 1.d0
      a(3,4)      = 0.d0
      a(3,5)      = 0.d0

      do i = 4, imax - 3
        a(i,1)    = 0.1702929d0
        a(i,2)    = 0.6522474d0
        a(i,3)    = 1.d0
        a(i,4)    = 0.6522474d0
        a(i,5)    = 0.1702929d0
      end do

      a(imax-2,1) = 0.d0
      a(imax-2,2) = 0.d0
      a(imax-2,3) = 1.d0
      a(imax-2,4) = 0.d0
      a(imax-2,5) = 0.d0

      a(imax-1,1) = 0.d0
      a(imax-1,2) = 0.d0
      a(imax-1,3) = 1.d0
      a(imax-1,4) = 0.d0
      a(imax-1,5) = 0.d0

      a(imax,1)   = 0.d0
      a(imax,2)   = 0.d0
      a(imax,3)   = 1.d0
      a(imax,4)   = 0.d0
      a(imax,5)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsf(rhsx,rhsy,rhsz,k)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      complex*16 rhsx(ptsx,jmax), rhsy(ptsx,jmax), rhsz(ptsx,jmax)

      do j = 1, jmax
        rhsx(1,j) = wx(1,j,k)
        rhsx(2,j) = wx(2,j,k)
        rhsx(3,j) = wx(3,j,k)
        rhsy(1,j) = wy(1,j,k)
        rhsy(2,j) = wy(2,j,k)
        rhsy(3,j) = wy(3,j,k)
        rhsz(1,j) = wz(1,j,k)
        rhsz(2,j) = wz(2,j,k)
        rhsz(3,j) = wz(3,j,k)

        do i = 4, ptsx - 3
          rhsx(i,j) = 0.989185600d0 *   wx(i,j,k) +
     &                0.660590000d0 * ( wx(i+1,j,k) + wx(i-1,j,k) ) +
     &                0.166677400d0 * ( wx(i+2,j,k) + wx(i-2,j,k) ) +
     &                0.000679925d0 * ( wx(i+3,j,k) + wx(i-3,j,k) )
          rhsy(i,j) = 0.989185600d0 *   wy(i,j,k) +
     &                0.660590000d0 * ( wy(i+1,j,k) + wy(i-1,j,k)) +
     &                0.166677400d0 * ( wy(i+2,j,k) + wy(i-2,j,k)) +
     &                0.000679925d0 * ( wy(i+3,j,k) + wy(i-3,j,k))
          rhsz(i,j) = 0.989185600d0 *   wz(i,j,k) +
     &                0.660590000d0 * ( wz(i+1,j,k) + wz(i-1,j,k)) +
     &                0.166677400d0 * ( wz(i+2,j,k) + wz(i-2,j,k)) +
     &                0.000679925d0 * ( wz(i+3,j,k) + wz(i-3,j,k))
        end do

        rhsx(ptsx-2,j) = wx(ptsx-2,j,k)
        rhsx(ptsx-1,j) = wx(ptsx-1,j,k)
        rhsx(ptsx,j)   = wx(ptsx,j,k)
        rhsy(ptsx-2,j) = wy(ptsx-2,j,k)
        rhsy(ptsx-1,j) = wy(ptsx-1,j,k)
        rhsy(ptsx,j)   = wy(ptsx,j,k)
        rhsz(ptsx-2,j) = wz(ptsx-2,j,k)
        rhsz(ptsx-1,j) = wz(ptsx-1,j,k)
        rhsz(ptsx,j)   = wz(ptsx,j,k)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filter_trid

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      integer i, j, k, i_ini
      real*8 afil(ptsx), bfil(ptsx), cfil(ptsx)
      complex*16 rhsx(ptsx,jmax,kfour), rhsy(ptsx,jmax,kfour), 
     &           rhsz(ptsx,jmax,kfour)
      common/filt/ afil, bfil, cfil

      i_ini = 1
      if (my_rank.eq.0) i_ini = 2

      call rhsf_trid(rhsx,rhsy,rhsz)
      call tridseqx(afil,bfil,cfil,rhsx)
      call boundary_exchange_derivs(rhsx)
      call tridseqx(afil,bfil,cfil,rhsy)
      call boundary_exchange_derivs(rhsy)
      call tridseqx(afil,bfil,cfil,rhsz)
      call boundary_exchange_derivs(rhsz)
      do i = i_ini, ptsx
        do j = 2, jmax
          do k = 1, kfour
            wx(i,j,k) = rhsx(i,j,k)
            wy(i,j,k) = rhsy(i,j,k)
            wz(i,j,k) = rhsz(i,j,k)
          end do
        end do
      end do

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lhs_tridf(a,b,c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include 'par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 0.d0
      b(2)      = 1.d0
      c(2)      = 0.d0

      a(3)      = 0.d0
      b(3)      = 1.d0
      c(3)      = 0.d0

      do i = 4, ptsx - 3
        a(i)    = alphaf
        b(i)    = 1.d0
        c(i)    = alphaf
      end do

      a(ptsx-2) = 0.d0
      b(ptsx-2) = 1.d0
      c(ptsx-2) = 0.d0

      a(ptsx-1) = 0.d0
      b(ptsx-1) = 1.d0
      c(ptsx-1) = 0.d0

      a(ptsx)   = 0.d0
      b(ptsx)   = 1.d0
      c(ptsx)   = 0.d0

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsf_trid(rhsx,rhsy,rhsz)

      ! calculate the RHS for the wx filter 
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      complex*16 rhsx(ptsx,jmax,kfour), rhsy(ptsx,jmax,kfour), 
     &           rhsz(ptsx,jmax,kfour)

      do k = 1, kfour
      do j = 1, jmax
        rhsx(1,j,k) = (   15.d0 * wx(1,j,k) +  4.d0 * wx(2,j,k)
     &                  -  6.d0 * wx(3,j,k) +  4.d0 * wx(4,j,k)
     &                  -         wx(5,j,k) ) / 16.d0
        rhsx(2,j,k) = (           wx(1,j,k) + 12.d0 * wx(2,j,k)
     &                  +  6.d0 * wx(3,j,k) -  4.d0 * wx(4,j,k)
     &                  +         wx(5,j,k) ) / 16.d0 
        rhsx(3,j,k) = ( -         wx(1,j,k) +  4.d0 * wx(2,j,k)
     &                  + 10.d0 * wx(3,j,k) +
     &                     4.d0 * wx(4,j,k) -         wx(5,j,k) )
     &                  / 16.d0 

        rhsy(1,j,k) = (   15.d0 * wy(1,j,k) +  4.d0 * wy(2,j,k)
     &                  -  6.d0 * wy(3,j,k) +  4.d0 * wy(4,j,k)
     &                  -         wy(5,j,k) ) / 16.d0 
        rhsy(2,j,k) = (           wy(1,j,k) + 12.d0 * wy(2,j,k)
     &                  +  6.d0 * wy(3,j,k) -  4.d0 * wy(4,j,k)
     &                  +         wy(5,j,k) ) / 16.d0 
        rhsy(3,j,k) = ( -         wy(1,j,k) +  4.d0 * wy(2,j,k)
     &                  + 10.d0 * wy(3,j,k) +  4.d0 * wy(4,j,k)
     &                  -         wy(5,j,k) ) / 16.d0

        rhsz(1,j,k) = (   15.d0 * wz(1,j,k) +  4.d0 * wz(2,j,k)
     &                  -  6.d0 * wz(3,j,k) +  4.d0 * wz(4,j,k)
     &                  -         wz(5,j,k) ) / 16.d0
        rhsz(2,j,k) = (           wz(1,j,k) + 12.d0 * wz(2,j,k)
     &                  +  6.d0 * wz(3,j,k) -  4.d0 * wz(4,j,k)
     &                  +         wz(5,j,k) ) / 16.d0
        rhsz(3,j,k) = ( -         wz(1,j,k) +  4.d0 * wz(2,j,k)
     &                  + 10.d0 * wz(3,j,k) +  4.d0 * wz(4,j,k)
     &                  -         wz(5,j,k) ) / 16.d0

        do i = 4, ptsx - 3
          rhsx(i,j,k) = af *   wx(i,j,k) +
     &                  bf * ( wx(i+1,j,k) + wx(i-1,j,k) ) +
     &                  cf * ( wx(i+2,j,k) + wx(i-2,j,k) ) +
     &                  df * ( wx(i+3,j,k) + wx(i-3,j,k) )
          rhsy(i,j,k) = af *   wy(i,j,k) +
     &                  bf * ( wy(i+1,j,k) + wy(i-1,j,k)) +
     &                  cf * ( wy(i+2,j,k) + wy(i-2,j,k)) +
     &                  df * ( wy(i+3,j,k) + wy(i-3,j,k))
          rhsz(i,j,k) = af *   wz(i,j,k) +
     &                  bf * ( wz(i+1,j,k) + wz(i-1,j,k)) +
     &                  cf * ( wz(i+2,j,k) + wz(i-2,j,k)) +
     &                  df * ( wz(i+3,j,k) + wz(i-3,j,k))
        end do

        rhsx(ptsx-2,j,k) = ( -         wx(ptsx,j,k)  
     &                       +  4.d0 * wx(ptsx-1,j,k)
     &                       + 10.d0 * wx(ptsx-2,j,k)
     &                       +  4.d0 * wx(ptsx-3,j,k)
     &                       -         wx(ptsx-4,j,k) ) / 16.d0 
        rhsx(ptsx-1,j,k) = (           wx(ptsx,j,k)
     &                       + 12.d0 * wx(ptsx-1,j,k)
     &                       +  6.d0 * wx(ptsx-2,j,k)
     &                       -  4.d0 * wx(ptsx-3,j,k)
     &                       +         wx(ptsx-4,j,k) ) / 16.d0 
        rhsx(ptsx,j,k)   = (   15.d0 * wx(ptsx,j,k)
     &                       +  4.d0 * wx(ptsx-1,j,k)
     &                       -  6.d0 * wx(ptsx-2,j,k)
     &                       +  4.d0 * wx(ptsx-3,j,k)
     &                       -         wx(ptsx-4,j,k) ) / 16.d0 

        rhsy(ptsx-2,j,k) = ( -         wy(ptsx,j,k)
     &                       +  4.d0 * wy(ptsx-1,j,k)
     &                       + 10.d0 * wy(ptsx-2,j,k)
     &                       +  4.d0 * wy(ptsx-3,j,k)
     &                       -         wy(ptsx-4,j,k) ) / 16.d0 
        rhsy(ptsx-1,j,k) = (           wy(ptsx,j,k)
     &                       + 12.d0 * wy(ptsx-1,j,k)
     &                       +  6.d0 * wy(ptsx-2,j,k)
     &                       -  4.d0 * wy(ptsx-3,j,k)
     &                       +         wy(ptsx-4,j,k) ) / 16.d0 
        rhsy(ptsx,j,k)   = (   15.d0 * wy(ptsx,j,k)
     &                       +  4.d0 * wy(ptsx-1,j,k)
     &                       -  6.d0 * wy(ptsx-2,j,k)
     &                       +  4.d0 * wy(ptsx-3,j,k)
     &                       -         wy(ptsx-4,j,k) ) / 16.d0 

        rhsz(ptsx-2,j,k) = ( -         wz(ptsx,j,k)
     &                       +  4.d0 * wz(ptsx-1,j,k)
     &                       + 10.d0 * wz(ptsx-2,j,k)
     &                       +  4.d0 * wz(ptsx-3,j,k)
     &                       -         wz(ptsx-4,j,k) ) / 16.d0 
        rhsz(ptsx-1,j,k) = (           wz(ptsx,j,k)
     &                       + 12.d0 * wz(ptsx-1,j,k)
     &                       +  6.d0 * wz(ptsx-2,j,k)
     &                       -  4.d0 * wz(ptsx-3,j,k)
     &                       +         wz(ptsx-4,j,k) ) / 16.d0 
        rhsz(ptsx,j,k)   = (   15.d0 * wz(ptsx,j,k)
     &                       +  4.d0 * wz(ptsx-1,j,k)
     &                       -  6.d0 * wz(ptsx-2,j,k)
     &                       +  4.d0 * wz(ptsx-3,j,k)
     &                       -         wz(ptsx-4,j,k) ) / 16.d0 

      end do
      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c            end of filter calculations                 c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
