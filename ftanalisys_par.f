ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c             fourier analysis in time prg              c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fourier analisys 19012009

      implicit none
      include 'par.for'
      character*15 nome, nome2
      integer i, j, k, f, N, MM, jj
      parameter ( N = 16, MM = 4 )
      real*8 a(N), b(N), A1x(imax/2+1,jmax,N)
      real*8 yy, x(3), y(3), cof(3), ymax, kk, xad
      real*8 amaxx(imax/2+1,N), amax(imax/2+1,N), alfaux(imax/2+1,N)
      complex*16 ux_f(imax/2+1,jmax,kfour,N)

      call initval(ux_f,N)

      do k = 1, kfour

        do i = 1, imax / 2 + 1
          do j = 1, jmax
            do f = 1, N
              a(f) = dreal(ux_f(i,j,k,f))
              b(f) = dimag(ux_f(i,j,k,f))
            end do
            call fft1d ( a, b,.FALSE., N, MM)
            do f = 1, N / 2
c             if (k.eq.1.and.f.eq.1) then
c               kk = 4.d0
c              else if (k.eq.1.or.f.eq.1) then
c               kk = 2.d0
c              else
                kk = 1.d0
c             end if 
              A1x(i,j,f) = kk * dsqrt(a(f) * a(f) + b(f) * b(f))
            end do
          end do

          do f = 1, N/2
            yy = 1.d-30
            jj = 3
            do j = 1, jmax
              if (dabs(A1x(i,j,f)).gt.yy) then 
                jj = j
                yy = dabs(A1x(i,j,f))
              end if
            end do

            x(1) = dble(jj-2) * dy
            x(2) = dble(jj-1) * dy
            x(3) = dble(jj)   * dy

            y(1) = dabs(A1x(i,jj-1,f))
            y(2) = dabs(A1x(i,jj,f))
            y(3) = dabs(A1x(i,jj+1,f))

            call polcoe(x,y,cof)

            ymax       = - 0.5d0 * cof(2) / cof(3)
            amaxx(i,f) = cof(3) * ymax * ymax + cof(2) * ymax + cof(1)
            amax(i,f)  = dlog(amaxx(i,f))
          end do

        end do

        call derfax(alfaux, amax, N)

        do f = 1, N / 2

          write(nome,'(a,i0.2,a,i0.2,a)')'mode',f-1,'_', k-1,'.dat'
          write(nome2,'(a,i0.2,a,i0.2)')'mode',f-1,'_', k-1

          open (1, file = nome, status = 'unknown')
          write(1,*) 'VARIABLES="x","U_max","alfa_ux"'
          write(1,*) 'ZONE T="',nome2,'", I=',imax / 2 - 3
          do i = 5, imax / 2 + 1
            xad = x0 + dble(i-1) * dx * 2.d0
c           write(1,3) dsqrt(re*xad), amaxx(i,f), alfaux(i,f)
            write(1,3) xad, amaxx(i,f), alfaux(i,f)
          end do
          close (unit=1)

        end do

      end do

    3 format(1x,3d17.9)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initval(ux_f, N)

      ! subroutine to read the values
      implicit none
      include 'par.for'
      character*20 nm
      integer i, j, k, f, N, inter, shift, my_rank, t
      complex*16 ux(ptsx,jmax,kfour), ux_t(imax,jmax,kfour,N),
     &           ux_f(imax/2+1,jmax,kfour,N)

      inter = 2**( msh - 1 ) * ( stencil - 2 )
      do my_rank = 0, np - 1
        shift = my_rank * (ptsx - inter - 1)
        do f = 1, N
          write(nm,'(a,i0.2,a,i0.2)')'pert_',my_rank,'_',f
          open(unit=1,file=nm,form='unformatted')
          read(1) ux
          close(unit=1)
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx
                ux_t(i+shift,j,k,f) = ux(i,j,k)
              end do
            end do
          end do
        end do
      end do
      do f = 1, N
        do k = 1, kfour
          do j = 1, jmax
            do i = 1, imax/2+1
              ux_f(i,j,k,f) = ux_t(i*2-1,j,k,f)
            end do
          end do
        end do
      end do
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine polcoe(x, y, cof)

      implicit none
      integer i, j, k, n
      real*8 x(3), y(3), cof(3), b, ff, phi, s(3)

      n = 3
      do i = 1, n
        s(i)   = 0.d0
        cof(i) = 0.d0
      end do
      s(n) = -x(1)
      
      do i = 2, n
        do j = n+1-i, n-1
          s(j) = s(j) - x(i) * s(j+1)
        end do
        s(n) = s(n) - x(i)
      end do
      do j = 1, n
        phi = dble(n)
        do k = n-1, 1, -1
          phi = dble(k) * s(k+1) + x(j) * phi
        end do
        ff = y(j) / phi
        b = 1.d0
        do k = n, 1, -1
          cof(k) = cof(k) + b    * ff
          b      = s(k)   + x(j) * b
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derfax(ddx, fc, N)

c first derivatives calculation in x direction
      implicit none
      include 'par.for'
      integer i, f, N
      real*8 fc(imax/2+1,N), ddx(imax/2+1,N)

      do f = 1, N / 2
        do i = 4, imax / 2 + 1 - 3
          ddx(i,f) = ( -         fc(i-3,f) +  9.d0 * fc(i-2,f)
     &                 - 45.d0 * fc(i-1,f) + 45.d0 * fc(i+1,f)
     &                 -  9.d0 * fc(i+2,f) +         fc(i+3,f) )
     &               / ( 120.d0 * dx )
        end do
      end do
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fft1d(a, b, INV, N, eme)

      implicit none
      include 'par.for'
      logical INV
      integer j, nby2, l, k, me, mm, i, lpk, eme, N 
      real*8 a(n), b(n), tr, ti, sr, si, ninv
      real*8 ur(n/2), ui(n/2), wr, wi, sign, pi, t, pibyk
      data pi/3.141592653589793d0/

c***********************************************************************
c
c   ******************************************************************
c   *                                                                *
c   *                   Subroutine fft1d                             *
c   *                                                                *
c   ******************************************************************
c   Written by G. A. Fenton, Jan. 6, 1988.
c
c                 Single Precision Version
c
c    Computes the discrete Fourier transform of an observed series X(t)
c    which can be expressed as
c
c             X(t_i) = A(t_i) + iB(t_i)
c
c    where A(t_i) and B(t_i) are the real and imaginary parts of X(t) at the
c    time t_i. The vectors A and B are overwritten with the calculated fourier
c    coefficients A(k) and B(k) corresponding to the frequencies 2*pi*k/(N*s)
c    where s is the time increment. This routine can also perform the inverse
c    fourier transform by setting the flag INV.
c    If INV = .TRUE. then the original series X(t) is reconstructed from
c    the fourier coefficients A(k) and B(k).
c    Input and output variables are;
c
c      A   real vector of length N. On input A contains the real part of the
c          observed series. On output A contains the fourier coefficients A(k)
c          corresponding to the frequencies 2*pi*k/(N*s) where s is the time
c          increment.
c
c      B   real vector of length N. On input B contains the imaginary
c          components of the observed series. On ouput, B contains the fourier
c          coefficients B(k).
c
c      m   the value of m is such that N = 2**m
c
c      N   the number of data points entered. N must be a power of 2.
c
c    INV   logical flag which is false if the fourier coefficients are to be
c          calculated, true if the inverse transform is required.
c
c   ur,ui  temporary double precision vectors of length N/2
c
c   Note: this routine is modified from a routine presented by D.E. Newland
c         (originally written by J.W. Cooley et al.)
c
c
c............................................................................

      sign = 1.d0
      if ( .not. INV ) then
         sign = -1.d0
         ninv = 1.d0/dble(n)
         do j = 1,N
            a(j) = a(j) * ninv
            b(j) = b(j) * ninv
         end do
      endif

c.........................................................reorder vector
      nby2 = n/2
      j = 1

      do 30 l = 1, n-1
         if( l .lt. j ) then
             t    = a(j)
             a(j) = a(l)
             a(l) = t
             t    = b(j)
             b(j) = b(l)
             b(l) = t
         endif
         k = nby2
  20     if( k .ge. j ) go to 30
             j = j - k
             k = k / 2
             go to 20
  30  j = j + k
      

c.........................................Calculate fourier coefficients
      me = 1
      do mm = 1,eme

         k     = me
         pibyk = pi / dble(k)
         me    = 2 * me
         wr    = dcos(pibyk)
         wi    = sign * dsin(pibyk)

         ur(1) = 1.d0
         ui(1) = 0.d0
         do i = 2, k
            ur(i) = ur(i-1) * wr - ui(i-1) * wi
            ui(i) = ur(i-1) * wi + ui(i-1) * wr
         end do

         do j = 1,k
            sr = ur(j)
            si = ui(j)
            do l = j, n, me
               lpk    = l + k
               tr     = a(lpk) * sr - b(lpk) * si
               ti     = b(lpk) * sr + a(lpk) * si
               a(lpk) = a(l) - tr
               b(lpk) = b(l) - ti
               a(l)   = a(l) + tr
               b(l)   = b(l) + ti
            end do
         end do

      end do

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c      end of fourier analysis in time prg              c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
