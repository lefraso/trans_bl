ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c        program that calculate the base flow           c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program Falkner_Skan 19012009

      ! this program generates Falkner_Skan boundary layer
      implicit none
      include 'par.for'
      integer i, j, k, cont, jjmax, jj
      real*8 eig(3,20*jmax), y(3), m, ue, betafs, dya
      real*8 f, fp, fpp, xad, h, lx, eta, yad, uad, vad, vort
      real*8 ux(imax,jmax), uy(imax,jmax), wz(imax,jmax)

      ! definition of betafs, Falkner-Scan parameter if = 0 then bondary layer=blasius     
      betafs = 0.d0
      jjmax  = (jmax - 1) * 20 + 1
      m      = betafs / (2.d0 - betafs)
      dya    = dy * dsqrt(Re * x0)

      do i = 1, imax
        xad = x0 + dble(i - 1) * dx
        ue  = 1.d0 * (xad)**m
        lx  = dsqrt(ue/xad)
        h   = dya * lx / 20.d0

        ! calculates the f2 at wall
        call rtsaf(fpp, h, betafs)

        ! initial values for derivatives
        y(1) = 0.d0
        y(2) = 0.d0
        y(3) = fpp
        cont = 1

        do j = 1, 3
          eig(j,1) = y(j)
        end do

        do k = 1, jjmax
          call rk4(y, h, betafs)
          do j = 1, 3
            eig(j,k+1) = y(j)
          end do
          if (cont.eq.k) then
            f    = eig(1,k)
            fp   = eig(2,k)
            fpp  = eig(3,k)
            eta  = dble(cont-1) * h
            yad  = eta / lx
            uad  = fp * ue
            vad  = lx * ( (1.d0 - betafs) * eta * fp - f ) /
     &             (2.d0 - betafs)
            vort = lx * ue * fpp 
            jj = (k - 1) / 20 + 1
            ux(i,jj) = uad
            uy(i,jj) = vad  / dsqrt(Re)
            wz(i,jj) = vort * dsqrt(Re)
            cont = cont + 20
          endif
        end do

      end do

      open(1,file='blasius.bin',form='unformatted')
      write(1) ux,uy,wz
      close (unit=1)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine rtsaf(rtsafe, h, betafs)

       implicit none
       integer j
       real*8 rtsafe, x1, x2, xacc, h, betafs
       real*8 df, dxx, dxold, f, fh, fl, temp, xh, xl

       x1   = 0.d0
       x2   = 1.5d0
       xacc = 2.d-16

       call funcd(x1, fl, df, h, betafs)
       call funcd(x2, fh, df, h, betafs)

       if(fl.lt.0.d0) then
         xl = x1
         xh = x2
        else
         xh = x1
         xl = x2
       endif
       rtsafe = 0.5d0 * ( x1 + x2 )
       dxold  = dabs(x2 - x1)
       dxx    = dxold
       call funcd(rtsafe, f, df, h, betafs)
       do j = 1, 150
         if ( ((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f).ge.0.d0
     &     .or. abs(2.d0*f).gt.abs(dxold*df) ) then
           dxold  = dxx
           dxx    = 0.5d0 * (xh-xl)
           rtsafe = xl+dxx
          else
           dxold  = dxx
           dxx    = f / df
           temp   = rtsafe
           rtsafe = rtsafe-dxx
         endif

         call funcd(rtsafe,f,df,h,betafs)
         if(abs(f).lt.xacc) return
         if(f.lt.0.d0)then
           xl = rtsafe
          else
           xh = rtsafe
         endif

       end do

       write(*,*)'rtsafe exceeding maximum iteractions'

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine funcd(x,fc,dfc,h,betafs)

      implicit none
      real*8 x, dx, fc2, betafs
      real*8 y(3), fc, dfc, ys(3), h

      y(1) = 0.d0
      y(2) = 0.d0
      y(3) = x

      dx   = 1.d-7

      call rkdumb(y,h,betafs)
      fc = y(2)-1.d0

      y(1) = 0.d0
      y(2) = 0.d0
      y(3) = x + dx

      call rkdumb(y,h,betafs)
      fc2 = y(2)-1.d0

      dfc = (fc2-fc)/dx

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rkdumb(y,h,betafs)

      implicit none
      include 'par.for'
      integer i, npt
      real*8 y(3), h, h2, betafs

      npt = 20*jmax
      h2  = h

      if (betafs.ge.0.1d0) h2 = 6.02d0/dble(npt)
      if (betafs.ge.0.2d0) h2 = 5.2d0/dble(npt)
      if (betafs.ge.0.4d0) h2 = 4.6d0/dble(npt)
      if (betafs.ge.0.5d0) h2 = 4.5d0/dble(npt)

      do i = 1, npt
        call rk4(y,h2,betafs)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rk4(y,h,betafs)

      implicit none
      integer i
      real*8 h, hh, h6, betafs
      real*8 y(3), dydx(3), yt(3), dyt(3), dym(3)

      hh = h * 0.5d0
      h6 = h / 6.d0

      call derivs(y,dydx,betafs)
      do i = 1, 3
         yt(i) = y(i) + hh * dydx(i)
      end do

      call derivs(yt,dyt,betafs) 
      do i = 1, 3
         yt(i) = y(i) + hh * dyt(i)
      end do

      call derivs(yt,dym,betafs) 
      do i = 1, 3
         yt(i)  = y(i) + h * dym(i)
         dym(i) = dyt(i) + dym(i)
      end do

      call derivs(yt,dyt,betafs) 
      do i = 1, 3
        y(i) = y(i) + h6 * ( dydx(i) + dyt(i) + 2.d0*dym(i) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derivs (y,dydx,betafs)

      implicit none
      include 'par.for'
      real*8 y(3), dydx(3), betafs

      dydx(1) =  y(2)
      dydx(2) =  y(3)
      dydx(3) = -y(1) * y(3) / (2.d0-betafs) -
     &           betafs * ( 1.d0 - y(2) * y(2) )/( 2.d0 - betafs )

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c    end program that calculate the base flow           c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
