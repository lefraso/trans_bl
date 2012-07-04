ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c         nonlinear terms calculation                   c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lterms(a,b,c)

      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8 uxb(ptsx,jmax),uyb(ptsx,jmax),wzb(ptsx,jmax)
      complex*16 a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &           c(ptsx,jmax,kfour)
      common/blas/ uxb,uyb,wzb

      do j = 1, jmax
        do i = 1, ptsx
          a(i,j,1) =   uyb(i,j) * wx(i,j,1) - wy(i,j,1) * uxb(i,j)
          b(i,j,1) =   uxb(i,j) * (wz(i,j,1)-wzb(i,j)) + 
     &                 wzb(i,j) * (ux(i,j,1)-uxb(i,j))
          c(i,j,1) = - uyb(i,j) * (wz(i,j,1)-wzb(i,j)) - 
     &                 wzb(i,j) * (uy(i,j,1)-uyb(i,j))
        end do
      end do

      do k = 2, kfour
        do j = 1, jmax
          do i = 1, ptsx
            a(i,j,k) =   uyb(i,j) * wx(i,j,k) - wy(i,j,k) * uxb(i,j)
            b(i,j,k) =   uxb(i,j) * wz(i,j,k) + ux(i,j,k) * wzb(i,j)
            c(i,j,k) = - uyb(i,j) * wz(i,j,k) - uy(i,j,k) * wzb(i,j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nlterms(a,b,c)

      ! calculate the non linear terms of transport
      ! equations, in terms are ux, uy, uz, wx, wy and wz and  
      ! out terms are a, b and c
      !		a = uy * wx - ux * wy
      !		b = ux * wz - uz * wx
      !		c = uz * wy - uy * wz
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, j, k
      real*8    uxp(ptsx,jmax,kphys), wxp(ptsx,jmax,kphys),
     &          uyp(ptsx,jmax,kphys), wyp(ptsx,jmax,kphys),
     &          uzp(ptsx,jmax,kphys), wzp(ptsx,jmax,kphys),
     &           ap(ptsx,jmax,kphys),  bp(ptsx,jmax,kphys),
     &           cp(ptsx,jmax,kphys)
      complex*16  a(ptsx,jmax,kfour),   b(ptsx,jmax,kfour),
     &            c(ptsx,jmax,kfour)

      ! fft transforms from fourier to physical space
      call f_to_p(uxp,ux)
      call f_to_p(uyp,uy)
      call f_to_p(uzp,uz)
      call f_to_p(wxp,wx)
      call f_to_p(wyp,wy)
      call f_to_p(wzp,wz)

      ! nonlinear terms calculation in physical space
      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
           ap(i,j,k) = uyp(i,j,k) * wxp(i,j,k) - uxp(i,j,k) * wyp(i,j,k)
           bp(i,j,k) = uxp(i,j,k) * wzp(i,j,k) - uzp(i,j,k) * wxp(i,j,k)
           cp(i,j,k) = uzp(i,j,k) * wyp(i,j,k) - uyp(i,j,k) * wzp(i,j,k)
          end do
        end do
      end do

      ! fft back from physical to fourier space
      call p_to_f(ap,a)
      call p_to_f(bp,b)
      call p_to_f(cp,c)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine f_to_p(datap,dataf)

      ! transform from Fourier space to Physical space
      ! dataf-> Fourier space (in)
      ! datap -> Physical space (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp
      real*8 datap(ptsx,jmax,kphys)
      complex*16 dataf(ptsx,jmax,kfour)

      theta = -6.28318530717959d0/dble(kphys)
      c1    = 0.5d0
      wpr   = -2.0d0*dsin(0.5d0*theta)**2
      wpi   = dsin(theta)
      
      do j = 1, jmax
        do i = 1, ptsx
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
      subroutine p_to_f(datap,dataf)

      ! transform from Physical space to Fourier space
      ! datap -> Physical space data (in)
      ! dataf -> fourier space data (out)
      implicit none
      include 'par.for'
      integer p, p1, p2, p3, p4, n2p3, i, j, k
      real*8 c1, h1i, h1r, h2i, h2r, wis, wrs, theta,
     &       wi, wpi, wpr, wr, wtemp
      real*8 datap(ptsx,jmax,kphys)
      complex*16 dataf(ptsx,jmax,kfour)

      theta = 6.28318530717959d0/dble(kphys)
      c1    = 0.5d0
      wpr   = -2.0d0*dsin(0.5d0*theta)**2
      wpi   = dsin(theta)

      do j = 1, jmax
        do i = 1, ptsx
          call four1(datap,+1,i,j)
          wr   = 1.0d0 + wpr
          wi   = wpi
          n2p3 = kphys+3
          do p = 2, kphys/4
            p1            = 2*p-1
            p2            = p1+1
            p3            = n2p3-p2
            p4            = p3+1
            wrs           = wr
            wis           = wi
            h1r           = c1*(datap(i,j,p1)+datap(i,j,p3))
            h1i           = c1*(datap(i,j,p2)-datap(i,j,p4))
            h2r           = c1*(datap(i,j,p2)+datap(i,j,p4))
            h2i           = -c1*(datap(i,j,p1)-datap(i,j,p3))
            datap(i,j,p1) = h1r+wrs*h2r-wis*h2i
            datap(i,j,p2) = h1i+wrs*h2i+wis*h2r
            datap(i,j,p3) = h1r-wrs*h2r+wis*h2i
            datap(i,j,p4) = -h1i+wrs*h2i+wis*h2r
            wtemp         = wr
            wr            = wr*wpr-wi*wpi+wr
            wi            = wi*wpr+wtemp*wpi+wi
          end do
          h1r          = datap(i,j,1)
          datap(i,j,1) = (h1r+datap(i,j,2))/dble(kphys)
          datap(i,j,2) = h1r-datap(i,j,2)
          do k = 2, kphys
            datap(i,j,k) = 2.d0*datap(i,j,k)/dble(kphys)
          end do
          do k = 0, kfour - 1
            if (dabs(datap(i,j,2*k+1)).lt.1d-14) 
     &         datap(i,j,2*k+1) = 0.d0
            if (dabs(datap(i,j,2*k+2)).lt.1d-14) 
     &         datap(i,j,2*k+2) = 0.d0
            dataf(i,j,k+1)=dcmplx(datap(i,j,2*k+1),datap(i,j,2*k+2))
          end do
          dataf(i,j,1) = dcmplx(dreal(dataf(i,j,1)),0.d0)
          do k = 2 * kfour + 1, kphys
            datap(i,j,k) = 0.d0
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine four1(dataff,isig,ii,jj)

      implicit none
      include 'par.for'
      integer isig, i, istep, z, mm, mmax, ii, jj
      real*8 tempi, tempr, dataff(ptsx,jmax,kphys)
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
        mm = kphys/2
    1   if ((mm.ge.2).and.(z.gt.mm)) then
          z = z - mm
          mm = mm / 2
          goto 1
        endif
        z = z + mm
      end do 
      mmax = 2
    2 if (kphys.gt.mmax) then
        istep = 2 * mmax
        theta = 6.28318530717959d0/dble(isig*mmax)
        wpr   = -2.d0*dsin(0.5d0*theta)**2
        wpi   = dsin(theta)
        wr    = 1.d0
        wi    = 0.d0
        do mm = 1, mmax, 2
          do i = mm, kphys, istep
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
c         end of nonlinear terms calculation            c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
