ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c       subroutines transforms binary data in asc       c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program formated 12032009

      implicit none
      include 'par.for'
      character*15 nome
      integer i, j, k, t, my_rank, inter, shift
      real*8 x, y, z,
     &           uxp(imax,jmax,kphys), wxp(imax,jmax,kphys),
     &           uyp(imax,jmax,kphys), wyp(imax,jmax,kphys),
     &           uzp(imax,jmax,kphys), wzp(imax,jmax,kphys),
     &           thp(imax,jmax,kphys)
      complex*16 uxt(imax,jmax,kfour), wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour), wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour), wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour), 
     &            ux(ptsx,jmax,kfour),  wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),  wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),  wz(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour)
      real*8 uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax)


      ! Disturbances variables data
      inter = 2**( msh - 1 ) * ( stencil - 2 )
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(2,file=nome,form='unformatted')
        read(2) t
        if (my_form.eq.2) then
          read(2) ux,uy,uz,wx,wy,wz,th
         else
          read(2) ux,uy,uz,wx,wy,wz
        end if
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

      if (my_form.eq.0.or.my_form.eq.4) then
        open(1,file='baseflow2D/basens.bin',form='unformatted')
        read(1) uxbt, uybt, wzbt
        close(unit=1)
        do j = 1, jmax
          do i = 1, imax
            uxt(i,j,1) = uxt(i,j,1) - uxbt(i,j)
            uyt(i,j,1) = uyt(i,j,1) - uybt(i,j)
            wzt(i,j,1) = wzt(i,j,1) - wzbt(i,j)
          end do
        end do
      end if

      call f_to_p(uxp,uxt)
      call f_to_p(uyp,uyt)
      call f_to_p(uzp,uzt)
      call f_to_p(wxp,wxt)
      call f_to_p(wyp,wyt)
      call f_to_p(wzp,wzt)
      call f_to_p(thp,tht)

      if (my_form.eq.2) then
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'spatial.dat',status = 'unknown')
        write(3,*) 'VARIABLES="x","y","z","u","v","w",
     &             "wx","wy","wz","theta"'
        write(3,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,', F=POINT'
        
        z = -1.d0 * dz
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(3,6)x, y, z, uxp(i,j,kphys), uyp(i,j,kphys), 
     &                uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
     &                wzp(i,j,kphys), thp(i,j,kphys)
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
            do i = 1, imax
              x = x0 + dble(i-1) * dx
            write(3,6)x, y, z, uxp(i,j,k), uyp(i,j,k), uzp(i,j,k),
     &                wxp(i,j,k), wyp(i,j,k), wzp(i,j,k), thp(i,j,k)
            end do
          end do
        end do
        close (unit=3)

      else
        ! writes data to spacial space to be open by tecplot
        open (3, file = 'spatial.dat',status = 'unknown')
        write(3,*) 'VARIABLES="x","y","z","u","v","w","wx","wy","wz"'
        write(3,*) 'ZONE I=',imax,', J=',jmax,', K=',kphys+1,',F=POINT'
        
        z = -1.d0 * dz
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          do i = 1, imax
            x = x0 + dble(i-1) * dx
            write(3,5)x, y, z, uxp(i,j,kphys), uyp(i,j,kphys), 
     &                uzp(i,j,kphys), wxp(i,j,kphys), wyp(i,j,kphys),
     &                wzp(i,j,kphys)
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
            do i = 1, imax
              x = x0 + dble(i-1) * dx
            write(3,5)x, y, z, uxp(i,j,k), uyp(i,j,k), uzp(i,j,k),
     &                wxp(i,j,k), wyp(i,j,k), wzp(i,j,k)
            end do
          end do
        end do
        close (unit=3)
      end if

    5 format(1x,3d14.6,6d17.9)
    6 format(1x,3d14.6,7d17.9)

      stop
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
      integer isig, i, istep, z, mm, mmax, ii, jj
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
        mm=kphys/2
    1   if ((mm.ge.2).and.(z.gt.mm)) then
          z = z-mm
          mm = mm/2
          goto 1
        endif
        z = z+mm
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
c   end subroutines transforms binary data in asc       c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
