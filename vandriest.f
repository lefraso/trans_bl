      program teste_vandriest

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8  S_xxp(ptsx,jmax,kphys), S_xyp(ptsx,jmax,kphys),
     &        S_xzp(ptsx,jmax,kphys), S_yyp(ptsx,jmax,kphys),
     &        S_yzp(ptsx,jmax,kphys), S_zzp(ptsx,jmax,kphys),
     &        nu_tp(ptsx,jmax,kphys), damp, cte, y, y_plus, ratio
      complex*16 S_xx(ptsx,jmax,kfour), S_xy(ptsx,jmax,kfour),
     &           S_xz(ptsx,jmax,kfour), S_yy(ptsx,jmax,kfour),
     &           S_yz(ptsx,jmax,kfour), S_zz(ptsx,jmax,kfour), 
     &           dvdx(ptsx,jmax,kfour)
      common/tensor/ S_xxp, S_xyp, S_xzp, S_yyp, S_yzp, S_zzp

      ! S_xx = dudx
      call derparx(S_xx, ux)
      ! S_yy = dvdy
      call deryfv(S_yy, uy)
      ! S_xy = dudy (variable economy)
      call dery(S_xy, ux)
      ! S_yz = dwdy (variable economy)
      call dery(S_yz, uz)
      call derparx(dvdx, uy)
      ! S_xz = dwdx (variable economy)
      call derparx(S_xz, uz)
      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            S_xy(i,j,k) = 0.5d0 * (S_xy(i,j,k) + dvdx(i,j,k))
            S_xz(i,j,k) = 0.5d0 * (S_xz(i,j,k) + v_kb(k) * ux(i,j,k))
            S_yz(i,j,k) = 0.5d0 * (S_yz(i,j,k) + v_kb(k) * uy(i,j,k))
            S_zz(i,j,k) = v_kb(k) * uz(i,j,k)
          end do
        end do
      end do

      call f_to_p(S_xxp,S_xx)
      call f_to_p(S_yyp,S_yy)
      call f_to_p(S_zzp,S_zz)
      call f_to_p(S_xyp,S_xy)
      call f_to_p(S_xzp,S_xz)
      call f_to_p(S_yzp,S_yz)

      do k = 1, kphys
        do j = 1, jmax
          if (stf.eq.1.d0) then
            y = dble(j-1) * dy
           else
            y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
          end if
          do i = 1, ptsx
            y_plus = dsqrt( Re * dreal(wz(i,1,1)) / dsqrt(fac_y) ) * y ! FUNCIONA PARA ESCOAMENTO 2D
            ratio  = min(y_plus/25.d0, 100.d0)
            ratio  = ratio ** 3
            damp   = (1.d0 - dexp(-ratio))
            cte    = C_s*C_s * delta_les*delta_les * damp
            nu_tp(i,j,k) = cte * 
     &                 dsqrt(2.d0 * ( S_xxp(i,j,k)*S_xxp(i,j,k) +
     &                                S_yyp(i,j,k)*S_yyp(i,j,k) +
     &                                S_zzp(i,j,k)*S_zzp(i,j,k) +
     &                        2.d0 * (S_xyp(i,j,k)*S_xyp(i,j,k) +
     &                                S_xzp(i,j,k)*S_xzp(i,j,k) +
     &                                S_yzp(i,j,k)*S_yzp(i,j,k))))
          end do
        end do
      end do

      return
      end


