ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c         LES calculation                               c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine les_terms(F_nux, F_nuy, F_nuz)
      
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8 S_xxp(ptsx,jmax,kphys),    S_xyp(ptsx,jmax,kphys),
     &       S_xzp(ptsx,jmax,kphys),    S_yyp(ptsx,jmax,kphys),
     &       S_yzp(ptsx,jmax,kphys),    S_zzp(ptsx,jmax,kphys), 
     &       nu_tp(ptsx,jmax,kphys), dnu_tdxp(ptsx,jmax,kphys),
     &    dnu_tdyp(ptsx,jmax,kphys), dnu_tdzp(ptsx,jmax,kphys), 
     &       lapup(ptsx,jmax,kphys),    lapvp(ptsx,jmax,kphys), 
     &       lapwp(ptsx,jmax,kphys),   F_nuxp(ptsx,jmax,kphys), 
     &      F_nuyp(ptsx,jmax,kphys),   F_nuzp(ptsx,jmax,kphys),
     &      auxp(ptsx,jmax,kphys)
      complex*16 d2udx2(ptsx,jmax,kfour),  d2udy2(ptsx,jmax,kfour),
     &           d2vdx2(ptsx,jmax,kfour),  d2vdy2(ptsx,jmax,kfour), 
     &           d2wdx2(ptsx,jmax,kfour),  d2wdy2(ptsx,jmax,kfour),
     &             lapu(ptsx,jmax,kfour),    lapv(ptsx,jmax,kfour), 
     &             lapw(ptsx,jmax,kfour), 
     &          dnu_tdy(ptsx,jmax,kfour), dnu_tdz(ptsx,jmax,kfour),
     &            F_nux(ptsx,jmax,kfour),   F_nuy(ptsx,jmax,kfour),
     &            F_nuz(ptsx,jmax,kfour),   nu_t(ptsx,jmax,kfour)
      common/tensor/ S_xxp, S_xyp, S_xzp, S_yyp, S_yzp, S_zzp

      call calc_Sij

      call calc_nu_t_smag(nu_tp)
c     call calc_nu_t_wale(nu_tp)

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            auxp(i,j,k) = nu_tp(i,j,k)
          end do
        end do
      end do
      call p_to_f(auxp, nu_t) ! converts nu_tp from Physical space to nu_t in Fourier space
      call derparxr3d(dnu_tdxp, nu_tp)
      call deryr3d(dnu_tdyp, nu_tp)

      call derparxx(d2udx2, ux)
      call derparxx(d2vdx2, uy)
      call derparxx(d2wdx2, uz)
      call deryy(d2udy2, ux)
      call deryy(d2vdy2, uy)
      call deryy(d2wdy2, uz)

      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            lapu(i,j,k)    = d2udx2(i,j,k) + d2udy2(i,j,k) +
     &                               v_k2b2(k) * ux(i,j,k)
            lapv(i,j,k)    = d2vdx2(i,j,k) + d2vdy2(i,j,k) +
     &                               v_k2b2(k) * uy(i,j,k)
            lapw(i,j,k)    = d2wdx2(i,j,k) + d2wdy2(i,j,k) +
     &                               v_k2b2(k) * uz(i,j,k)
            dnu_tdz(i,j,k) = v_kb(k) * nu_t(i,j,k)
          end do
        end do
      end do
       
      call f_to_p(lapup, lapu)
      call f_to_p(lapvp, lapv)
      call f_to_p(lapwp, lapw)
      call f_to_p(dnu_tdzp, dnu_tdz)

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            F_nuxp(i,j,k) =  nu_tp(i,j,k) *    lapup(i,j,k)
!    &             + 2.d0 * (S_xxp(i,j,k) * dnu_tdxp(i,j,k)
!    &             +         S_xyp(i,j,k) * dnu_tdyp(i,j,k)
!    &             +         S_xzp(i,j,k) * dnu_tdzp(i,j,k) )
            F_nuyp(i,j,k) =  nu_tp(i,j,k) *    lapvp(i,j,k)
!    &             + 2.d0 * (S_xyp(i,j,k) * dnu_tdxp(i,j,k)
!    &             +         S_yyp(i,j,k) * dnu_tdyp(i,j,k)
!    &             +         S_yzp(i,j,k) * dnu_tdzp(i,j,k) )
            F_nuzp(i,j,k) =  nu_tp(i,j,k) *    lapwp(i,j,k)
!    &             + 2.d0 * (S_xzp(i,j,k) * dnu_tdxp(i,j,k)
!    &             +         S_yzp(i,j,k) * dnu_tdyp(i,j,k)
!    &             +         S_zzp(i,j,k) * dnu_tdzp(i,j,k) )
          end do
        end do
      end do

      call p_to_f(F_nuxp, F_nux)
      call p_to_f(F_nuyp, F_nuy)
      call p_to_f(F_nuzp, F_nuz)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_Sij

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8  S_xxp(ptsx,jmax,kphys), S_xyp(ptsx,jmax,kphys),
     &        S_xzp(ptsx,jmax,kphys), S_yyp(ptsx,jmax,kphys),
     &        S_yzp(ptsx,jmax,kphys), S_zzp(ptsx,jmax,kphys),
     &        dudxp(ptsx,jmax,kphys), dudyp(ptsx,jmax,kphys), 
     &        dudzp(ptsx,jmax,kphys), dvdxp(ptsx,jmax,kphys),
     &        dvdyp(ptsx,jmax,kphys), dvdzp(ptsx,jmax,kphys),
     &        dwdxp(ptsx,jmax,kphys), dwdyp(ptsx,jmax,kphys),
     &        dwdzp(ptsx,jmax,kphys)
      complex*16 dud(ptsx,jmax,kfour), dvd(ptsx,jmax,kfour), 
     &           dwd(ptsx,jmax,kfour)
      common/tensor/ S_xxp, S_xyp, S_xzp, S_yyp, S_yzp, S_zzp
      common/derivatives/ dudxp, dvdxp, dwdxp, dudyp, dvdyp, dwdyp, 
     &                    dudzp, dvdzp, dwdzp

      call derparx(dud, ux)
      call derparx(dvd, uy)
      call derparx(dwd, uz)

      call f_to_p(dudxp,dud)
      call f_to_p(dvdxp,dvd)
      call f_to_p(dwdxp,dwd)

      call dery(dud, ux)
      call deryfv(dvd, uy)
      call dery(dwd, uz)

      call f_to_p(dudyp,dud)
      call f_to_p(dvdyp,dvd)
      call f_to_p(dwdyp,dwd)

      do k = 1, kfour
        do j = 1, jmax
          do i = 1, ptsx
            dud(i,j,k) = v_kb(k) * ux(i,j,k)
            dvd(i,j,k) = v_kb(k) * uy(i,j,k)
            dwd(i,j,k) = v_kb(k) * uz(i,j,k)
          end do
        end do
      end do

      call f_to_p(dudzp,dud)
      call f_to_p(dvdzp,dvd)
      call f_to_p(dwdzp,dwd)

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
            S_xxp(i,j,k) = dudxp(i,j,k)
            S_yyp(i,j,k) = dvdyp(i,j,k)
            S_zzp(i,j,k) = dwdzp(i,j,k)
            S_xyp(i,j,k) = 0.5d0 * (dudyp(i,j,k) + dvdxp(i,j,k))
            S_xzp(i,j,k) = 0.5d0 * (dudzp(i,j,k) + dwdxp(i,j,k))
            S_yzp(i,j,k) = 0.5d0 * (dvdzp(i,j,k) + dwdyp(i,j,k))
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_nu_t_smag(nu_tp)

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
      complex*16 dvdx(ptsx,jmax,kfour)
      common/tensor/ S_xxp, S_xyp, S_xzp, S_yyp, S_yzp, S_zzp

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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_nu_t_wale(nu_tp)

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8 dudxp(ptsx,jmax,kphys), dudyp(ptsx,jmax,kphys), 
     &       dudzp(ptsx,jmax,kphys), dvdxp(ptsx,jmax,kphys),
     &       dvdyp(ptsx,jmax,kphys), dvdzp(ptsx,jmax,kphys),
     &       dwdxp(ptsx,jmax,kphys), dwdyp(ptsx,jmax,kphys),
     &       dwdzp(ptsx,jmax,kphys), nu_tp(ptsx,jmax,kphys), 
     &       S_xxp(ptsx,jmax,kphys), S_xyp(ptsx,jmax,kphys),
     &       S_xzp(ptsx,jmax,kphys), S_yyp(ptsx,jmax,kphys),
     &       S_yzp(ptsx,jmax,kphys), S_zzp(ptsx,jmax,kphys), 
     &       var, ss, eqnA, s11d, s22d, s33d, s12d, 
     &       s13d, s23d, sdsd
      common/tensor/ S_xxp, S_xyp, S_xzp, S_yyp, S_yzp, S_zzp
      common/derivatives/ dudxp, dvdxp, dwdxp, dudyp, dvdyp, dwdyp, 
     &                    dudzp, dvdzp, dwdzp

      do k = 1, kphys
        do j = 1, jmax
          do i = 1, ptsx
        
            var = ( dudxp(i,j,k)*dudxp(i,j,k) + 
     &              dvdyp(i,j,k)*dvdyp(i,j,k) + 
     &              dwdzp(i,j,k)*dwdzp(i,j,k) )

            ss  = ( var + 
     &          2.d0*S_xyp(i,j,k)*S_xyp(i,j,k)+
     &          2.d0*S_xzp(i,j,k)*S_xzp(i,j,k)+
     &          2.d0*S_yzp(i,j,k)*S_yzp(i,j,k))

            eqnA = var +
     &        2.d0*dudyp(i,j,k)*dvdxp(i,j,k) +
     &        2.d0*dudzp(i,j,k)*dwdxp(i,j,k) +
     &        2.d0*dvdzp(i,j,k)*dwdyp(i,j,k)
            s11d = dudxp(i,j,k)*dudxp(i,j,k) +
     &             dudyp(i,j,k)*dvdxp(i,j,k) +
     &             dudzp(i,j,k)*dwdxp(i,j,k) - eqnA/3.d0
            s22d = dvdxp(i,j,k)*dudyp(i,j,k) +
     &             dvdyp(i,j,k)*dvdyp(i,j,k) + 
     &             dvdzp(i,j,k)*dwdyp(i,j,k) - eqnA/3.d0
            s33d = dwdxp(i,j,k)*dudzp(i,j,k) +
     &             dwdyp(i,j,k)*dvdzp(i,j,k) +
     &             dwdzp(i,j,k)*dwdzp(i,j,k) - eqnA/3.d0
            s12d =0.5d0*( dudxp(i,j,k)*dudyp(i,j,k) + 
     &                    dudyp(i,j,k)*dvdyp(i,j,k) + 
     &                    dudzp(i,j,k)*dwdyp(i,j,k) +
     &                    dvdxp(i,j,k)*dudxp(i,j,k) +
     &                    dvdyp(i,j,k)*dvdxp(i,j,k) +
     &                    dvdzp(i,j,k)*dwdxp(i,j,k) )
            s13d =0.5d0*( dudxp(i,j,k)*dudzp(i,j,k) +
     &                    dudyp(i,j,k)*dvdzp(i,j,k) +
     &                    dudzp(i,j,k)*dwdzp(i,j,k) +
     &                    dwdxp(i,j,k)*dudxp(i,j,k) +
     &                    dwdyp(i,j,k)*dvdxp(i,j,k) +
     &                    dwdzp(i,j,k)*dwdxp(i,j,k) )
            s23d =0.5d0*( dvdxp(i,j,k)*dudzp(i,j,k) +
     &                    dvdyp(i,j,k)*dvdzp(i,j,k) +
     &                    dvdzp(i,j,k)*dwdzp(i,j,k) +
     &                    dwdxp(i,j,k)*dudyp(i,j,k) +
     &                    dwdyp(i,j,k)*dvdyp(i,j,k) +
     &                    dwdzp(i,j,k)*dwdyp(i,j,k) )
            sdsd = (  s11d*s11d +      s22d*s22d +      s33d*s33d  +
     &           2.d0*s12d*s12d + 2.d0*s13d*s13d + 2.d0*s23d*s23d  )
            nu_tp(i,j,k) = C_w * C_w * delta_les * delta_les * 
     &                     sdsd**1.5 / (ss**2.5+sdsd**1.25)
c           if (nu_tp(i,j,k).lt.0.d0) nu_tp(i,j,k) = 0.d0
          end do
        end do
      end do

      return
      end
