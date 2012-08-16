
 module fs

  use constants

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 contains

  subroutine general(fpp_fs, eta_end_fs, eta_adp, i)

   implicit none
   include '../comm.fs'
   real(kind=8), intent(inout) :: fpp_fs, eta_end_fs, eta_adp
   integer, intent(in)         :: i
   integer                     :: j, cont, N_eta
   real(kind=8)                :: y(6), w_new, w_zero_n, zk, zkpw, delta_w, x, E, h, m!, eta_adp

   h       = deta
   x       = fpp_fs                ! chute inicial
   E       = 1.d0
   cont    = 1

   do while (dabs(E).gt.gtol_fs)
 
    w_new = 100.d0
    N_eta = (eta_adp) / h + 1

    do while (dabs(w_new/x).gt.tol_fs)

     y(1) = 0.d0
     y(2) = 0.d0
     y(3) = x
     y(4) = 0.d0
     y(5) = 0.d0
     y(6) = 1.d0
 
     do j = 2 , N_eta
      call rk4(y, h, i)
     enddo

     w_new = ( - y(5) * y(2) - ( y(6) * y(3) ) + y(5) ) / ( y(5)**2 + y(6)**2 )
     x     = x + w_new

    enddo
 
    write(*,*) x
    E = ( 1.d0 - y(2) )**2 + y(3)**2
!   print *, cont, beta_fs, x, eta_adp, E
    print *, cont, beta_fs(i), x, eta_adp, E
!   eta_adp = eta_zero + eta_variation * cont
    eta_adp = eta_adp + eta_variation ! * cont
    cont    = cont + 1
 
   enddo

!  m          =  beta_fs / ( 2.d0 - beta_fs )         ! pag 164 do livro do Schlichting, m significa o parametro de adimensionalização pressão/gradiente.
   m          =  beta_fs(i) / ( 2.d0 - beta_fs(i) )   ! pag 164 do livro do Schlichting, m significa o parametro de adimensionalização pressão/gradiente.
   fpp_fs     = x
   eta_end_fs = eta_adp - eta_variation

  end subroutine general

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
  subroutine baseflow_fs(fpp_fs)

   implicit none
   include '../comm.fs'
   real(kind=8), intent(inout) :: fpp_fs(imax)
   real(kind=8)                :: y(6), h, m, ue, lx, xad, dya, eig(imax,jmax), eta, f, fp, fpp, dy
   real(kind=8)                :: ux(imax,jmax), uy(imax,jmax), wz(imax,jmax), uad, vad, vort, yad
   real(kind=8)                :: vf(jmax), vfp(jmax), vfpp(jmax), u_eta(jmax), deta0, x, yy
   real(kind=8)                :: a, b, c, d
   integer                     :: i, j, ind, flag, pt, k

   do i = 1 , imax
 
    m = beta_fs(i) / ( 2.d0 - beta_fs(i) )
 
    xad = x0 + dble(i-1) * dx
    ue  = xad ** m
    eta = 0.d0

    y(1)    = 0.d0
    y(2)    = 0.d0
    y(3)    = fpp_fs(i)
    y(4)    = 0.d0
    y(5)    = 0.d0
    y(6)    = 1.d0
    ux(i,1) = 0.d0
    uy(i,1) = 0.d0
    wz(i,1) = ue * dsqrt(0.5d0 * (m + 1.d0) * ue * Re / xad) * y(3)

    deta0 = dy0 * dsqrt(0.5d0 * Re * ue * (m + 1.d0) / xad)

    do j = 2 , jmax
     h  = deta0 * (stf**(j-2))
     pt = int( h  / 1.d-4 )
     h  = h / dble(pt)
     do k = 1 , pt
      call rk4(y, h, i)
     enddo
     if(stf.ne.1.d0) then
       eta = deta0 * (stf**(j-1) - 1.d0) / (stf - 1.d0)
     else
       eta = deta0 * (j-1)
     endif
     yad     = eta / dsqrt(0.5d0 * Re * ue * (m + 1.d0) / xad)
     uad     = ue * y(2)
     vad     = - dsqrt( 0.5d0 * ue* ( m + 1.d0 ) / (xad * Re)) * y(1) &
               - 0.5d0 * ( m - 1.d0 ) * yad * ue * y(2) / xad
     vort    = ue * dsqrt(0.5d0 * (m + 1.d0) * Re * ue / xad) * y(3) !                                   &
!            + 0.125d0 * (m - 1.d0)**2 * yad**2 * dsqrt(2.d0 * (m + 1.d0) / (Re * xad)) * y(3) / xad**2  &
!            + 0.25d0 * (3.d0 * m**2 - 4.d0 * m + 1.d0) * yad / (xad**2 * dsqrt(Re)) * y(2)              &
!            + 0.25d0 * (m - 1.d0) * dsqrt(2.d0 * (m + 1.d0) /(Re * xad**3)) * y(1)
     ux(i,j) = uad
     uy(i,j) = vad
     wz(i,j) = vort
!    write(*,"(1x,4d17.9)") eta, y(1), y(2), y(3)
!    write(3,"(1x,4d17.9)") yad, uad, vad, vort
    enddo
   enddo

   open(1,file='base_fs.bin', form='unformatted')
   write(1) stf
   write(1) ux, uy, wz
   close(unit=1)

   open (2, file = 'base_fs.dat',status = 'unknown')
   write(2,*) 'VARIABLES="x","y","velu","vely","vortz"'
   write(2,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

   do j = 1, jmax
     if(stf.eq.1.d0) then
      yy = dble(j-1) * dy0
     else
      yy = dy0 * (stf**(j-1) - 1.d0) / (stf - 1.d0)
     endif
     do i = 1, imax
       x = x0 + dble(i-1) * dx
       write(2,*)x,yy,ux(i,j),uy(i,j),wz(i,j)
     end do
   end do
   close(unit=2)
 
  end subroutine baseflow_fs
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
  subroutine delta_calculation(fpp_fs, eta_end_fs, delta, i)

   implicit none
   real(kind=8), intent(inout) :: fpp_fs, eta_end_fs, delta
   integer, intent(in)         :: i
   real(kind=8)                :: y(6), h
   real(kind=8)                :: ux, ux_new, variation, delta_eta
   integer                     :: j, Neta

   Neta      = 1500
   delta_eta = eta_end_fs / dble(Neta - 1)
   h         = delta_eta

   y(1) = 0.d0
   y(2) = 0.d0
   y(3) = fpp_fs
   y(4) = 0.d0
   y(5) = 0.d0
   y(6) = 1.d0

   ux = 0.d0

   do j = 2, Neta
    call rk4(y, h, i)
    if(y(2) .ge. 0.99d0) then
     delta = dble(j-1) * delta_eta
     write(*,*) '99% = ', delta, delta / dsqrt(2.d0)
     return
    endif
   enddo
 
  end subroutine delta_calculation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine rk4(y, h, i)
 
   implicit none
   real(kind=8), intent(inout) :: y(6)
   real(kind=8), intent(in)    :: h
   integer, intent(in)         :: i
   real(kind=8)                :: hh, h6, dydx(6), yt(6), dyt(6)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs(y, dydx, i)
   yt = y + hh * dydx           ! vector operation
 
   call derivs(yt, dyt, i)
   yt = y + hh * dyt            ! vector operation
   dydx = dydx + 2.d0 * dyt     ! vector operation
 
   call derivs(yt, dyt, i)
   yt = y + h * dyt             ! vector operation
   dydx = dydx + 2.d0 * dyt     ! vector operation
 
   call derivs(yt, dyt, i)
   y = y + h6 * ( dyt + dydx )  ! vector operation

  end subroutine rk4

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine derivs (y, dydx, i)

   implicit none
   include '../comm.fs'
   real(kind=8), intent(in)    :: y(6)
   real(kind=8), intent(out)   :: dydx(6)
   integer, intent(in)         :: i

   dydx(1) = y(2)
   dydx(2) = y(3)
   dydx(3) = ( - y(1) * y(3) ) - (beta_fs(i) * (1.d0 - y(2)**2))
   dydx(4) = y(5)
   dydx(5) = y(6)
   dydx(6) = ( - y(4) * y(3) - y(1) * y(6) ) + ( 2.d0 * beta_fs(i) * y(2) * y(5) )

  end subroutine derivs

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine rk42(y, h, i)
 
   implicit none
   real(kind=8), intent(inout) :: y(6)
   real(kind=8), intent(in)    :: h
   integer, intent(in)         :: i
   real(kind=8)                :: hh, h6, dydx(6), yt(6), dyt(6)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs2(y, dydx, i)
   yt = y + hh * dydx           ! vector operation
 
   call derivs2(yt, dyt, i)
   yt = y + hh * dyt            ! vector operation
   dydx = dydx + 2.d0 * dyt     ! vector operation
 
   call derivs2(yt, dyt, i)
   yt = y + h * dyt             ! vector operation
   dydx = dydx + 2.d0 * dyt     ! vector operation
 
   call derivs2(yt, dyt, i)
   y = y + h6 * ( dyt + dydx )  ! vector operation

  end subroutine rk42
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine derivs2 (y, dydx, i)

   implicit none
   include '../comm.fs'
   real(kind=8), intent(in)    :: y(6)
   real(kind=8), intent(out)   :: dydx(6)
   integer, intent(in)         :: i
   real(kind=8)                :: m

   m = beta_fs(i) / (2.d0 - beta_fs(i))   ! pag 164 do livro do Schlichting, m significa o parametro de adimensionalização pressão/gradiente.

   dydx(1) =   y(2)
   dydx(2) =   y(3)
   dydx(3) = - y(1) * y(3) * 0.5d0 * ( m + 1.d0 ) - m * ( 1.d0 - y(2) * y(2) )
   dydx(4) =   y(5)
   dydx(5) =   y(6)
   dydx(6) =   ( - y(4) * y(3) - y(1) * y(6)) + (2.d0 * beta_fs(i) * y(2) * y(5))

  end subroutine derivs2

 end module fs
