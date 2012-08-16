 program fs

  use constants
  use fs
  include '../comm.fs'
  real(kind = 8) :: fpp(imax), eta_end, eta_adp, delta, h
  integer        :: i

  ! reads beta_fs from a file
  open(1,file='../beta_fs.dist',form='formatted')
  read(1,*) beta_fs
  close(unit=1)

  i = 1
  fpp(i) = etapp_guess
  eta_adp = eta_zero
  call general(fpp(i), eta_end, eta_adp, i)
  do i = 2, imax
   fpp(i)  = fpp(i-1)
   eta_adp = eta_adp - 0.5d0
   call general(fpp(i), eta_end, eta_adp, i)
  enddo
 
! do i = 1, imax
!  call delta_calculation(fpp(i), eta_end, delta, i)
 
!  if(stf.eq.1.d0) then
!   h = 4.d0 * delta / dble(jmax-1)
!  else
!   h = (4.d0 * delta) * (stf - 1.d0) / ( stf**(jmax-1) - 1.d0 )
!  endif
! enddo

  call baseflow_fs(fpp)

 end program fs
