 program fs

! use constants
  use fs
  real(kind = 8) :: fpp(imax), thp(imax), eta_adp
  integer        :: i

  ! reads beta_fs from a file
  open(1,file='../beta_fs.dist',form='formatted')
  read(1,*) beta_fs
  close(unit=1)

  fpp(1)  = etapp_guess
  eta_adp = eta_zero
  call fpp_finder(fpp(1), eta_adp, 1)
  do i = 2, imax
   fpp(i)  = fpp(i-1)
   eta_adp = eta_adp - 0.25d0
   call fpp_finder(fpp(i), eta_adp, i)
   call delta_calculation(fpp(i), eta_adp, i)
  enddo
 
  if (my_form.eq.2) then
    do i = 1, imax
      call thetap_finder(thp(i), fpp(i), i)
    enddo
  endif

  call baseflow_fs(fpp, thp)

 end program fs
