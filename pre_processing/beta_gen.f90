 program beta_gen

  use constants
! use fs
  real(kind = 8) :: x, beta_fs
  integer        :: i

  ! writes beta_fs in a file
  open (1, file = '../beta_fs.dist',status = 'unknown')
  do i = 1, imax
   x = i * dx
!  beta_fs = x**2 - 5 * x + 4
   beta_fs = 0.d0
!  write(1,*) x, beta_fs
   write(1,*) beta_fs
  enddo
  close(unit=1)
  
 end program beta_gen

