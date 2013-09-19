
      subroutine ns_residual(t)

      ! NS-residual residual calculation 

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include 'mpif.h'
      
      character*40 file_name
  
      integer :: i, j, t, it
      integer :: status(MPI_status_size)

      double precision :: uxwz(ptsx,jmax)               
      double precision :: uywz(ptsx,jmax)               
      double precision :: res(ptsx,jmax)               
 
      double precision :: duxwzdx(ptsx,jmax)               
      double precision :: duywzdy(ptsx,jmax)               
      double precision :: d2wzdx2(ptsx,jmax)               
      double precision :: d2wzdy2(ptsx,jmax)               
      
      double precision :: local_max(3), global_max(3), x, y 

      ! Non-linear calculation 
      do j = 1, jmax
       do i = 1, ptsx
        uxwz(i,j) = ux(i,j) * wz(i,j)
        uywz(i,j) = uy(i,j) * wz(i,j)
       enddo
      enddo  

      ! x-derivative calculation
      call derparx(duxwzdx,uxwz)  
      call derparxx(d2wzdx2,wz)  

      ! y-derivative calculation
      call dery(duywzdy,uywz)
      call deryy(d2wzdy2,wz)

      ! residual calculation  
      local_max = 0.d0
      do j = 1, jmax 
       do i = 1, ptsx 
         res(i,j) = duxwzdx(i,j) + duywzdy(i,j) + 
     &    ( - 1.d0 / Re ) * ( 
     &    d2wzdx2(i,j) + fac_y * d2wzdy2(i,j)          
     &    )
          if (local_max(3) .lt. abs(res(i,j))) then
           local_max(1) = i
           local_max(2) = j
           local_max(3) = abs(res(i,j)) 
          endif 
       enddo
      enddo
     
      if (my_rank.gt.0) then
        ! Send the maximum residual to the node 1
        call MPI_Send(local_max, 3, MPI_DOUBLE_PRECISION, 0, 51,
     &                MPI_COMM_WORLD, ierr)
       else
        global_max = local_max
        do it = 1, numproc
          ! Receive residual from all processing elements and find the global maximum
          call MPI_Recv(local_max, 3, MPI_DOUBLE_PRECISION, it, 51,
     &                  MPI_COMM_WORLD, status, ierr)
          if (local_max(3) .gt. global_max(3)) then
           global_max = local_max  
          endif        
        end do
        write(*,*) t, int(global_max(1)), int(global_max(2)),
     &             global_max(3) 
      end if

      if ( mod(t,1000).eq.0 ) then
      write(file_name,'(a,i0.5,a,i0.2,a)')
     &                      'residual_',t,'_',my_rank,'.dat'
      open (3, file = file_name, status = 'unknown')
      write(3,*) 'VARIABLES="x","y","residual"'
      write(3,*) 'ZONE I=',ptsx,', J=',jmax,', F=POINT'
      do j = 1, jmax
       if(stf.eq.1.d0) then
        y = dble(j-1) * dy
       else
        y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)
       endif
       y = y / sqrt(fac_y)
       do i = 1, ptsx
         x = x0 + dble(i-1+shift) * dx
         write(3,*) x, y, abs(res(i,j))
       end do
      end do
      close(3)
      endif

      end subroutine ns_residual 
