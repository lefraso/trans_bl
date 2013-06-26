cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c                main subroutines                                c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program trans_bl_20130503

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'mpif.h'
      logical run
      integer p
      real*8 tempo_inicial, tempo_final

      ! Start up MPI
      call MPI_Init(ierr)

      ! Find out process rank
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

      ! Find out number of processes
      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

      ! p nodes in the environment. So the number of the last node is numproc -->  
      numproc = p - 1

      ! calculates the intersection
      ! inter = 2**(number of meshes used in multigrid solver - 1) * (stencil-2)
      inter = 2**( msh - 1 ) * ( stencil - 2 )

      ! calculate the i_shift from one computer to another wiht the domain
      ! decomposition
      shift = my_rank * (ptsx - inter - 1)

      ! to calculate the total time simulation
      if (my_rank.eq.0) then
        tempo_inicial = mpi_wtime() ! calcula o tempo inicial
      end if

      ! verification of adopted points in x and y directions
      run = .true.
      call verifica(run)
      if (run) call solver

      ! to calculate the total time simulation
      if (my_rank.eq.0) then
        tempo_final = mpi_wtime() !calcula o tempo final
        tempo_final = tempo_final - tempo_inicial
        write(*,4) tempo_final
      end if
   4  format(1x,1d17.9)

      call MPI_Finalize(ierr)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine verifica(run)

      implicit none
      include 'par.for'
      include 'comm.par'
      logical run
      integer nodesx, chute1, chute2
      real*8 n_npr

      ! calculates the value of points in the x direction of each node
      nodesx = ( imax + (inter + 1) * numproc ) / (numproc+1)
      n_npr  = dble( imax + (inter + 1) * numproc ) / dble(numproc+1)

      ! Stop the program if the number of total points is not exactly
      ! divided by nodes and if the  number of points in the x
      ! direction can not be used in multigrid program and if the
      ! number of points in the y direction can not be used in
      ! multigrid program

      if (numproc.ne.np-1) then
        write(*,*) 'Numero de processadores e incompativel'
        run = .false.
      end if

      if (ptsx.ne.nodesx) then
        write(*,*) 'ALTERAR O VALOR DE PTSX NO ARQUIVO par.for PARA:'
        write(*,*) nodesx 
        run = .false.
      end if

      if ((nodesx.ne.n_npr).or.(mod(nodesx-1,2**(msh-1)).ne.0)) then
        chute1 = ( (nodesx-1) / 2**(msh-1) ) * 2**(msh-1) + 1
        chute2 = ( ( (nodesx-1) / 2**(msh-1) ) + 1 ) * 2**(msh-1) + 1
        write(*,*)' Number of points in the x direction can
     &              not be used in multigrid solver'
        write(*,*)' The number of points in the x direction
     &              should be:'
        write(*,*)  (chute1 - 1) * (numproc+1) - numproc
     &              * inter + 1,'   or'
        write(*,*)  (chute2 - 1) * (numproc+1) - numproc
     &              * inter + 1
        run = .false.
      end if

      if (nodesx .lt. 2*inter) then
        write(*,*) ' Number of points in the x dirextion can note be
     &               used in multigrid solver'
        write(*,*) ' The number of points in the x direction should
     &               be:'
        write(*,*) ((nodesx * 2) - 1) * (numproc+1) - numproc
     &             * inter + 2
        run = .false.
      end if

      if((dble(jmax-1)/dble(2**(msh-1))).ne.((jmax-1)/2**(msh-1))) then
        write(*,*)' Number of points in the y direction can
     &              not be used in multigrid solver'
        write(*,*)' The number of points in the y direction
     &              should be ', 2**(msh-1) *
     &              ((jmax - 1) / 2**(msh-1)) + 1,' or ',
     &              2**(msh-1) * (((jmax - 1) / 2**(msh-1)) + 1) + 1
        run = .false.
      end if

      if (jmax.lt.2**msh) then
        write(*,*) ' Number of points in the y direction can note be
     & used in multigrid solver'
        write(*,*) ' The number of points in the y direction should
     & be:'
        write(*,*) jmax * 2 - 1
        run = .false.
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solver

      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'mpif.h'
      integer i, j, k, t, fanal, t0, i_ini
      real*8 dt2, dt6
      complex*16 dv1x(ptsx,jmax,kfour), dv2x(ptsx,jmax,kfour),
     &           dv1y(ptsx,jmax,kfour), dv2y(ptsx,jmax,kfour),
     &           dv1z(ptsx,jmax,kfour), dv2z(ptsx,jmax,kfour),
     &           dv1t(ptsx,jmax,kfour), dv2t(ptsx,jmax,kfour),
     &            wx1(ptsx,jmax,kfour),  wy1(ptsx,jmax,kfour),
     &            wz1(ptsx,jmax,kfour),  th1(ptsx,jmax,kfour)

      i_ini = 1
      if (my_rank.eq.0) i_ini = 3
      select case (my_form)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! CASE 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       case(0)
         call init_fi(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, wx1, wy1, 
     &                wz1, dt2, dt6, t0)
        
        fanal = stpp / 16
        
        do t = t0, tt + 15 * fanal
c       do t = t0, tt
        
          ! first Runge-Kutta step
          call drv_fi(dv1x, dv1y, dv1z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx1(i,j,k) = wx(i,j,k)
                wy1(i,j,k) = wy(i,j,k)
                wz1(i,j,k) = wz(i,j,k)
                wx(i,j,k)  = wx1(i,j,k) + dv1x(i,j,k) * dt2
                wy(i,j,k)  = wy1(i,j,k) + dv1y(i,j,k) * dt2
                wz(i,j,k)  = wz1(i,j,k) + dv1z(i,j,k) * dt2
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
            call ts2d_pert(t, 0.5d0)
c           call ts3d_pert(t, 0.5d0)
          end if
          call loop(1d-5)
        
          ! second Runge-Kutta step
          call drv_fi(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt2
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt2
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt2
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          call loop(1d-5)
        
          ! third Runge-Kutta step
          call drv_fi(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
            call ts2d_pert(t, 1.d0)
c           call ts3d_pert(t, 1.d0)
          end if
          call loop(1d-5)
        
          ! fourth Runge-Kutta step
          call drv_fi(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k) = wx1(i,j,k)
     &                    + dt6 * ( dv1x(i,j,k) + dv2x(i,j,k) )
                wy(i,j,k) = wy1(i,j,k)
     &                    + dt6 * ( dv1y(i,j,k) + dv2y(i,j,k) )
                wz(i,j,k) = wz1(i,j,k)
     &                    + dt6 * ( dv1z(i,j,k) + dv2z(i,j,k) )
      !!!!!! filtro para eliminação das oscilações de Gibbs
c               wx(i,j,k) =  wx1(i,j,k) + gibbs(k)
c    &                    *  dt6 * ( dv1x(i,j,k) + dv2x(i,j,k) )
c               wy(i,j,k) =  wy1(i,j,k) + gibbs(k)
c    &                    *  dt6 * ( dv1y(i,j,k) + dv2y(i,j,k) )
c               wz(i,j,k) =  wz1(i,j,k) + gibbs(k)
c    &                    *  dt6 * ( dv1z(i,j,k) + dv2z(i,j,k) )
              end do
            end do
          end do
          call loop(1d-6)

          write(*,*)my_rank, t, ux(ptsx,jmax/2,2)
          if (mod(t,stpp).eq.0) call escreve(t)
        
          if(t.ge.tt) call escreve2(t,fanal)
        
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! CASE 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       case(1)
         call init_gv(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, wx1, wy1, 
     &                wz1, dt2, dt6, t0)
        
        fanal = stpp / 16
        
        do t = t0, tt + 15 * fanal
c       do t = t0, tt
        
          ! first Runge-Kutta step
          call drv_gv(dv1x, dv1y, dv1z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx1(i,j,k) = wx(i,j,k)
                wy1(i,j,k) = wy(i,j,k)
                wz1(i,j,k) = wz(i,j,k)
                wx(i,j,k)  = wx1(i,j,k) + dv1x(i,j,k) * dt2
                wy(i,j,k)  = wy1(i,j,k) + dv1y(i,j,k) * dt2
                wz(i,j,k)  = wz1(i,j,k) + dv1z(i,j,k) * dt2
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
c           call gv_pert(t)
            call ts2d_pert(t, 0.5d0)
c           call ts3d_pert(t, 0.5d0)
          end if
          call loop(1d-5)
        
          ! second Runge-Kutta step
          call drv_gv(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt2
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt2
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt2
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          call loop(1d-5)
        
          ! third Runge-Kutta step
          call drv_gv(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
c           call gv_pert(t)
            call ts2d_pert(t, 1.d0)
c           call ts3d_pert(t, 1.d0)
          end if
          call loop(1d-5)
        
          ! fourth Runge-Kutta step
          call drv_gv(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k) = wx1(i,j,k)
     &                    + dt6 * ( dv1x(i,j,k) + dv2x(i,j,k) )
                wy(i,j,k) = wy1(i,j,k)
     &                    + dt6 * ( dv1y(i,j,k) + dv2y(i,j,k) )
                wz(i,j,k) = wz1(i,j,k)
     &                    + dt6 * ( dv1z(i,j,k) + dv2z(i,j,k) )
              end do
            end do
          end do
          call loop(1d-6)
        
          write(*,*)my_rank, t, ux(ptsx,jmax/2,2)
          if (mod(t,stpp).eq.0) call escreve(t)
        
          if(t.ge.tt) call escreve2(t,fanal)
        
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! CASE 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       case(2)
         call init_gvth(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, dv1t, dv2t, 
     &                  wx1, wy1, wz1, th1, dt2, dt6, t0)
        
        fanal = stpp / 16
        
        do t = t0, tt + 15 * fanal
c       do t = t0, tt
        
          ! first Runge-Kutta step
          call drv_gvth(dv1x, dv1y, dv1z, dv1t)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx1(i,j,k) = wx(i,j,k)
                wy1(i,j,k) = wy(i,j,k)
                wz1(i,j,k) = wz(i,j,k)
                th1(i,j,k) = th(i,j,k)
                wx(i,j,k)  = wx1(i,j,k) + dv1x(i,j,k) * dt2
                wy(i,j,k)  = wy1(i,j,k) + dv1y(i,j,k) * dt2
                wz(i,j,k)  = wz1(i,j,k) + dv1z(i,j,k) * dt2
                th(i,j,k)  = th1(i,j,k) + dv1t(i,j,k) * dt2
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
c           call gv_pert(t)
            call ts2d_pert(t, 0.5d0)
c           call ts3d_pert(t, 0.5d0)
          end if
          call loop(1d-5)
        
          ! second Runge-Kutta step
          call drv_gvth(dv2x, dv2y, dv2z, dv2t)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt2
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt2
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt2
                th(i,j,k)   = th1(i,j,k) + dv2t(i,j,k) * dt2
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
                dv1t(i,j,k) = dv1t(i,j,k) + 2.d0 * dv2t(i,j,k)
              end do
            end do
          end do
          call loop(1d-5)
        
          ! third Runge-Kutta step
          call drv_gvth(dv2x, dv2y, dv2z, dv2t)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt
                th(i,j,k)   = th1(i,j,k) + dv2t(i,j,k) * dt
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
                dv1t(i,j,k) = dv1t(i,j,k) + 2.d0 * dv2t(i,j,k)
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
c           call gv_pert(t)
            call ts2d_pert(t, 1.d0)
c           call ts3d_pert(t, 1.d0)
          end if
          call loop(1d-5)
        
          ! fourth Runge-Kutta step
          call drv_gvth(dv2x, dv2y, dv2z, dv2t)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k) = wx1(i,j,k)
     &                    + dt6 * ( dv1x(i,j,k) + dv2x(i,j,k) )
                wy(i,j,k) = wy1(i,j,k)
     &                    + dt6 * ( dv1y(i,j,k) + dv2y(i,j,k) )
                wz(i,j,k) = wz1(i,j,k)
     &                    + dt6 * ( dv1z(i,j,k) + dv2z(i,j,k) )
                th(i,j,k) = th1(i,j,k)
     &                    + dt6 * ( dv1t(i,j,k) + dv2t(i,j,k) )
              end do
            end do
          end do
          call loop(1d-6)
        
          write(*,*)my_rank, t, ux(ptsx,jmax/2,2)
          if (mod(t,stpp).eq.0) call escreve(t)
        
          if(t.ge.tt) call escreve2(t,fanal)
        
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! CASE 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       case(4)
         call init_fi(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, wx1, wy1, 
     &                wz1, dt2, dt6, t0)
        
        fanal = stpp / 16
        
        do t = t0, tt + 15 * fanal
c       do t = t0, tt
        
          ! first Runge-Kutta step
          call drv_les(dv1x, dv1y, dv1z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx1(i,j,k) = wx(i,j,k)
                wy1(i,j,k) = wy(i,j,k)
                wz1(i,j,k) = wz(i,j,k)
                wx(i,j,k)  = wx1(i,j,k) + dv1x(i,j,k) * dt2
                wy(i,j,k)  = wy1(i,j,k) + dv1y(i,j,k) * dt2
                wz(i,j,k)  = wz1(i,j,k) + dv1z(i,j,k) * dt2
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
            call ts2d_pert(t, 0.5d0)
c           call ts3d_pert(t, 0.5d0)
          end if
          call loop(1d-5)
        
          ! second Runge-Kutta step
          call drv_les(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt2
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt2
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt2
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          call loop(1d-5)
        
          ! third Runge-Kutta step
          call drv_les(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k)   = wx1(i,j,k) + dv2x(i,j,k) * dt
                wy(i,j,k)   = wy1(i,j,k) + dv2y(i,j,k) * dt
                wz(i,j,k)   = wz1(i,j,k) + dv2z(i,j,k) * dt
                dv1x(i,j,k) = dv1x(i,j,k) + 2.d0 * dv2x(i,j,k)
                dv1y(i,j,k) = dv1y(i,j,k) + 2.d0 * dv2y(i,j,k)
                dv1z(i,j,k) = dv1z(i,j,k) + 2.d0 * dv2z(i,j,k)
              end do
            end do
          end do
          ! disturbance introductions
          if (my_rank.eq.0) then
            call ts2d_pert(t, 1.d0)
c           call ts3d_pert(t, 1.d0)
          end if
          call loop(1d-5)
        
          ! fourth Runge-Kutta step
          call drv_les(dv2x, dv2y, dv2z)
          do k = 1, kfour
            do j = 2, jmax - 1
              do i = i_ini, ptsx
                wx(i,j,k) = wx1(i,j,k)
     &                    + dt6 * ( dv1x(i,j,k) + dv2x(i,j,k) )
                wy(i,j,k) = wy1(i,j,k)
     &                    + dt6 * ( dv1y(i,j,k) + dv2y(i,j,k) )
                wz(i,j,k) = wz1(i,j,k)
     &                    + dt6 * ( dv1z(i,j,k) + dv2z(i,j,k) )
              end do
            end do
          end do
          call loop(1d-6)
        
          write(*,*)my_rank, t, ux(ptsx,jmax/2,2)
          if (mod(t,stpp).eq.0) call escreve(t)
        
          if(t.ge.tt) call escreve2(t,fanal)
        
        end do

      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_fi(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, wx1, wy1, 
     &                   wz1, dt2, dt6, t0)

      ! initialize the program main variables
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.coef'
      include 'comm.multi'
      include 'mpif.h'
      include 'comm.fs'
      character*23 nome
      integer i, j, k, t0
      real*8 dt2, dt6, a(imax,5), luf(imax,5), bdfc(2,imax), ep,
     &       uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax),
     &       afil(ptsx), bfil(ptsx), cfil(ptsx), stf_v, beta_fs_v,
     &       ueptsx(ptsx)
      complex*16 dv1x(ptsx,jmax,kfour), dv2x(ptsx,jmax,kfour),
     &           dv1y(ptsx,jmax,kfour), dv2y(ptsx,jmax,kfour),
     &           dv1z(ptsx,jmax,kfour), dv2z(ptsx,jmax,kfour),
     &            wx1(ptsx,jmax,kfour),  wy1(ptsx,jmax,kfour),
     &            wz1(ptsx,jmax,kfour),   d2uydx2(ptsx,kfour),
     &             fx(ptsx,jmax,kfour),   fy(ptsx,jmax,kfour),
     &             fz(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb
      common/derw/ d2uydx2
      common/fil/ luf
      common/force/ fx, fy, fz
      common/bd/ bdfc
      common/filt/ afil, bfil, cfil

      t0  = 1
      dt2 = 0.5d0 * dt
      dt6 = dt / 6.d0

      ! mounts the lu matrix for the filter
      call lhsf(a)
      call ludecomp(a, luf)
      call lhs_tridf(afil,bfil,cfil)

      ! x function of the disturbance strip for TS disturbances
      call var_ts

!!    Gibbs filter for filtering Fourier functions
      call gibbs_function

      ! calculates the delta function for the immersed boundary method
      call delta

      ! all the variables are set to zero
      do k = 1, kfour
        do i = 1, ptsx
          d2uydx2(i,k) = dcmplx(0.d0,0.d0)
          do j = 1, jmax
            ux(i,j,k)   = dcmplx(0.d0,0.d0)
            uy(i,j,k)   = dcmplx(0.d0,0.d0)
            uz(i,j,k)   = dcmplx(0.d0,0.d0)
            wx(i,j,k)   = dcmplx(0.d0,0.d0)
            wy(i,j,k)   = dcmplx(0.d0,0.d0)
            wz(i,j,k)   = dcmplx(0.d0,0.d0)
            wx1(i,j,k)  = dcmplx(0.d0,0.d0)
            wy1(i,j,k)  = dcmplx(0.d0,0.d0)
            wz1(i,j,k)  = dcmplx(0.d0,0.d0)
            dv1x(i,j,k) = dcmplx(0.d0,0.d0)
            dv2x(i,j,k) = dcmplx(0.d0,0.d0)
            dv1y(i,j,k) = dcmplx(0.d0,0.d0)
            dv2y(i,j,k) = dcmplx(0.d0,0.d0)
            dv1z(i,j,k) = dcmplx(0.d0,0.d0)
            dv2z(i,j,k) = dcmplx(0.d0,0.d0)
            fx(i,j,k)   = dcmplx(0.d0,0.d0)
            fy(i,j,k)   = dcmplx(0.d0,0.d0)
            fz(i,j,k)   = dcmplx(0.d0,0.d0)
          end do
        end do
      end do

      ! reads the boundary layer profile
      write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
      open (1, file = nome,form = 'unformatted')
      read(1) uxb, uyb, wzb
      close(unit=1)
      ! gives the values of the boundary layer 
      ! profile for each node
      do j = 1, jmax
        do i = 1, ptsx
          ux(i,j,1) = uxb(i,j)
          uy(i,j,1) = uyb(i,j)
          wz(i,j,1) = wzb(i,j)
        end do
      end do


      ! reads the derivative and Poisson coefficients
      open(1,file='pre_processing/coefs.bin',form='unformatted')
      read(1) fp_fd_coef
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      close(unit=1)

      open(1,file='coefs.dat',status='unknown')
      write(1,*) 'fp_fd'
      write(1,*) fp_fd_coef
      write(1,*) 'sp_fd'
      write(1,*) sp_fd_coef
      write(1,*) 'cp_fd'
      write(1,*) cp_fd_coef
      write(1,*) 'pp_fd'
      write(1,*) pp_fd_coef
      write(1,*) 'lp_fd'
      write(1,*) lp_fd_coef
      write(1,*) 'fp_sd'
      write(1,*) fp_sd_coef
      write(1,*) 'sp_sd'
      write(1,*) sp_sd_coef
      write(1,*) 'cp_sd'
      write(1,*) cp_sd_coef
      write(1,*) 'pp_sd'
      write(1,*) pp_sd_coef
      write(1,*) 'lp_sd'
      write(1,*) lp_sd_coef
      write(1,*) '1sp_poi'
      write(1,*) sp_poi_coef(:,1)
      write(1,*) '2sp_poi'
      write(1,*) sp_poi_coef(:,2)
      write(1,*) '3sp_poi'
      write(1,*) sp_poi_coef(:,3)
      write(1,*) '4sp_poi'
      write(1,*) sp_poi_coef(:,4)
      write(1,*) '1cp_poi'
      write(1,*) cp_poi_coef(:,1)
      write(1,*) '2cp_poi'
      write(1,*) cp_poi_coef(:,2)
      write(1,*) '3cp_poi'
      write(1,*) cp_poi_coef(:,3)
      write(1,*) '4cp_poi'
      write(1,*) cp_poi_coef(:,4)
      write(1,*) '1pp_poi'
      write(1,*) pp_poi_coef(:,1)
      write(1,*) '2pp_poi'
      write(1,*) pp_poi_coef(:,2)
      write(1,*) '3pp_poi'
      write(1,*) pp_poi_coef(:,3)
      write(1,*) '4pp_poi'
      write(1,*) pp_poi_coef(:,4)
      write(1,*) '1lp_poi'
      write(1,*) lp_poi_coef(:,1)
      write(1,*) '2lp_poi'
      write(1,*) lp_poi_coef(:,2)
      write(1,*) '3lp_poi'
      write(1,*) lp_poi_coef(:,3)
      write(1,*) '4lp_poi'
      write(1,*) lp_poi_coef(:,4)
      write(1,*) 'watw'
      write(1,*) w_at_w_coef
      write(1,*) 'dwydy'
      write(1,*) dwydy_coef
      close(unit=1)

      ! mounts the lhs for the derivative calculation
      call derivs_k

      ! pressure gradient variable calculation
      do i = 1, ptsx
        ueptsx(i) = uxb(i,jmax)
      end do
      call derparxue(duexmdx,ueptsx)
      ! pressure gradient variable calculation

      ! variables used in the buffer domains
      do i = 1, i0
        ep        = dble(i-1) / dble(i0-1)
        bdfc(1,i) = ((6.d0 * ep - 15.d0) * ep + 10.d0) * ep**3
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i0, i3
        bdfc(1,i) = 1.d0
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i3, i4
        ep        = dble(i-i3) / dble(i4-i3)
c       bdfc(1,i) = 1.d0 + ((- 6.d0 * ep + 15.d0) * ep - 10.d0) * ep**3
        bdfc(1,i) = (1.d0 - ep**50)**4 * dexp(-ep**4 / 10.d0)
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i4, imax
        bdfc(1,i) = 0.d0
        bdfc(2,i) = bdfc(1,i)
      end do

      call create_ctes

      ! if the program has stoped, it can be continued by putting
      ! start = 1 in the par.for and recompiling the program
      if (start.eq.1) then 
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(3,file=nome,form='unformatted')
        read(3) t0
        read(3) ux,uy,uz,wx,wy,wz
        close(3)
        t0 = t0 + 1
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_gv(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, wx1, wy1, 
     &                   wz1, dt2, dt6, t0)

      ! initialize the program main variables
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.coef'
      include 'comm.multi'
      include 'mpif.h'
      include 'comm.fs'
      character*23 nome
      integer i, j, k, t0, igv
      real*8 dt2, dt6, a(imax,5), luf(imax,5), bdfc(2,imax), ep,
     &       uxb(ptsx,jmax), uyb(ptsx,jmax), wzb(ptsx,jmax), kc,
     &       afil(ptsx), bfil(ptsx), cfil(ptsx), stf_v, beta_fs_v,
     &       Go(imax), fc, varg(ptsx,jmax), y, m, xad, 
     &       dvargdx(ptsx,jmax), ueptsx(ptsx), x
      complex*16 dv1x(ptsx,jmax,kfour), dv2x(ptsx,jmax,kfour),
     &           dv1y(ptsx,jmax,kfour), dv2y(ptsx,jmax,kfour),
     &           dv1z(ptsx,jmax,kfour), dv2z(ptsx,jmax,kfour),
     &            wx1(ptsx,jmax,kfour),  wy1(ptsx,jmax,kfour),
     &            wz1(ptsx,jmax,kfour),   d2uydx2(ptsx,kfour)
!     complex*16 uxo(ptsx,jmax,11), wxo(ptsx,jmax,11),
!    &           uyo(ptsx,jmax,11), wyo(ptsx,jmax,11),
!    &           uzo(ptsx,jmax,11), wzo(ptsx,jmax,11)
      common/blas/ uxb, uyb, wzb
      common/derw/ d2uydx2
      common/fil/ luf
      common/vc/ varg, dvargdx
      common/bd/ bdfc
      common/filt/ afil, bfil, cfil

      t0  = 1
      dt2 = 0.5d0 * dt
      dt6 = dt / 6.d0

      ! mounts the lu matrix for the filter
      call lhsf(a)
      call ludecomp(a, luf)
      call lhs_tridf(afil,bfil,cfil)

      ! x function of the disturbance strip for TS disturbances
      call var_ts
      call var_gv

      ! all the variables are set to zero
      do k = 1, kfour
        do i = 1, ptsx
          d2uydx2(i,k) = dcmplx(0.d0,0.d0)
          do j = 1, jmax
            ux(i,j,k)   = dcmplx(0.d0,0.d0)
            uy(i,j,k)   = dcmplx(0.d0,0.d0)
            uz(i,j,k)   = dcmplx(0.d0,0.d0)
            wx(i,j,k)   = dcmplx(0.d0,0.d0)
            wy(i,j,k)   = dcmplx(0.d0,0.d0)
            wz(i,j,k)   = dcmplx(0.d0,0.d0)
            wx1(i,j,k)  = dcmplx(0.d0,0.d0)
            wy1(i,j,k)  = dcmplx(0.d0,0.d0)
            wz1(i,j,k)  = dcmplx(0.d0,0.d0)
            dv1x(i,j,k) = dcmplx(0.d0,0.d0)
            dv2x(i,j,k) = dcmplx(0.d0,0.d0)
            dv1y(i,j,k) = dcmplx(0.d0,0.d0)
            dv2y(i,j,k) = dcmplx(0.d0,0.d0)
            dv1z(i,j,k) = dcmplx(0.d0,0.d0)
            dv2z(i,j,k) = dcmplx(0.d0,0.d0)
          end do
        end do
      end do

      ! reads the boundary layer profile
      write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
      open (1, file = nome,form = 'unformatted')
      read(1) uxb, uyb, wzb
      close(unit=1)

      open(1,file='pre_processing/coefs.bin',form='unformatted')
      read(1) fp_fd_coef
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      close(unit=1)

      open(1,file='coefs.dat',status='unknown')
      write(1,*) 'fp_fd'
      write(1,*) fp_fd_coef
      write(1,*) 'sp_fd'
      write(1,*) sp_fd_coef
      write(1,*) 'cp_fd'
      write(1,*) cp_fd_coef
      write(1,*) 'pp_fd'
      write(1,*) pp_fd_coef
      write(1,*) 'lp_fd'
      write(1,*) lp_fd_coef
      write(1,*) 'fp_sd'
      write(1,*) fp_sd_coef
      write(1,*) 'sp_sd'
      write(1,*) sp_sd_coef
      write(1,*) 'cp_sd'
      write(1,*) cp_sd_coef
      write(1,*) 'pp_sd'
      write(1,*) pp_sd_coef
      write(1,*) 'lp_sd'
      write(1,*) lp_sd_coef
      write(1,*) '1sp_poi'
      write(1,*) sp_poi_coef(:,1)
      write(1,*) '2sp_poi'
      write(1,*) sp_poi_coef(:,2)
      write(1,*) '3sp_poi'
      write(1,*) sp_poi_coef(:,3)
      write(1,*) '4sp_poi'
      write(1,*) sp_poi_coef(:,4)
      write(1,*) '1cp_poi'
      write(1,*) cp_poi_coef(:,1)
      write(1,*) '2cp_poi'
      write(1,*) cp_poi_coef(:,2)
      write(1,*) '3cp_poi'
      write(1,*) cp_poi_coef(:,3)
      write(1,*) '4cp_poi'
      write(1,*) cp_poi_coef(:,4)
      write(1,*) '1pp_poi'
      write(1,*) pp_poi_coef(:,1)
      write(1,*) '2pp_poi'
      write(1,*) pp_poi_coef(:,2)
      write(1,*) '3pp_poi'
      write(1,*) pp_poi_coef(:,3)
      write(1,*) '4pp_poi'
      write(1,*) pp_poi_coef(:,4)
      write(1,*) '1lp_poi'
      write(1,*) lp_poi_coef(:,1)
      write(1,*) '2lp_poi'
      write(1,*) lp_poi_coef(:,2)
      write(1,*) '3lp_poi'
      write(1,*) lp_poi_coef(:,3)
      write(1,*) '4lp_poi'
      write(1,*) lp_poi_coef(:,4)
      write(1,*) 'watw'
      write(1,*) w_at_w_coef
      write(1,*) 'dwydy'
      write(1,*) dwydy_coef
      close(unit=1)

      ! mounts the lhs for the derivative calculation
      call derivs_k
      
!     pressure gradient variable calculation
      do i = 1, ptsx
        ueptsx(i) = uxb(i,jmax)
      end do
      call derparxue(duexmdx,ueptsx)
!     pressure gradient variable calculation

      ! definition of varg to be used in the program (curvature term)
      do i = 1, imax
c       Go(i) = Re**(0.25d0) * dsqrt(0.1d0/3.2d0)
c       Go(i) = 2.d0
        Go(i) = 0.d0
      end do
      do i = 1, ptsx
        x = 1.d0 + dble(i-1+shift)*dx
        kc = Go(i+shift) * Go(i+shift) / dsqrt(Re)
        do j = 1, jmax
          varg(i,j) = kc
          select case(cc)
          case(2) 
           varg(i,j) = varg(i,j) *
     &                 2.d0*dsin((5.d0*x+3.d0)*pi/48.d0)
          case(3) 
           varg(i,j) = varg(i,j) *
     &              (-1.d0)*dtanh(3.d0*(x-8.d0))
          case(4) 
           varg(i,j) = varg(i,j) *
     &               0.5d0*(1.d0-dtanh(3.d0*(x-10.24d0)))
          end select        
        end do
      end do

      ! this is used for variable curvature
      call derparxr(dvargdx,varg) 

      ! variables used in the buffer domains
      do i = 1, i0
        ep        = dble(i-1)/dble(i0-1)
        bdfc(1,i) = ((6.d0*ep-15.d0)*ep+10.d0)*ep**3
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i0, i3
        bdfc(1,i) = 1.d0
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i3, i4
        ep        = dble(i-i3)/dble(i4-i3)
c       bdfc(1,i) = 1.d0+((-6.d0*ep+15.d0)*ep-10.d0)*ep**3
        bdfc(1,i) = (1.d0-ep**50)**4*dexp(-ep**4/10.d0)
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i4, imax
        bdfc(1,i) = 0.d0
        bdfc(2,i) = bdfc(1,i)
      end do
      ! buffer domain type 2
      igv = 2 * (i1 + i2)
      do i = 1, igv
        ep        = dble(i-1)/dble(igv-1)
        bdfc(2,i) = ((6.d0*ep-15.d0)*ep+10.d0)*ep**3
      end do

      call create_ctes

      ! if the program has stoped, it can be continued by putting
      ! start = 1 in the par.for and recompiling the program
c     if (start.eq.1) then 
c       write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
c       write(nome,'(a,i0.2,a)')'data_',my_rank,'o.bin'
c       open(3,file=nome,form='unformatted')
c       read(3) t0
c       read(3) ux,uy,uz,wx,wy,wz
c       read(3) uxo,uyo,uzo,wxo,wyo,wzo
c       close(3)
c       t0 = t0 + 1
c       do k = 1, 11
c         do i = 1, ptsx
c           do j = 1, jmax
c             ux(i,j,k) = uxo(i,j,k)
c             uy(i,j,k) = uyo(i,j,k)
c             uz(i,j,k) = uzo(i,j,k)
c             wx(i,j,k) = wxo(i,j,k)
c             wy(i,j,k) = wyo(i,j,k)
c             wz(i,j,k) = wzo(i,j,k)
c           end do
c         end do
c       end do
c     end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_gvth(dv1x, dv2x, dv1y, dv2y, dv1z, dv2z, dv1t, 
     &                     dv2t, wx1, wy1, wz1, th1, dt2, dt6, t0)

      ! initialize the program main variables
      implicit none
      include 'par.for'
      include 'comm.par'
      include 'comm.var'
      include 'comm.coef'
      include 'comm.multi'
      include 'mpif.h'
      include 'comm.fs'
      character*15 nome
      integer i, j, k, t0, igv
      real*8 dt2, dt6, a(imax,5), luf(imax,5), bdfc(2,imax), ep,
     &       thb(ptsx,jmax), uxb(ptsx,jmax), uyb(ptsx,jmax), 
     &       wzb(ptsx,jmax), kc, afil(ptsx), bfil(ptsx), cfil(ptsx), 
     &       stf_v, beta_fs_v, Go(imax), fc, varg(ptsx,jmax), y, m, 
     &       xad,dvargdx(ptsx,jmax), ueptsx(ptsx)
      complex*16 dv1x(ptsx,jmax,kfour), dv2x(ptsx,jmax,kfour),
     &           dv1y(ptsx,jmax,kfour), dv2y(ptsx,jmax,kfour),
     &           dv1z(ptsx,jmax,kfour), dv2z(ptsx,jmax,kfour),
     &           dv1t(ptsx,jmax,kfour), dv2t(ptsx,jmax,kfour),
     &            wx1(ptsx,jmax,kfour),  wy1(ptsx,jmax,kfour),
     &            wz1(ptsx,jmax,kfour),   d2uydx2(ptsx,kfour),
     &            th1(ptsx,jmax,kfour)
      common/blas/ uxb, uyb, wzb
      common/blast/ thb
      common/derw/ d2uydx2
      common/fil/ luf
      common/vc/ varg, dvargdx
      common/bd/ bdfc
      common/filt/ afil, bfil, cfil

      t0  = 1
      dt2 = 0.5d0 * dt
      dt6 = dt / 6.d0

      ! mounts the lu matrix for the filter
      call lhsf(a)
      call ludecomp(a, luf)
      call lhs_tridf(afil,bfil,cfil)

      ! x function of the disturbance strip for TS disturbances
      call var_ts
      call var_gv

      ! definition of varg to be used in the program (curvature term)
      do i = 1, imax
        Go(i) = Re**(0.25d0) * dsqrt(0.1d0/3.2d0)
      end do
      do i = 1, ptsx
        kc = Go(i+shift) * Go(i+shift) / dsqrt(Re)
        do j = 1, jmax
          if(stf.eq.1.d0) then
           y = dble(j-1) * dy
          else
           y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
          endif        
          varg(i,j) = kc
        end do
      end do
      call derparxr(dvargdx,varg) ! this is used for variable curvature

      ! all the variables are set to zero
      do k = 1, kfour
        do i = 1, ptsx
          d2uydx2(i,k) = dcmplx(0.d0,0.d0)
          do j = 1, jmax
            ux(i,j,k)   = dcmplx(0.d0,0.d0)
            uy(i,j,k)   = dcmplx(0.d0,0.d0)
            uz(i,j,k)   = dcmplx(0.d0,0.d0)
            wx(i,j,k)   = dcmplx(0.d0,0.d0)
            wy(i,j,k)   = dcmplx(0.d0,0.d0)
            wz(i,j,k)   = dcmplx(0.d0,0.d0)
            th(i,j,k)   = dcmplx(0.d0,0.d0)
            wx1(i,j,k)  = dcmplx(0.d0,0.d0)
            wy1(i,j,k)  = dcmplx(0.d0,0.d0)
            wz1(i,j,k)  = dcmplx(0.d0,0.d0)
            th1(i,j,k)  = dcmplx(0.d0,0.d0)
            dv1x(i,j,k) = dcmplx(0.d0,0.d0)
            dv2x(i,j,k) = dcmplx(0.d0,0.d0)
            dv1y(i,j,k) = dcmplx(0.d0,0.d0)
            dv2y(i,j,k) = dcmplx(0.d0,0.d0)
            dv1z(i,j,k) = dcmplx(0.d0,0.d0)
            dv2z(i,j,k) = dcmplx(0.d0,0.d0)
            dv1t(i,j,k) = dcmplx(0.d0,0.d0)
            dv2t(i,j,k) = dcmplx(0.d0,0.d0)
          end do
        end do
      end do

      ! reads the boundary layer profile
      write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
      open (1, file = nome,form = 'unformatted')
      read(1) uxb, uyb, wzb, thb
      close(unit=1)

      open(1,file='pre_processing/coefs.bin',form='unformatted')
      read(1) fp_fd_coef
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      read(1) ! integration in the y direction, used in baseflow2D 
      close(unit=1)

      open(1,file='coefs.dat',status='unknown')
      write(1,*) 'fp_fd'
      write(1,*) fp_fd_coef
      write(1,*) 'sp_fd'
      write(1,*) sp_fd_coef
      write(1,*) 'cp_fd'
      write(1,*) cp_fd_coef
      write(1,*) 'pp_fd'
      write(1,*) pp_fd_coef
      write(1,*) 'lp_fd'
      write(1,*) lp_fd_coef
      write(1,*) 'fp_sd'
      write(1,*) fp_sd_coef
      write(1,*) 'sp_sd'
      write(1,*) sp_sd_coef
      write(1,*) 'cp_sd'
      write(1,*) cp_sd_coef
      write(1,*) 'pp_sd'
      write(1,*) pp_sd_coef
      write(1,*) 'lp_sd'
      write(1,*) lp_sd_coef
      write(1,*) '1sp_poi'
      write(1,*) sp_poi_coef(:,1)
      write(1,*) '2sp_poi'
      write(1,*) sp_poi_coef(:,2)
      write(1,*) '3sp_poi'
      write(1,*) sp_poi_coef(:,3)
      write(1,*) '4sp_poi'
      write(1,*) sp_poi_coef(:,4)
      write(1,*) '1cp_poi'
      write(1,*) cp_poi_coef(:,1)
      write(1,*) '2cp_poi'
      write(1,*) cp_poi_coef(:,2)
      write(1,*) '3cp_poi'
      write(1,*) cp_poi_coef(:,3)
      write(1,*) '4cp_poi'
      write(1,*) cp_poi_coef(:,4)
      write(1,*) '1pp_poi'
      write(1,*) pp_poi_coef(:,1)
      write(1,*) '2pp_poi'
      write(1,*) pp_poi_coef(:,2)
      write(1,*) '3pp_poi'
      write(1,*) pp_poi_coef(:,3)
      write(1,*) '4pp_poi'
      write(1,*) pp_poi_coef(:,4)
      write(1,*) '1lp_poi'
      write(1,*) lp_poi_coef(:,1)
      write(1,*) '2lp_poi'
      write(1,*) lp_poi_coef(:,2)
      write(1,*) '3lp_poi'
      write(1,*) lp_poi_coef(:,3)
      write(1,*) '4lp_poi'
      write(1,*) lp_poi_coef(:,4)
      write(1,*) 'watw'
      write(1,*) w_at_w_coef
      write(1,*) 'dwydy'
      write(1,*) dwydy_coef
      close(unit=1)

      ! mounts the lhs for the derivative calculation
      call derivs_k

!     pressure gradient variable calculation
      do i = 1, ptsx
        ueptsx(i) = uxb(i,jmax)
      end do
      call derparxue(duexmdx,ueptsx)
!     pressure gradient variable calculation

      ! variables used in the buffer domains
      do i = 1, i0
        ep        = dble(i-1)/dble(i0-1)
        bdfc(1,i) = ((6.d0*ep-15.d0)*ep+10.d0)*ep**3
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i0, i3
        bdfc(1,i) = 1.d0
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i3, i4
        ep        = dble(i-i3)/dble(i4-i3)
c       bdfc(1,i) = 1.d0+((-6.d0*ep+15.d0)*ep-10.d0)*ep**3
        bdfc(1,i) = (1.d0-ep**50)**4*dexp(-ep**4/10.d0)
        bdfc(2,i) = bdfc(1,i)
      end do
      do i = i4, imax
        bdfc(1,i) = 0.d0
        bdfc(2,i) = bdfc(1,i)
      end do
      ! buffer domain type 2
      igv = 2 * (i1 + i2)
      do i = 1, igv
        ep        = dble(i-1)/dble(igv-1)
        bdfc(2,i) = ((6.d0*ep-15.d0)*ep+10.d0)*ep**3
      end do

      call create_ctes

      ! if the program has stoped, it can be continued by putting
      ! start = 1 in the par.for and recompiling the program
      if (start.eq.1) then 
        write(nome,'(a,i0.2,a)')'data_',my_rank,'.bin'
        open(3,file=nome,form='unformatted')
        read(3) t0
        read(3) ux,uy,uz,wx,wy,wz,th
        close(3)
        t0 = t0 + 1
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv_fi(dvx, dvy, dvz)

      ! calculate the derivatives for RK method
      implicit none
      include 'par.for'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      complex*16 d2wxdx2(ptsx,jmax,kfour),d2wxdy2(ptsx,jmax,kfour),
     &           d2wydx2(ptsx,jmax,kfour),d2wydy2(ptsx,jmax,kfour),
     &           d2wzdx2(ptsx,jmax,kfour),d2wzdy2(ptsx,jmax,kfour),
     &             dfydx(ptsx,jmax,kfour),  dfzdx(ptsx,jmax,kfour),
     &             dfxdy(ptsx,jmax,kfour),  dfzdy(ptsx,jmax,kfour),
     &              dadx(ptsx,jmax,kfour),   dady(ptsx,jmax,kfour),
     &              dbdx(ptsx,jmax,kfour),   dcdy(ptsx,jmax,kfour),
     &               dvx(ptsx,jmax,kfour),    dvy(ptsx,jmax,kfour),
     &               dvz(ptsx,jmax,kfour),     fy(ptsx,jmax,kfour),
     &                fx(ptsx,jmax,kfour),     fz(ptsx,jmax,kfour),
     &                 a(ptsx,jmax,kfour),      b(ptsx,jmax,kfour),
     &                 c(ptsx,jmax,kfour)
      common/force/ fx, fy, fz

      ! linear or non-linear product calculations
c     call lterms_fi(a, b, c)
      call nlterms_fi(a, b, c)

      ! derivative calculations
      call derparx(dadx, a)
      call derparx(dbdx, b)

      call dery(dady, a)
      call dery(dcdy, c)

      call derparxx(d2wxdx2, wx)
      call derparxx(d2wydx2, wy)
      call derparxx(d2wzdx2, wz)

      call deryy(d2wxdy2, wx)
      call deryy(d2wydy2, wy)
      call deryy(d2wzdy2, wz)

      ! immersed boundary method forcing terms and derivative calculations
!     call cvirt
!     call derparx(dfydx, fy)
!     call derparx(dfzdx, fz)
!     call dery(dfxdy, fx)
!     call dery(dfzdy, fz)

      do k = 1, kfour
        do j = 2, jmax
          do i = 1, ptsx

            dvx(i,j,k) = - dady(i,j,k) + v_kb(k) * b(i,j,k)
     &                   + ( d2wxdx2(i,j,k) + d2wxdy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wx(i,j,k) ) / Re
c    &                   + v_kb(k) * fy(i,j,k) / fac_y - dfzdy(i,j,k)

            dvy(i,j,k) = - v_kb(k) * c(i,j,k) + dadx(i,j,k)
     &                   + ( d2wydx2(i,j,k) + d2wydy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wy(i,j,k) ) / Re
c    &                   + dfzdx(i,j,k) - v_kb(k) * fx(i,j,k)

            dvz(i,j,k) = - dbdx(i,j,k) + dcdy(i,j,k)
     &                   + ( d2wzdx2(i,j,k) + d2wzdy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wz(i,j,k) ) / Re
c    &                   + dfxdy(i,j,k)  -  dfydx(i,j,k) / fac_y

          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv_gv(dvx, dvy, dvz)

      ! calculate the derivatives for RK method
      implicit none
      include 'par.for'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8 kc, varg(ptsx,jmax), dvargdx(ptsx,jmax)
      complex*16 d2wxdx2(ptsx,jmax,kfour),d2wxdy2(ptsx,jmax,kfour),
     &           d2wydx2(ptsx,jmax,kfour),d2wydy2(ptsx,jmax,kfour),
     &           d2wzdx2(ptsx,jmax,kfour),d2wzdy2(ptsx,jmax,kfour),
     &              dadx(ptsx,jmax,kfour),   dady(ptsx,jmax,kfour),
     &              dbdx(ptsx,jmax,kfour),   dcdy(ptsx,jmax,kfour),
     &               dvx(ptsx,jmax,kfour),    dvy(ptsx,jmax,kfour),
     &               dvz(ptsx,jmax,kfour),   dddx(ptsx,jmax,kfour),
     &                 a(ptsx,jmax,kfour),      b(ptsx,jmax,kfour),
     &                 c(ptsx,jmax,kfour),      d(ptsx,jmax,kfour)
      common/vc/ varg, dvargdx

      ! linear or non-linear product calculations
c     call lterms_gv(a, b, c, d)
      call nlterms_gv(a, b, c, d)

      ! derivative calculations
      call derparx(dadx, a)
      call derparx(dbdx, b)
      call derparx(dddx, d)

      call dery(dady, a)
      call dery(dcdy, c)

      call derparxx(d2wxdx2, wx)
      call derparxx(d2wydx2, wy)
      call derparxx(d2wzdx2, wz)

      call deryy(d2wxdy2, wx)
      call deryy(d2wydy2, wy)
      call deryy(d2wzdy2, wz)

      do k = 1, kfour
        do j = 2, jmax
          do i = 1, ptsx

            dvx(i,j,k) = - dady(i,j,k) + v_kb(k) * b(i,j,k)
     &                   - v_kb(k) * d(i,j,k) * varg(i,j)
     &                   + ( d2wxdx2(i,j,k) + d2wxdy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wx(i,j,k) ) / Re

            dvy(i,j,k) = - v_kb(k) * c(i,j,k) + dadx(i,j,k)
     &                   + ( d2wydx2(i,j,k) + d2wydy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wy(i,j,k) ) / Re

            dvz(i,j,k) = - dbdx(i,j,k) + dcdy(i,j,k)
     &                   + varg(i,j) * dddx(i,j,k)
     &                   + d(i,j,k) * dvargdx(i,j)
     &                   + ( d2wzdx2(i,j,k) + d2wzdy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wz(i,j,k) ) / Re

          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv_gvth(dvx, dvy, dvz, dvt)

      ! calculate the derivatives for RK method
      implicit none
      include 'par.for'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      real*8 kc, varg(ptsx,jmax), dvargdx(ptsx,jmax)
      complex*16 d2wxdx2(ptsx,jmax,kfour),d2wxdy2(ptsx,jmax,kfour),
     &           d2wydx2(ptsx,jmax,kfour),d2wydy2(ptsx,jmax,kfour),
     &           d2wzdx2(ptsx,jmax,kfour),d2wzdy2(ptsx,jmax,kfour),
     &           d2thdx2(ptsx,jmax,kfour),d2thdy2(ptsx,jmax,kfour),
     &            duthdx(ptsx,jmax,kfour), dvthdy(ptsx,jmax,kfour),
     &              dadx(ptsx,jmax,kfour),   dady(ptsx,jmax,kfour),
     &              dbdx(ptsx,jmax,kfour),   dcdy(ptsx,jmax,kfour),
     &              dddx(ptsx,jmax,kfour),    wth(ptsx,jmax,kfour),
     &               dvx(ptsx,jmax,kfour),    dvy(ptsx,jmax,kfour),
     &               dvz(ptsx,jmax,kfour),    dvt(ptsx,jmax,kfour),
     &               uth(ptsx,jmax,kfour),    vth(ptsx,jmax,kfour),
     &                 a(ptsx,jmax,kfour),      b(ptsx,jmax,kfour),
     &                 c(ptsx,jmax,kfour),      d(ptsx,jmax,kfour)
      common/vc/ varg, dvargdx

      ! linear or non-linear product calculations
c     call lterms_th(a,b,c,d,uth,vth,wth)
      call nlterms_th(a,b,c,d,uth,vth,wth)

      ! derivative calculations
      call derparx(dadx,a)
      call derparx(dbdx,b)
      call derparx(dddx,d)
      call derparx(duthdx,uth)

      call dery(dady,a)
      call dery(dcdy,c)
      call dery(dvthdy,vth)

      call derparxx(d2wxdx2,wx)
      call derparxx(d2wydx2,wy)
      call derparxx(d2wzdx2,wz)
      call derparxx(d2thdx2,th)

      call deryy(d2wxdy2,wx)
      call deryy(d2wydy2,wy)
      call deryy(d2wzdy2,wz)
      call deryy(d2thdy2,th)

      do k = 1, kfour
        do j = 2, jmax
          do i = 1, ptsx

            dvx(i,j,k) = - dady(i,j,k) + v_kb(k) * b(i,j,k)
     &                   - v_kb(k) * d(i,j,k) * varg(i,j)
     &                   + ( d2wxdx2(i,j,k) + d2wxdy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wx(i,j,k) ) / Re

            dvy(i,j,k) = - v_kb(k) * c(i,j,k) + dadx(i,j,k)
     &                   + ( d2wydx2(i,j,k) + d2wydy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wy(i,j,k) ) / Re

            dvz(i,j,k) = - dbdx(i,j,k) + dcdy(i,j,k)
     &                   + varg(i,j) * dddx(i,j,k)
c    &                   + d(i,j,k) * dvargdx(i,j)
     &                   + ( d2wzdx2(i,j,k) + d2wzdy2(i,j,k) * fac_y
     &                   +  v_k2b2(k) * wz(i,j,k) ) / Re

c VERIFICAR COMO FUNCIONA O FAC_Y PARA EQUAÇÃO DE TRANSPORTE DE THETA
            dvt(i,j,k) = - duthdx(i,j,k) - dvthdy(i,j,k)
     &                   - v_kb(k) * wth(i,j,k)
     &                   + ( d2thdx2(i,j,k) + d2thdy2(i,j,k)
     &                   +  v_k2b2(k) * th(i,j,k) )/(Re*Pr)

          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drv_les(dvx, dvy, dvz)

      ! calculate the derivatives for RK method
      implicit none
      include 'par.for'
      include 'comm.var'
      include 'comm.fourier'
      integer i, j, k
      complex*16 d2wxdx2(ptsx,jmax,kfour),d2wxdy2(ptsx,jmax,kfour),
     &           d2wydx2(ptsx,jmax,kfour),d2wydy2(ptsx,jmax,kfour),
     &           d2wzdx2(ptsx,jmax,kfour),d2wzdy2(ptsx,jmax,kfour),
     &             dfydx(ptsx,jmax,kfour),  dfzdx(ptsx,jmax,kfour),
     &             dfxdy(ptsx,jmax,kfour),  dfzdy(ptsx,jmax,kfour),
     &              dadx(ptsx,jmax,kfour),   dady(ptsx,jmax,kfour),
     &              dbdx(ptsx,jmax,kfour),   dcdy(ptsx,jmax,kfour),
     &               dvx(ptsx,jmax,kfour),    dvy(ptsx,jmax,kfour),
     &               dvz(ptsx,jmax,kfour),     fy(ptsx,jmax,kfour),
     &                fx(ptsx,jmax,kfour),     fz(ptsx,jmax,kfour),
     &                 a(ptsx,jmax,kfour),      b(ptsx,jmax,kfour),
     &                 c(ptsx,jmax,kfour)!,
c    &             F_nux(ptsx,jmax,kfour),  F_nuy(ptsx,jmax,kfour),
c    &             F_nuz(ptsx,jmax,kfour)
      common/force/ fx, fy, fz

      ! linear or non-linear product calculations
c     call lterms_fi(a, b, c)
      call nlterms_fi(a, b, c)

      ! derivative calculations
      call derparx(dadx, a)
      call derparx(dbdx, b)

      call dery(dady, a)
      call dery(dcdy, c)

      call derparxx(d2wxdx2, wx)
      call derparxx(d2wydx2, wy)
      call derparxx(d2wzdx2, wz)

      call deryy(d2wxdy2, wx)
      call deryy(d2wydy2, wy)
      call deryy(d2wzdy2, wz)

c      call les_terms(F_nux, F_nuy, F_nuz)
       call les_terms(fx, fy, fz) ! only LES terms without immersed boundary

      ! immersed boundary method forcing terms and derivative calculations
c     call cvirt
c     do k = 1, kfour
c       do j = 1, jmax
c         do i = 1, ptsx
c           fx(i,j,k) = fx(i,j,k) + F_nux(i,j,k)
c           fy(i,j,k) = fy(i,j,k) + F_nuy(i,j,k)
c           fz(i,j,k) = fz(i,j,k) + F_nuz(i,j,k)
c         end do
c       end do
c     end do

      call derparx(dfydx, fy)
      call derparx(dfzdx, fz)
      call dery(dfxdy, fx)
      call dery(dfzdy, fz)

      do k = 1, kfour
        do j = 2, jmax
          do i = 1, ptsx

            dvx(i,j,k) = - dady(i,j,k) + v_kb(k) * b(i,j,k)
     &                   + ( d2wxdx2(i,j,k) + d2wxdy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wx(i,j,k) ) / Re
     &                   + v_kb(k) * fy(i,j,k) / fac_y - dfzdy(i,j,k)

            dvy(i,j,k) = - v_kb(k) * c(i,j,k) + dadx(i,j,k)
     &                   + ( d2wydx2(i,j,k) + d2wydy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wy(i,j,k) ) / Re
     &                   + dfzdx(i,j,k) - v_kb(k) * fx(i,j,k)

            dvz(i,j,k) = - dbdx(i,j,k) + dcdy(i,j,k)
     &                   + ( d2wzdx2(i,j,k) + d2wzdy2(i,j,k) * fac_y
     &                   + v_k2b2(k) * wz(i,j,k) ) / Re
     &                   + dfxdy(i,j,k)  -  dfydx(i,j,k) / fac_y

          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gv_pert(t)

      ! introduces GV perturbations in the flow between i1 and i2
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, k, t
      real*8 A, ep, fcx2(i2), fcllx2(i2)
      complex*16 d2uydx2(ptsx,kfour)
      common/equa2/ fcx2,fcllx2
      common/derw/ d2uydx2

c     A = 1d-4
      A = 7.5d-3 * dsqrt(fac_y)
c     A = 1.85d-4
      
c     if (t.eq.1) call wdata2(a)
      
      ! this is to initialize the amplitude of
      ! disturbance smoothly in 1 period (lamb=stpp*dt)
      if (t.lt.stpp) then
        ep = dble(t)/dble(stpp)
        A  = A*((6.d0*ep-15.d0)*ep+10.d0)*ep**3
      end if
      
      do i = i1 + 1, i2 - 1
       do k = 2, 2
        uy(i,1,k)    = dcmplx( A*fcx2(i), 0.d0 )
        d2uydx2(i,k) = dcmplx( A*fcllx2(i), 0.d0 )
       enddo 
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ts2d_pert(t, tempt)

      ! introduces perturbations in the flow between i1 and i2
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, t
      real*8 A, H, ept, tempt, fct
      real*8 fcx(i2), fcllx(i2)
      complex*16 d2uydx2(ptsx,kfour)
      common/derw/ d2uydx2
      common/equa/ fcx, fcllx

      ! here one must choose which amplitude and frequency 
      ! disturbances are going to be applied
      A = 1.d-5 * dsqrt(fac_y) ! 0,75 %
c     A = 1.d-6 * dsqrt(fac_y) ! 0,75 %
      H = 1.d0

      fct = A * dsin(H * ( dble(t) + tempt ) * omega * dt)

      ! initialize the amplitude of disturbance smoothly
      ! in 1 period (lamb = stpp * dt)
      if (t.lt.stpp) then
        ept = dble(t) / dble(stpp)
        fct = fct * ( ( 6.d0 * ept - 15.d0 ) * ept + 10.d0 ) * ept**3
      end if

      do i = i1 + 1, i2 - 1
        uy(i,1,1)    = dcmplx(fcx(i)   * fct,0.d0)
        d2uydx2(i,1) = dcmplx(fcllx(i) * fct,0.d0)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ts3d_pert(t, tempt)

      ! introduces perturbations in the flow between i1 and i2
      implicit none
      include 'par.for'
      include 'comm.var'
      integer i, t
      real*8 A, H, ept, tempt, fct
      real*8 fcx(i2), fcllx(i2)
      complex*16 d2uydx2(ptsx,kfour)
      common/derw/ d2uydx2
      common/equa/ fcx, fcllx

      ! here you must choose which amplitude and frequency
      ! disturbances are you going to be applied
      A = 1.d-5 * dsqrt(fac_y)
      H = 1.d0

c     if (t.eq.1) call wdata(A,H)

      fct = A * dsin(H * ( dble(t) + tempt ) * omega * dt)

      ! initialize the amplitude of disturbance smoothly
      ! in 1 period (lamb = stpp * dt)
      if (t.lt.stpp) then
        ept = dble(t) / dble(stpp)
        fct = fct * ( ( 6.d0 * ept - 15.d0 ) * ept + 10.d0 ) * ept**3
      end if

      do i = i1 + 1, i2 - 1
        uy(i,1,2)    = dcmplx(fcx(i)   * fct, 0.d0)
        d2uydx2(i,2) = dcmplx(fcllx(i) * fct, 0.d0)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine var_ts

      implicit none
      include 'par.for'
      integer i
      real*8 x, x1, fcx(i2), fcllx(i2)
      real*8 a, b, c, d, e, f, imed, div, fcmax
      common/equa/ fcx, fcllx

      do i = 1, i2
        fcx(i)   = 0.d0
        fcllx(i) = 0.d0
      end do

      imed = dble(i1 + i2) / 2.d0
      x1   = x0 + dble(i1 - 1) * dx
      div  = 1.d0 / ( dble(imed - i1) * dx )

      a = -  19683.d0 / 4096.d0 * div**8
      b =   177147.d0 / 4096.d0 * div**7
      c = -  19683.d0 /  128.d0 * div**6
      d =   137781.d0 /  512.d0 * div**5
      e = -  59049.d0 /  256.d0 * div**4
      f =    19683.d0 /  256.d0 * div**3

      fcmax = 0.d0
      do i = i1, (i1+i2)/2
        x      = dble(i - i1) * dx
        fcx(i) = ( ( ( ( ( a * x + b ) * x + c ) * x + d ) * x + e )
     &         * x + f ) * x**4
        fcmax  = max(fcmax, dabs(fcx(i)))
      end do

      do i = i1, (i1+i2)/2
        x      = dble(i - i1) * dx
        fcx(i) = ( ( ( ( ( a * x + b ) * x + c ) * x + d ) * x + e )
     &         * x + f ) * x**4 / fcmax
        if (dabs(fcx(i)).lt.1d-15) fcx(i) = 0.d0
        fcx(i2+i1-i) = - fcx(i)
        fcllx(i)     = ( ( ( ( ( 72.d0 * a   * x + 56.d0 * b ) * x
     &                         + 42.d0 * c ) * x + 30.d0 * d ) * x
     &                         + 20.d0 * e ) * x + 12.d0 * f ) * x
     &               * x / fcmax
        if (dabs(fcllx(i)).lt.1d-15) fcllx(i) = 0.d0
        fcllx(i2+i1-i) = - fcllx(i)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine var_gv

      implicit none
      include 'par.for'
      integer i
      real*8 x, fcx2(i2), fcllx2(i2)
      common/equa2/ fcx2,fcllx2

      do i = 1, i2
        fcx2(i)   = 0.d0
        fcllx2(i) = 0.d0
      end do

      do i = i1, i2
        x         = pi*dble(i-i1)/dble(i2-i1)
        fcx2(i)   = (dsin(x))**3
        fcllx2(i) = 6.d0*(dcos(x))**2*dsin(x)-3.d0*(dsin(x))**3
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loop(erro)

      ! actualize the variables for each RK-step
      implicit none
      include 'par.for'
      real*8 erro
      complex*16 lapv(ptsx,kfour)
      complex*16 duydy(ptsx,jmax,kfour), dwzdx(ptsx,jmax,kfour)

      call bzone
c     if (erro.lt.1.d-5) call filter
      if (erro.lt.1.d-5) call filter_trid
c     call filter_trid
      call outuy(dwzdx)
      call poi_uy(dwzdx, erro)
      call poi_ux(duydy)
      call poi_uz(duydy)
      call wx_wall(lapv)
      call wz_wall(lapv)
      call bzonew

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wdata(A, H)

      implicit none
      include 'par.for'
      real*8 A, H

      open (1, file = 'datats.dat',status = 'unknown')
      write(1,*) 'Data used in the simulation'
      write(1,*)
      write(1,*) 'Imax =',imax
      write(1,*) 'ptsx =',ptsx
      write(1,*) 'dx =',dx
      write(1,*) 'i1 =',i1
      write(1,*) 'i2 =',i2
      write(1,*) 'i3 =',i3
      write(1,*) 'i4 =',i4
      write(1,*)
      write(1,*) 'Jmax =',jmax
      write(1,*) 'dy =',dy
      write(1,*)
      write(1,*) 'Kmax =',kfour
      write(1,*) 'Kphys =',kphys
      write(1,*)
      write(1,*) 'tt =',tt
      write(1,*) 'dt =',dt
      write(1,*) 'steps/period =',stpp
      write(1,*) 'CFL =',dt/dx
      write(1,*)
      write(1,*) 'Beta(spanwise wavenumber) =',beta
      write(1,*) 'Alpha(streamwise wavenumber) =',alpha
      write(1,*) 'Omega(frequency) =',omega
      write(1,*)
      write(1,*) 'A =',A
      write(1,*) 'H =',H
      close (unit=1)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine create_ctes

      implicit none
      include 'par.for'
      include 'comm.multi'   
      include 'comm.fourier'   
 
      integer lvl, j, k 
      real*8 aux
 
      ! Multigrid spatial calculations
 
      ! dy0 at each multigrid level
      do lvl = 1 , msh
       if(stf.ne.1.d0) then
        v_dy0(lvl) = dy * ( ( stf**(2**(lvl-1))-1.d0) / (stf-1.d0) )
       else
        v_dy0(lvl) = dy * 2.d0**(lvl-1)
       endif 
      enddo

      do lvl = 1 , msh
       v_stf(lvl) = stf ** ( 2**(lvl-1) )
       v_dx2(lvl) = (dx * dble(2**(lvl-1)))**2
      enddo

      ! dy at each space
       if (stf.ne.1.d0) then
        do j = 1 , jmax - 1
         v_qdy(j)  = 1.d0 / (v_dy0(1) * v_stf(1)**(j-1))
        enddo
       else
        do j = 1 , jmax - 1
         v_qdy(j)  = 1.d0 / (v_dy0(1))
        enddo
      endif

      do lvl = 1 , msh
       if (stf.ne.1.d0) then
        do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
         aux           = (v_dy0(lvl) * v_stf(lvl)**(j-1))
         v_qdy2(j,lvl) = 1.d0 / (aux**2)
        enddo
       else
        do j = 1 , ( jmax - 1 ) / 2**(lvl-1)
         aux           = (v_dy0(lvl))
         v_qdy2(j,lvl) = 1.d0 / (aux**2)
        enddo
       endif
      enddo

      v_ptsx(1) = ptsx
      v_ptsy(1) = jmax
      do lvl = 2 , msh
       v_ptsx(lvl) = (v_ptsx(lvl-1) + 1) / 2
       v_ptsy(lvl) = (v_ptsy(lvl-1) + 1) / 2
      enddo

      do k = 1 , kfour
       v_k2b2(k) = - dble(k-1) * dble(k-1) * beta * beta
       v_kb(k)   = - im * dble(k-1) * beta
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c                end of main subroutines                         c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
