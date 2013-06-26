cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                          main subroutines                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program baseflow_16082012

      implicit none
      include '../par.for'
      include '../comm.par'
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

      ! calculate the i_shift from one computer to another with the domain
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
      include '../par.for'
      include '../comm.par'
      logical run
      integer nodesx, chute1, chute2
      real*8 n_npr

      ! calculates the value of points in the x direction for each node
      nodesx = ( imax + (inter + 1) * numproc ) / (numproc + 1)
      n_npr  = dble( imax + (inter + 1) * numproc ) / dble(numproc + 1)

      ! Stop the program if the number of total points is not exactly
      ! divided by nodes and if the number of points in the x
      ! direction can not be used in multigrid program and if the
      ! number of points in the y direction can not be used in
      ! multigrid program

      if (ptsx.ne.nodesx) then
        write(*,*) 'ALTERAR O VALOR DE PTSX NO ARQUIVO par.for PARA:'
        write(*,*) nodesx
        run = .false.
      end if

      if ((nodesx.ne.n_npr).or.(mod(nodesx-1,2**(msh-1)).ne.0)) then
        chute1 = (   (nodesx - 1) / 2**(msh - 1) )       * 2**(msh - 1)
        chute2 = ( ( (nodesx - 1) / 2**(msh - 1) ) + 1 ) * 2**(msh - 1)
        write(*,*)' Number of points in the x direction can
     &              not be used in multigrid solver'
        write(*,*)' The number of points in the x direction
     &              should be:'
        write(*,*)  chute1 * (numproc + 1) - numproc * inter + 1,'   or'
        write(*,*)  chute2 * (numproc + 1) - numproc * inter + 1
        run = .false.
      end if

      if (nodesx .lt. 2*inter) then
        write(*,*) ' Number of points in the x direction can not be
     &               used in multigrid solver'
        write(*,*) ' The number of points in the x direction should
     &               be:'
        write(*,*) ((nodesx * 2) - 1) * (numproc + 1) - numproc
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
        write(*,*) ' Number of points in the y direction can not be
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
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include 'mpif.h'
      integer i, j, t, t0, i_ini
      real*8 dv1(ptsx,jmax),  wz1(ptsx,jmax), 
     &       dvt1(ptsx,jmax), th1(ptsx,jmax)

      i_ini = 1
      if (my_rank.eq.0) i_ini = 2

      if (my_form.eq.2) then 
        call init_theta(dv1, wz1, dvt1, th1)
        do t = 1, tt_base
  
          call drv_theta(dv1, dvt1)
          do j = 2, jmax
            do i = i_ini, ptsx
              wz1(i,j) = wz(i,j)
              th1(i,j) = th(i,j)
              wz(i,j)  = wz1(i,j) + dv1(i,j)  * dt_base
              th(i,j)  = th1(i,j) + dvt1(i,j) * dt_base
            end do
          end do
          call loop(1d-4)
  
          call drv_theta(dv1, dvt1)
          do j = 2, jmax
            do i = i_ini, ptsx
              wz(i,j)  = wz1(i,j) + dv1(i,j)  * dt_base
              th(i,j)  = th1(i,j) + dvt1(i,j) * dt_base
            end do
          end do
          call loop(1d-5)
  
          write(*,*)my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax)
          if (mod(t,500).eq.0) call escreve_theta
  
        end do
        call escreve_theta
       else ! my_form = 0 and my_form = 1

        call init(dv1, wz1)
        do t = 1, tt_base
  
          call drv(dv1)
          do j = 2, jmax
            do i = i_ini, ptsx
              wz1(i,j) = wz(i,j)
              wz(i,j)  = wz1(i,j) + dv1(i,j) * dt_base
            end do
          end do
          call loop(1d-4)
  
          call drv(dv1)
          do j = 2, jmax
            do i = i_ini, ptsx
              wz(i,j)  = wz1(i,j) + dv1(i,j) * dt_base
            end do
          end do
          call loop(1d-5)
  
          write(*,*)my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax)
          if (mod(t,500).eq.0) call escreve
  
        end do
        call escreve
      end if


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init(dv1, wz1)

      ! initialize the program main variables
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include '../comm.multi'
      include '../comm.fs'
      include 'mpif.h'
      integer i, j
      real*8 m, stf_verif, dv1(ptsx,jmax), wz1(ptsx,jmax),
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax), xad,
     &       dwzdx(ptsx,jmax), ue, afil(ptsx), bfil(ptsx), cfil(ptsx),
     &       ueptsx(ptsx)
      common/dwdx/ dwzdx
      common/filt/ afil, bfil, cfil

      call lhs_tridf(afil,bfil,cfil)

      ! all the variables are set to zero
      do i = 1, ptsx
        do j = 1, jmax
          ux(i,j)  = 0.d0
          uy(i,j)  = 0.d0
          wz(i,j)  = 0.d0
          wz1(i,j) = 0.d0
          dv1(i,j) = 0.d0
        end do
      end do

      ! reads the boundary layer profile
      open(1,file='../pre_processing/base_fs.bin',form='unformatted')
      read(1) uxbt, uybt, wzbt
      close(unit=1)

      ! gives the values of the boundary layer profile for each node
      do j = 1, jmax
        do i = 1, ptsx
          ux(i,j) = uxbt(i+shift,j)
          uy(i,j) = uybt(i+shift,j)
          wz(i,j) = wzbt(i+shift,j)
        end do
      end do

      ! reads the derivative and Poisson coefficients
      open(1,file='../pre_processing/coefs.bin',form='unformatted')
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
      read(1) sp_integ_coef
      read(1) cp_integ_coef
      read(1) pp_integ_coef
      read(1) lp_integ_coef
      close(unit=1)

      ! mounts the lhs for the derivative calculation
      call derivs_k

      ! defines boundary layer parameters
      ! reads beta_fs from a file
      open(1,file='../beta_fs.dist',form='formatted')
      read(1,*) beta_fs
      close(unit=1)

      if (my_rank.eq.0) then
        ue = ux(1,jmax)
      end if
      call MPI_BCAST(ue, 1, mpi_double_precision, 0, mpi_comm_world,
     &               ierr)
      do i = 1, ptsx
        xad        = dble(i+shift-1)*dx + x0
        m          = beta_fs(i+shift) / (2.d0 - beta_fs(i+shift))
        ux(i,jmax) = ue * xad**m
        ueptsx(i)  = ux(i,jmax)
      end do
      call derparxue(duexmdx,ueptsx)
      ! defines boundary layer parameters

      call create_ctes

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init_theta(dv1, wz1, dvt1, th1)

      ! initialize the program main variables
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include '../comm.multi'
      include '../comm.fs'
      include 'mpif.h'
      integer i, j
      real*8 m, stf_verif, dv1(ptsx,jmax), wz1(ptsx,jmax),
     &       dvt1(ptsx,jmax), th1(ptsx,jmax), thbt(imax,jmax),
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax), xad,
     &       dwzdx(ptsx,jmax), ue, afil(ptsx), bfil(ptsx), cfil(ptsx)
      common/dwdx/ dwzdx
      common/filt/ afil, bfil, cfil

      call lhs_tridf(afil,bfil,cfil)

      ! all the variables are set to zero
      do i = 1, ptsx
        do j = 1, jmax
          ux(i,j)   = 0.d0
          uy(i,j)   = 0.d0
          wz(i,j)   = 0.d0
          th(i,j)   = 0.d0
          wz1(i,j)  = 0.d0
          dv1(i,j)  = 0.d0
          dvt1(i,j) = 0.d0
          th1(i,j)  = 0.d0
        end do
      end do

      ! reads the boundary layer profile
      open(1,file='../pre_processing/base_fs.bin',form='unformatted')
      read(1) uxbt, uybt, wzbt, thbt
      close(unit=1)

      ! gives the values of the boundary layer profile for each node
      do j = 1, jmax
        do i = 1, ptsx
          ux(i,j) = uxbt(i+shift,j)
          uy(i,j) = uybt(i+shift,j)
          wz(i,j) = wzbt(i+shift,j)
          th(i,j) = thbt(i+shift,j)
        end do
      end do

      ! reads beta_fs from a file
      open(1,file='../beta_fs.dist',form='formatted')
      read(1,*) beta_fs
      close(unit=1)

      ! defines boundary layer parameters
      if (my_rank.eq.0) then
        ue = ux(1,jmax)
      end if
      call MPI_BCAST(ue, 1, mpi_double_precision, 0, mpi_comm_world,
     &               ierr)
      do i = 1, ptsx
        xad        = dble(i+shift-1)*dx + x0
        m          = beta_fs(i+shift) / (2.d0 - beta_fs(i+shift))
        ux(i,jmax) = ue * xad**m
        duexmdx(i) = ue * m * xad**(m - 1.d0)
      end do

      ! reads the derivative and Poisson coefficients
      open(1,file='../pre_processing/coefs.bin',form='unformatted')
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
      read(1) sp_integ_coef
      read(1) cp_integ_coef
      read(1) pp_integ_coef
      read(1) lp_integ_coef
      close(unit=1)

      ! mounts the lhs for the derivative calculation
      call derivs_k

      call create_ctes

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine drv(dv)

      ! calculate the derivatives for RK method
      implicit none
      include '../par.for'
      include 'comm.var'
      integer i, j
      real*8 d2wzdx2(ptsx,jmax), d2wzdy2(ptsx,jmax),
     &           uwz(ptsx,jmax),     vwz(ptsx,jmax),
     &        duwzdx(ptsx,jmax),  dvwzdy(ptsx,jmax),
     &            dv(ptsx,jmax)

      do j = 1, jmax
        do i = 1, ptsx
          uwz(i,j) = ux(i,j) * wz(i,j)
          vwz(i,j) = uy(i,j) * wz(i,j)
        end do
      end do
      call derparx(duwzdx,uwz)
      call dery(dvwzdy, vwz)

      call derparxx(d2wzdx2, wz)
      call deryy(d2wzdy2, wz)

      do j = 2, jmax
        do i = 1, ptsx
          dv(i,j) = - duwzdx(i,j) - dvwzdy(i,j)
     &              + ( d2wzdx2(i,j) + d2wzdy2(i,j) * fac_y ) / Re
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine drv_theta(dv, dvt)

      ! calculate the derivatives for RK method
      implicit none
      include '../par.for'
      include 'comm.var'
      integer i, j
      real*8 d2wzdx2(ptsx,jmax), d2wzdy2(ptsx,jmax),
     &           uwz(ptsx,jmax),     vwz(ptsx,jmax),
     &        duwzdx(ptsx,jmax),  dvwzdy(ptsx,jmax),
     &            dv(ptsx,jmax),     uth(ptsx,jmax),
     &           vth(ptsx,jmax), d2thdx2(ptsx,jmax),
     &       d2thdy2(ptsx,jmax),     dvt(ptsx,jmax),
     &        duthdx(ptsx,jmax),  dvthdy(ptsx,jmax)

      do j = 1, jmax
        do i = 1, ptsx
          uwz(i,j) = ux(i,j) * wz(i,j)
          vwz(i,j) = uy(i,j) * wz(i,j)
          uth(i,j) = ux(i,j) * th(i,j)
          vth(i,j) = uy(i,j) * th(i,j)
        end do
      end do
      call derparx(duwzdx,uwz)
      call derparx(duthdx,uth)

      call dery(dvwzdy, vwz)
      call dery(dvthdy, vth)

      call derparxx(d2wzdx2, wz)
      call derparxx(d2thdx2, th)

      call deryy(d2wzdy2, wz)
      call deryy(d2thdy2, th)

      do j = 2, jmax
        do i = 1, ptsx
          dv(i,j) = - duwzdx(i,j) - dvwzdy(i,j)
     &              + ( d2wzdx2(i,j) + d2wzdy2(i,j) ) / Re

          dvt(i,j) = - duthdx(i,j) - dvthdy(i,j)
     &               + ( d2thdx2(i,j) + d2thdy2(i,j))/(Re*Pr)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine loop(erro)

      ! actualize the variables for each RK-step
      implicit none
      include '../par.for'
      real*8 erro

c     if (erro.lt.1d-4) call filter_trid
      call outuy
      call poi_uy(erro)
      call poi_ux
      call wz_wall

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine create_ctes

      implicit none
      include '../par.for'
      include '../comm.multi'
      integer lvl, j
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
 
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve

      ! write the results in binary form
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      character*20 nm
 
      write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
      open(1,file=nm,form='unformatted')
      write(1) ux, uy, wz
      close (unit=1)

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve_theta

      ! write the results in binary form
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      character*20 nm
 
      write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
      open(1,file=nm,form='unformatted')
      write(1) ux, uy, wz, th
      close (unit=1)

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c                      end of main subroutines                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
