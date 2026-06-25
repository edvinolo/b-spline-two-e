module time_propagation
    use kind_tools
    use constants_tools
    use input_tools
    use orbital_tools, only: basis
    use block_tools, only: block_diag_CS,block_CS
    use essential_states, only: essential_state
    use fields, only: pulse
    use propagator_CN
    use stdlib_math, only: arange
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    use omp_lib, only: omp_get_wtime
    implicit none

    real(dp), allocatable :: t(:)
    integer :: N_t
    type(pulse), allocatable :: pulses(:)
    logical :: pol_active(3)
    real(dp), allocatable :: A_t(:,:),E_t(:,:)

    integer :: n_states
    complex(dp), allocatable :: psi(:), psi_new(:)
    class(essential_state), allocatable :: ess_states(:)
    complex(dp), allocatable :: projs(:,:)

    integer :: N_t_vecs, i_vecs(2)
    complex(dp), allocatable :: psi_t(:,:)

    type(block_diag_CS) :: S

    real(dp) :: t_1,t_2

contains
    subroutine init_pulses()
        integer :: i,j

        allocate(pulses(n_pulses))
        pol_active = .false.
        associate(n_env => n_env_params, n_carr => n_carrier_params)
        do i = 1,n_pulses
            call pulses(i)%init(intensities(i),t_0(i),env_params(1:n_env(i),i),carrier_params(1:n_carr(i),i),&
                pol(i), env_names(i), carrier_names(i),i)
            do j = 1,3
                if (pulses(i)%pol(j)) pol_active(j) = .true.
            end do
        end do
        end associate
    end subroutine init_pulses

    subroutine init_time_grid()
        N_t = floor((time_limits(2) + dt -time_limits(1))/dt) + 1
        allocate(t(N_t))
        t = arange(time_limits(1), time_limits(2) + dt, dt)

        if (N_t /= size(t)) then
            write(stderr,*)
            write(stderr,*) "N_t and size(t) are different", N_t, size(t)
            write(stderr,*) "Aborting execution!"
            stop
        end if

        write(stdout,*)
        write(stdout,'(a)') "Time grid info:"
        write(stdout,'(a,g16.9)') 'Start: ', time_limits(1)
        write(stdout,'(a,g16.9)') 'Stop: ', time_limits(2) + dt
        write(stdout,'(a,es25.17e3)') 'dt: ', dt
        write(stdout,'(a,i0)') 'N_t: ', N_t
        write(stdout,*)

        allocate(A_t(N_t,3),E_t(N_t,3))
    end subroutine init_time_grid

    subroutine init_CN_propagator(H_0_block,S_block,dip)
        type(block_diag_CS), intent(inout) :: H_0_block
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(inout) :: dip(-1:1)

        S = S_block

        call init_CN(pulses,pol_active,dt,H_0_block,S_block,dip,full,B_subset)
    end subroutine init_CN_propagator

    subroutine init_psi(bas)
        type(basis), intent(in) :: bas

        n_states = bas%n_states
        allocate(psi(n_states),psi_new(n_states),source=(0.0_dp,0.0_dp))

        associate(psi_0 => ess_states(1))
            psi(psi_0%ptr(1):psi_0%ptr(2)) = psi_0%vector
        end associate

        allocate(projs(N_t,n_ess))

        if (store_vecs) call init_psi_t()
    end subroutine init_psi

    subroutine init_psi_t()
        i_vecs(1) = minloc(abs(t_vecs(1)-t),dim=1)
        i_vecs(2) = minloc(abs(t_vecs(2)-t),dim=1)
        N_t_vecs = i_vecs(2)-i_vecs(1) + 1

        allocate(psi_t(n_states,N_t_vecs),source=(0.0_dp,0.0_dp))
    end subroutine init_psi_t

    subroutine propagate()
        integer :: i

        t_1 = omp_get_wtime()
        write(stdout,*)
        write(stdout,'(a)') 'Propagating TDSE...'
        call compute_projections(1)
        call save_vec(1)
        do i = 1,N_t-1
            write(stdout,'(a,a,i0,a,i0)',advance='no') achar(13),'Timestep: ', i, '/', N_t
            flush(unit=stdout)
            ! write(stdout,*) abs(ess_states(1)%projection(psi,S%blocks(ess_states(1)%block),full))**2
            call CN_step(t(i),psi,psi_new,A_t(i,:),E_t(i,:))
            call swap(psi,psi_new)
            call compute_projections(i+1)
            call save_vec(i+1)
        end do
        A_t(N_t,:) = get_A_t(t(N_t))
        E_t(N_t,:) = get_E_t(t(N_t))
        write(stdout,'(a,a,i0,a,i0)') achar(13),'Timestep: ', N_t, '/', N_t

        t_2 = omp_get_wtime()
        write(stdout,'(a)') 'Done with TDSE propagation!'
        write(stdout,*) "Total time for propagation (s): ", t_2 - t_1
        write(stdout,*) "Total time for propgagation (min): ", (t_2 - t_1)/60.d0
        write(stdout,*) "Total time for propagation (h): ", (t_2 - t_1)/(60.d0*60.d0)
        write(stdout,*)
    end subroutine propagate

    subroutine swap(x,y)
        complex(dp), allocatable, intent(inout) :: x(:)
        complex(dp), allocatable, intent(inout) :: y(:)

        complex(dp), allocatable :: temp(:)

        call move_alloc(x,temp)
        call move_alloc(y,x)
        call move_alloc(temp,y)
    end subroutine swap

    subroutine compute_projections(i)
        integer, intent(in) :: i

        integer :: j

        associate(es => ess_states)
        do j = 1,n_ess
            projs(i,j) = es(j)%projection(psi,S%blocks(es(j)%block),full)
        end do
        end associate
    end subroutine compute_projections

    subroutine save_vec(i)
        integer, intent(in) :: i

        integer :: j

        if (store_vecs) then
            if ((i>=i_vecs(1)).and.(i<=i_vecs(2))) then
                j = i-i_vecs(1) + 1
                call zcopy(n_states,psi,1,psi_t(:,j),1)
            end if
        end if
    end subroutine save_vec

    subroutine write_results(res_dir)
        character(len=*), intent(in) :: res_dir

        call write_populations(res_dir)
        call write_vecs(res_dir)
        call write_fields(res_dir)
    end subroutine write_results

    subroutine write_populations(res_dir)
        character(len=*), intent(in) :: res_dir

        integer :: i,j,unit

        open(file = res_dir//"populations.out", newunit = unit, action = 'write')

        write(unit,'(a)') "# Populations for every timestep, first column time, the rest populations"
        do i = 1,N_t
            write(unit,'(es25.17e3)', advance = 'no') t(i)
            do j = 1,n_ess
                write(unit,'(es25.17e3)', advance='no') abs(projs(i,j))**2
            end do
            write(unit,'()')
        end do

        close(unit)
    end subroutine write_populations

    subroutine write_vecs(res_dir)
        character(len=*), intent(in) :: res_dir

        integer :: unit,i

        if (store_vecs) then
            open(file = res_dir//"vecs.dat", newunit = unit, action = 'write', form = 'unformatted')
            write(unit) n_states
            write(unit) N_t_vecs
            write(unit) 1

            write(unit) psi_t

            close(unit)

            open(file = res_dir//"vecs_t.out", newunit = unit, action = 'write')
            write(unit,'(a)') "# Times at which the wavefunction was stored during propagation, then A(t) and E(t)"
            do i = i_vecs(1),i_vecs(2)
                write(unit,'(es25.17e3)',advance='no') t(i)
                write(unit,'(a,es25.17e3,a,es25.17e3,a,es25.17e3)',advance='no') ' ', A_t(i,1),' ', A_t(i,2),' ', A_t(i,3)
                write(unit,'(a,es25.17e3,a,es25.17e3,a,es25.17e3)') ' ', E_t(i,1),' ', E_t(i,2),' ', E_t(i,3)
            end do
            close(unit)
        end if

    end subroutine write_vecs

    subroutine write_fields(res_dir)
        character(len=*), intent(in) :: res_dir

        integer :: unit,i

        open(file = res_dir//"fields.out", newunit = unit, action = 'write')
        write(unit,'(a)') "# t, A_x(t), A_y(t), A_z(t), E_x(t), E_y(t), E_z(t)"
        do i = 1,N_t
            write(unit,'(es25.17e3)',advance='no') t(i)
            write(unit,'(a,es25.17e3,a,es25.17e3,a,es25.17e3)',advance='no') ' ', A_t(i,1),' ', A_t(i,2),' ', A_t(i,3)
            write(unit,'(a,es25.17e3,a,es25.17e3,a,es25.17e3)') ' ', E_t(i,1),' ', E_t(i,2),' ', E_t(i,3)
        end do
        close(unit)
    end subroutine write_fields
end module time_propagation