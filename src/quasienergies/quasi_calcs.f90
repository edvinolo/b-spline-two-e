module quasi_calcs
    use kind_tools
    use constants_tools
    use stdlib_math, only: linspace, logspace
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS, block_diag_CS
    use orbital_tools
    use input_tools
    use PARDISO_tools, only: PARDISO_solver
    use ILU0_tools, only: ILU0
    use precond_tools, only: block_PC,block_PQ_PC,block_Jacobi_PC
    use GMRES_tools, only: zFGMRES
#if defined(WITH_FEAST)
    use eig_tools, only: drive_ARPACK_SI, FEAST
#else
    use eig_tools, only: drive_ARPACK_SI
#endif
    use quasi_circ
    use quasi_floquet
    use static, only: setup_H_static
    use omp_lib, only: omp_get_wtime
    implicit none

    ! Misc. integer variables
    integer :: n_states
    integer :: n_temp ! How many eigs to find when doing follow calcs

    ! Hamiltonian
    type(CSR_matrix) :: H
    type(CSR_matrix) :: S
    type(block_CS) :: H_block

    ! PARDISO solver object for shift and invert
    type(PARDISO_solver) :: solver

    ! zFGMRES solver object for iterative shift and invert
    type(zFGMRES) :: GMRES

    ! block_PC object for preconditioning GMRES
    class(block_PC), allocatable :: precond

    ! ILU0 object for test
    ! type(ILU0) :: PC_ILU0

    ! Result variables
    complex(dp), allocatable :: eigs(:,:)
    complex(dp), allocatable :: vecs(:,:,:)
    complex(dp), allocatable :: eigs_i(:)
    complex(dp), allocatable :: vecs_i(:,:)

    ! Timing variables
    real(dp) :: t_1_quasi,t_2_quasi

    ! Used for debugging the solver implementations. Could be removed when comfortable that things are working
    ! real(dp), allocatable :: dtest_vec(:,:)
    ! complex(dp), allocatable :: ztest_vec(:)

    complex(dp), external :: zdotu
contains
    subroutine omega_scan(H_0,S_block,D,bas,shift,res_dir)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        complex(dp), intent(in) :: shift
        character(len=*), intent(in) :: res_dir

        real(dp) :: V_0
        real(dp) :: omega(n_calc)
        real(dp) :: intensity

        integer :: i

        call start_timer()

        omega = linspace(calc_param(1),calc_param(2),n_calc)
        intensity = other_param

        call allocate_result(bas)
        if (track_proj) call allocate_temp_results()

        call S_block%to_CS(S,.false.)

        do i = 1,n_calc
            call print_params(i,intensity,omega(i),shift)
            V_0 = field_strength(intensity,omega(i))
            call compute_quasi(H_0,S_block,D,bas,i,V_0,omega(i),shift,n_quasi,eigs(:,i),vecs(:,:,i))
            if ((track_proj).and.(i>1)) call reorder_proj(i)
        end do

        call cleanup_solvers()

        call write_eigs(res_dir)
        call write_omega(res_dir,omega)
        call end_timer()
    end subroutine omega_scan

    subroutine intensity_scan(H_0,S_block,D,bas,shift,res_dir)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        complex(dp), intent(in) :: shift
        character(len=*), intent(in) :: res_dir

        real(dp) :: V_0
        real(dp) :: intensity(n_calc)
        real(dp) :: omega

        integer :: i

        call start_timer()

        intensity = logspace(log10(calc_param(1)),log10(calc_param(2)),n_calc)
        omega = other_param

        call allocate_result(bas)
        if (track_proj) call allocate_temp_results()

        call S_block%to_CS(S,.false.)

        do i = 1,n_calc
            call print_params(i,intensity(i),omega,shift)
            V_0 = field_strength(intensity(i),omega)
            call compute_quasi(H_0,S_block,D,bas,i,V_0,omega,shift,n_quasi,eigs(:,i),vecs(:,:,i))
            if ((track_proj).and.(i>1)) call reorder_proj(i)
        end do

        call cleanup_solvers()

        call write_eigs(res_dir)
        call write_intensity(res_dir,intensity)
        call end_timer()
    end subroutine intensity_scan

    subroutine static_scan(H_0,S_block,D,bas,shift,res_dir)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        complex(dp), intent(in) :: shift
        character(len=*), intent(in) :: res_dir

        real(dp) :: V_0(n_calc)
        real(dp) :: intensity(n_calc)
        real(dp) :: omega

        integer :: i

        call start_timer()

        intensity = logspace(log10(calc_param(1)),log10(calc_param(2)),n_calc)
        omega = other_param

        call allocate_result(bas)
        if (track_proj) call allocate_temp_results()

        call S_block%to_CS(S,.false.)

        do i = 1,n_calc
            call print_params(i,intensity(i),omega,shift)
            V_0(i) = field_strength(intensity(i),omega)
            call compute_quasi(H_0,S_block,D,bas,i,V_0(i),omega,shift,n_quasi,eigs(:,i),vecs(:,:,i))
            if ((track_proj).and.(i>1)) call reorder_proj(i)
        end do

        call cleanup_solvers()

        call write_eigs(res_dir)
        call write_static(res_dir,V_0)
        call end_timer()
    end subroutine static_scan

    subroutine omega_follow(H_0,S_Block,D,bas,shift,res_dir)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        complex(dp), intent(in) :: shift
        character(len=*), intent(in) :: res_dir

        real(dp) :: V_0
        real(dp) :: omega(n_calc)
        real(dp) :: intensity

        integer :: i,j,max_i

        call start_timer()

        omega = linspace(calc_param(1),calc_param(2),n_calc)
        intensity = other_param

        call allocate_result(bas)
        call allocate_temp_results()

        call S_block%to_CS(S,.false.)

        do i = 1,n_calc
            V_0 = field_strength(intensity,omega(i))

            if (i == 1) then 
                call print_params(i,intensity,omega(i),shift)
                call compute_quasi(H_0,S_block,D,bas,i,V_0,omega(i),shift,n_quasi,eigs(:,i),vecs(:,:,i))
            else
                do j = 1,n_quasi
                    call print_params(i,intensity,omega(i),eigs(j,i-1),j)
                    call compute_quasi(H_0,S_block,D,bas,i,V_0,omega(i),eigs(j,i-1),n_temp,eigs_i(:),vecs_i(:,:),vecs(:,j,i-1))
                    max_i =  find_max_proj(vecs(:,j,i-1))
                    eigs(j,i) = eigs_i(max_i)
                    vecs(:,j,i) = vecs_i(:,max_i)
                end do
            end if
        end do

        call cleanup_solvers()

        call write_eigs(res_dir)
        call write_omega(res_dir,omega)
        call end_timer()
    end subroutine omega_follow

    subroutine intensity_follow(H_0,S_Block,D,bas,shift,res_dir)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        complex(dp), intent(in) :: shift
        character(len=*), intent(in) :: res_dir

        real(dp) :: V_0
        real(dp) :: intensity(n_calc)
        real(dp) :: omega

        integer :: i,j,max_i

        call start_timer()

        intensity = logspace(log10(calc_param(1)),log10(calc_param(2)),n_calc)
        omega = other_param

        call allocate_result(bas)
        call allocate_temp_results()

        call S_block%to_CS(S,.false.)

        do i = 1,n_calc
            V_0 = field_strength(intensity(i),omega)

            if (i == 1) then
                call print_params(i,intensity(i),omega,shift) 
                call compute_quasi(H_0,S_block,D,bas,i,V_0,omega,shift,n_quasi,eigs(:,i),vecs(:,:,i))
            else
                do j = 1,n_quasi
                    call print_params(i,intensity(i),omega,eigs(j,i-1),j=j)
                    call compute_quasi(H_0,S_block,D,bas,i,V_0,omega,eigs(j,i-1),n_temp,eigs_i(:),vecs_i(:,:),vecs(:,j,i-1))
                    max_i =  find_max_proj(vecs(:,j,i-1))
                    eigs(j,i) = eigs_i(max_i)
                    vecs(:,j,i) = vecs_i(:,max_i)
                end do
            end if
        end do

        call cleanup_solvers()

        call write_eigs(res_dir)
        call write_intensity(res_dir,intensity)
        call end_timer()

    end subroutine intensity_follow

    subroutine compute_quasi(H_0,S_Block,D,bas,i,V_0_i,omega_i,shift_i,n_eigs_out,eigs_out,vecs_out,v0)
        type(block_diag_CS), intent(in) :: H_0
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(in) :: D
        type(basis), intent(in) :: bas
        integer, intent(in) :: i
        real(dp), intent(in) :: V_0_i
        real(dp), intent(in) :: omega_i
        complex(dp), intent(in) :: shift_i
        integer, intent(in) :: n_eigs_out
        complex(dp), intent(out) :: eigs_out(:)
        complex(dp), intent(out) :: vecs_out(:,:)
        complex(dp), optional, intent(in) :: v0(:)

        if (z_pol) then
            if (calc_type == 'static') then
                call setup_H_static(H_block, bas, H_0, S_block, D, shift_i, V_0_i)
            else
            ! H_0_shift = setup_H_0_floquet()
            end if
        else
            call setup_H_circ(H_block, bas, H_0, S_block, D, omega_i, shift_i, V_0_i)
        end if

        if (block_precond) then
            if (i == 1) then
                if (block_precond_type == 'PQ') then
                    allocate(block_PQ_PC :: precond)
                else if (block_precond_type == 'Jacobi') then
                    allocate(block_Jacobi_PC :: precond)
                else
                    write(stderr,*) ""
                    write(stderr,'(a,a)') "Error in Block_PC setup, unknown block_precond_type: ", block_precond_type
                    error stop
                end if
                call precond%setup(H_block)
            else
                call precond%update(H_block)
            end if
        end if

        call H_block%to_CS(H,.true.)

        if (direct_solver) then
            if (i == 1) call solver%setup(H%shape(1),H%nnz,H%data,H%index_ptr,H%indices)
            call solver%factor(H%data,H%index_ptr,H%indices)
            if (present(v0)) then
                call drive_ARPACK_SI(H,solver,S,full,shift_i,n_eigs_out,eigs_out(:),vecs_out(:,:),v0=v0)
            else
                call drive_ARPACK_SI(H,solver,S,full,shift_i,n_eigs_out,eigs_out(:),vecs_out(:,:))
            end if

        else if (block_precond) then
            if (i==1) then
                call GMRES%setup(H,full)
            else
                call GMRES%update(H)
            end if

            if (present(v0)) then
                call drive_ARPACK_SI(GMRES,precond,S,full,shift_i,n_eigs_out,eigs_out(:),vecs_out(:,:),v0=v0)
            else
                call drive_ARPACK_SI(GMRES,precond,S,full,shift_i,n_eigs_out,eigs_out(:),vecs_out(:,:))
            end if
        end if

    end subroutine compute_quasi

    subroutine cleanup_solvers()
        if (direct_solver) then
            call solver%cleanup()
        end if

        if (block_precond) then
            call GMRES%cleanup()
            call precond%cleanup()
        end if
    end subroutine cleanup_solvers

    function field_strength(intensity,omega) result(res)
        real(dp), intent(in) :: intensity
        real(dp), intent(in) :: omega
        real(dp) :: res

        res = 0
        if (gauge == 'l') res = E_0_au(intensity)
        if (gauge == 'v') res = A_0_au(intensity,omega)

        if (.not.z_pol) res = res/sqrt_2
    end function field_strength

    subroutine allocate_result(bas)
        type(basis), intent(in) :: bas

        n_states = compute_dim(bas)

        allocate(eigs(n_quasi,n_calc))
        allocate(vecs(n_states,n_quasi,n_calc), source = (0.0_dp,0.0_dp))

    end subroutine allocate_result

    subroutine allocate_temp_results()
        if ((calc_type == 'omega_follow').or.(calc_type == 'intensity_follow')) then
            n_temp = 2
        else
            n_temp = n_quasi
        end if
        allocate(eigs_i(n_temp),vecs_i(n_states,n_temp))
    end subroutine allocate_temp_results

    function compute_dim(bas) result(res)
        type(basis), intent(in) :: bas
        integer :: res

        res = 0
        if (z_pol) then
            !TODO
            if (calc_type == 'static') then
                res = sum(bas%syms%n_config)
            end if
        else
            res = sum(bas%syms%n_config)
        end if


    end function compute_dim

    subroutine write_eigs(res_dir)
        character(len=*), intent(in) :: res_dir

        ! character(len=200) :: fmt
        integer :: i,unit,j
        open(file = res_dir//"energies.out", newunit = unit, action = 'write')

        ! write(fmt,'(a,i0)')'(',
        ! fmt = trim(fmt//'(es22.15,es22.15))')
        ! write(fmt,'(a,i0,a)') '(', n_quasi,'("("es24.17,es24.17"j) "))' ! This should be readable by numpy.loadtxt
        ! fmt = '(a,es24.17,es24.17,a)'
        write(unit,'(a)') "# Energies in [a.u], each row is one point in parameter space"
        do i = 1,n_calc
            ! write(unit,fmt) eigs(:,i)
            do j = 1,n_quasi
                write(unit,'(a)', advance='no') '('

                if (real(eigs(j,i),kind = dp)<0) then
                    write(unit,'(es24.17e2)', advance='no') real(eigs(j,i),kind = dp)
                else
                    write(unit,'(es23.17e2)', advance='no') real(eigs(j,i),kind = dp)
                end if

                if (aimag(eigs(j,i))<0) then
                    write(unit,'(es24.17e2)', advance='no') aimag(eigs(j,i))
                else
                    write(unit,'(a)', advance='no') '+'
                    write(unit,'(es23.17e2)', advance='no') aimag(eigs(j,i))
                end if

                write(unit,'(a)', advance='no') 'j) '
            end do
            write(unit,*)''
        end do

        close(unit)
    end subroutine write_eigs

    subroutine write_omega(res_dir,omega)
        character(len=*), intent(in) :: res_dir
        real(dp), intent(in) :: omega(:)

        integer :: i,unit
        open(file = res_dir//"omega.out", newunit = unit, action = 'write')

        write(unit,'(a)') "# omega [a.u]"
        do i = 1,n_calc
            write(unit,'(es24.17)') omega(i)
        end do

        close(unit)
    end subroutine write_omega

    subroutine write_intensity(res_dir,intensity)
        character(len=*), intent(in) :: res_dir
        real(dp), intent(in) :: intensity(:)

        integer :: i,unit
        open(file = res_dir//"intensity.out", newunit = unit, action = 'write')

        write(unit,'(a)') "# intensity [W/cm^2]"
        do i = 1,n_calc
            write(unit,'(es24.17)') intensity(i)
        end do

        close(unit)
    end subroutine write_intensity

    subroutine write_static(res_dir,V_0)
        character(len=*), intent(in) :: res_dir
        real(dp), intent(in) :: V_0(:)

        integer :: i,unit
        open(file = res_dir//"static_field.out", newunit = unit, action = 'write')

        write(unit,'(a)') "# V_0 [a.u.]"
        do i = 1,n_calc
            write(unit,'(es24.17)') V_0(i)
        end do

        close(unit)
    end subroutine write_static

    subroutine print_params(iteration,intensity,omega,shift,j)
        integer, intent(in) :: iteration
        real(dp), intent(in) :: intensity
        real(dp), intent(in) :: omega
        complex(dp), intent(in) :: shift
        integer, optional, intent(in) :: j

        real(dp) :: U_p

        U_p = A_0_au(intensity,omega)**2/4

        write(6,*)
        write(6,'(a,i0,a,i0)') 'Iteration: ', iteration, '/', n_calc

        if (present(j)) then
            write(6,'(a,i0,a,i0)') 'Eigenvalue: ', j, '/', n_quasi
        end if

        write(6,'(a,es15.4,es15.4,a)') 'Shift: ', shift, ' [a.u.]'
        write(6,'(a,f10.7,a,f10.7,a)') 'Omega: ', omega, ' [a.u.], ', omega*au_to_eV, ' [eV]'
        write(6,'(a,es15.4,a)') 'Intensity: ', intensity, ' [W/cm2]'
        write(6,'(a,es15.4,a)') 'E_0: ', E_0_au(intensity), ' [a.u.]'
        write(6,'(a,es15.4,a)') 'A_0: ', A_0_au(intensity,omega), ' [a.u.]'
        write(6,'(a,es15.4,a,es15.4,a)') 'U_p: ', U_p, ' [a.u.], ', U_p*au_to_eV, ' [eV]'
        write(6,*)
    end subroutine print_params

    subroutine reorder_blocks(bas,H_0_block,S_block,D)
        type(basis), intent(inout) :: bas
        type(block_diag_CS), intent(inout) :: H_0_block
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(inout) :: D

        integer :: j,i

        ! write(6,*) bas%syms%m
        ! do i = 1,bas%n_sym
        !     write(6,*) D%blocks(i,:)%nnz
        ! end do

        do i = 1,n_precond
            j = precond_blocks(i)
            if (j == i) cycle

            ! Reorder the symmetry blocks
            bas%syms([i,j]) = bas%syms([j,i])

            ! Transform block diagonal matrices
            H_0_block%blocks([i,j]) = H_0_block%blocks([j,i])
            S_block%blocks([i,j]) = S_block%blocks([j,i])

            ! Transform off-diagonal blocks
            D%blocks([i,j],:) = D%blocks([j,i],:)
            D%blocks(:,[i,j]) = D%blocks(:,[j,i])
        end do

        ! If using upper storage, make sure that only blocks in the upper triangular part are present
        ! Transpose those blocks that are below diagonal
        if (.not.full) then
            do j = 1,bas%n_sym
                do i = j+1,bas%n_sym
                    if (D%blocks(i,j)%nnz /= 0) then
                        call D%blocks(i,j)%transp(D%blocks(j,i))
                        call D%blocks(i,j)%deall()
                    end if
                end do
            end do
        end if

        ! write(6,*) bas%syms%m
        ! do i = 1,bas%n_sym
        !     write(6,*) D%blocks(i,:)%nnz
        ! end do

    end subroutine reorder_blocks

    subroutine red_basis(bas,H_0_block,S_Block,D)
        type(basis), intent(inout) :: bas
        type(block_diag_CS), intent(inout) :: H_0_block
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(inout) :: D

        integer :: j,i,ptr

        type(CSR_matrix), allocatable :: H_temp(:),S_temp(:)
        type(CSR_matrix), allocatable :: D_temp(:,:)

        allocate(H_temp(n_relevant),S_temp(n_relevant),D_temp(n_relevant,n_relevant))

        bas%n_sym = n_relevant
        bas%syms = bas%syms(relevant_blocks)

        if (block_precond.and.(block_precond_type == 'PQ')) then
            allocate(precond_blocks(n_precond))
            ptr = 1
            do i = 1,n_relevant
                if (relevant_blocks(i) == precond_blocks_in(ptr)) then
                    precond_blocks(ptr) = i
                    ptr = ptr + 1
                end if
            end do
        end if

        H_0_block%block_shape = [n_relevant,n_relevant]
        do i = 1,n_relevant
            H_temp(i) = H_0_block%blocks(relevant_blocks(i))
        end do
        call move_alloc(H_temp,H_0_block%blocks)
        call H_0_block%compute_shape()

        S_block%block_shape = [n_relevant,n_relevant]
        do i = 1,n_relevant
            S_temp(i) = S_block%blocks(relevant_blocks(i))
        end do
        call move_alloc(S_temp,S_block%blocks)
        call S_block%compute_shape()

        D%block_shape = [n_relevant,n_relevant]
        do j = 1,n_relevant
            do i = 1,n_relevant
                D_temp(i,j) = D%blocks(relevant_blocks(i),relevant_blocks(j))
            end do
        end do
        call move_alloc(D_temp,D%blocks)
        call D%compute_shape()
    end subroutine red_basis

    ! Finds the vector in vecs_i with maximum S projection on vec.
    ! Puts the result in vecs_i(:,1) and the corresponding eigenvalue in eigs_i(1)
    function find_max_proj(vec) result(res)
        complex(dp), intent(in) :: vec(:)
        integer :: res

        integer :: n,i,max_i
        real(dp) :: proj, max_proj
        complex(dp), allocatable :: temp(:)

        max_i = -1
        max_proj = -1.0_dp
        n = size(eigs_i)

        allocate(temp(size(vec)))

        if (full) then
            call CSR_mv(S,vec,temp)
        else
            call CSR_mv_sym(S,vec,temp)
        end if

        do i = 1,n
            proj = abs(zdotu(n_states,vecs_i(:,i),1,temp,1))
            write(stdout,'(a,es11.4,es11.4,a,es11.4)') 'Eig: ', real(eigs_i(i),kind=dp),aimag(eigs_i(i)), ', proj: ', proj
            if (proj > max_proj) then
                max_i = i
                max_proj = proj
            end if
        end do

        res = max_i
    end function find_max_proj

    subroutine reorder_proj(i)
        integer, intent(in) :: i

        integer :: j,kk,max
        real(dp) :: proj,max_proj
        complex(dp), allocatable :: temp(:)
        logical :: free(n_quasi)

        free = .true.
        allocate(temp(n_states))

        do j=1,n_quasi
            if (full) then
                call CSR_mv(S,vecs(:,j,i),temp)
            else
                call CSR_mv_sym(S,vecs(:,j,i),temp)
            end if

            max_proj = -1.0_dp
            max = -1
            do kk=1,n_quasi
                proj = abs(zdotu(n_states,vecs(:,kk,i-1),1,temp,1))
                if (proj > max_proj) then
                    max = kk
                    max_proj = proj
                end if
            end do

            if (free(max)) then
                free(max) = .false.
                vecs_i(:,max) = vecs(:,j,i)
                eigs_i(max) = eigs(j,i)
            else
                free(j) = .false. ! Don't move it if already found vector that wants the max slot (TODO: choose the one with largest projection to go in max slot)
                vecs_i(:,j) = vecs(:,j,i)
                eigs_i(j) = eigs(j,i)
            end if
        end do

        vecs(:,:,i) = vecs_i
        eigs(:,i) = eigs_i

        write(stdout,*) ''
        write(stdout,*) 'Reordered eigenvalues:'
        do j=1,n_quasi
            write(stdout,*) eigs(j,i)
        end do
        write(stdout,*) ''
    end subroutine reorder_proj

    subroutine start_timer()
        t_1_quasi = omp_get_wtime()
    end subroutine start_timer

    subroutine end_timer()
        t_2_quasi = omp_get_wtime()
        write(stdout,*)
        write(stdout,'(a,es11.4)') "Total time for quasienergy calculations (s): ", t_2_quasi - t_1_quasi
        write(stdout,'(a,es11.4)') "Total time for quasienergy calculations (min): ", (t_2_quasi - t_1_quasi)/60_dp
        write(stdout,'(a,es11.4)') "Total time for quasienergy calculations (h): ", (t_2_quasi - t_1_quasi)/(60_dp*60_dp)
        write(stdout,*)
    end subroutine end_timer
end module quasi_calcs