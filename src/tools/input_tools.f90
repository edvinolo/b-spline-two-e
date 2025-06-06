module input_tools
    use kind_tools
    use dir_tools, only: is_dir
    use iso_fortran_env, only: stdout => output_unit ,stderr => error_unit
    implicit none

    ! Basis input variables
    character(len=:), allocatable :: basis_root_dir
    character(len=64) :: basis_output_dir

    integer :: k,m,Z,max_k,k_GL,max_L,max_l_1p
    double precision :: h_max,r_max,r_2_max

    integer :: CAP_order
    double precision :: CAP_r_0
    double complex :: CAP_eta

    character :: gauge
    logical :: z_pol
    logical :: full

    namelist /basis_input/ &
    & basis_output_dir,&
    & k,&
    & m,&
    & Z,&
    & h_max,&
    & r_max,&
    & r_2_max,&
    & k_GL,&
    & CAP_order,&
    & CAP_r_0,&
    & CAP_eta,&
    & max_L,&
    & max_l_1p,&
    & max_k,&
    & gauge,&
    & z_pol,&
    & full

    ! Quasienergy input variables
    character(len=:), allocatable :: basis_dir
    character(len=64) :: basis_input_dir
    character(len=:), allocatable :: quasi_root_dir
    character(len=64) :: quasi_output_dir

    character(len=:), allocatable :: calc
    character(len=64) :: calc_type

    real(dp) :: calc_param(2)
    real(dp) :: calc_start
    real(dp) :: calc_end
    real(dp) :: other_param
    integer(int32) :: n_calc
    integer :: n_quasi
    integer :: n_blocks(2)

    complex(dp) :: gs_energy

    logical :: store_vecs
    logical :: track_proj

    logical :: direct_solver

    logical :: block_precond
    character(len=64) :: block_precond_type
    logical :: couple_pq
    integer :: n_precond
    integer, allocatable :: precond_blocks(:), precond_blocks_in(:)

    logical :: reduce_basis
    integer :: n_relevant
    integer, allocatable :: relevant_blocks(:)

    namelist /quasi_input/ &
    & basis_input_dir,&
    & quasi_output_dir,&
    & calc_type,&
    & calc_param,&
    & other_param,&
    & n_calc,&
    & n_quasi,&
    & n_blocks,&
    & gs_energy,&
    & store_vecs,&
    & track_proj,&
    & direct_solver,&
    & block_precond,&
    & block_precond_type,&
    & couple_pq,&
    & n_precond,&
    & precond_blocks_in,&
    & reduce_basis,&
    & n_relevant,&
    & relevant_blocks

    ! Diagonalization input variables
    ! Some are shared with quasienergies so are not declared here
    character(len=:), allocatable :: diag_root_dir
    character(len=64) :: diag_output_dir

    integer :: n_syms
    integer, allocatable :: blocks(:)
    integer, allocatable :: n_eigs(:)
    complex(dp), allocatable :: target_energies(:)

    namelist /diag_input/ &
    & basis_input_dir,&
    & diag_output_dir,&
    & n_syms,&
    & blocks,&
    & n_eigs,&
    & target_energies

    integer, private :: i
contains
    function get_input_file() result(res)
        character(len=:), allocatable :: res

        integer :: n_command
        character(len=256) :: input_file

        n_command = command_argument_count()
        if(n_command>1) then
            write(stderr,*)
            write(stderr,*) "More than 1 command line arguments supplied. Only the first one is used, the rest is ignored."
            write(stderr,*)
        else if (n_command == 0) then
            write(stderr,*)
            write(stderr,*) "No command line arguments supplied. You need to tell the program the name of your input file."
            write(stderr,*)
            stop
        end if

        call get_command_argument(1,input_file)
        res = input_file
    end function get_input_file

    subroutine get_basis_input(input_file)
        character(len=*), intent(in) :: input_file

        logical :: bad_input

        open(unit=1,file=input_file,action='read')
        read(1,nml=basis_input)
        close(1)

        bad_input = .false.

        basis_root_dir = trim(basis_output_dir)
        if (.not.is_dir(basis_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The basis_output_dir you supplied: ", basis_root_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the basis_output_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if ((gauge /= 'l') .and. (gauge /= 'v')) then
            write(stderr,*)
            write(stderr,*) "Bad gauge input parameter: ", gauge
            write(stderr,*) "Must be l (length) or v (velocity)"
            write(stderr,*)
            bad_input = .true.
        end if

        write(stdout,*)
        write(stdout,*) "Basis setup input:"
        write(stdout,nml=basis_input)

        if (bad_input) call error_bad_input()
    end subroutine get_basis_input

    subroutine get_quasi_input(input_file)
        character(len=*), intent(in) :: input_file

        logical :: bad_input

        open(unit=1,file=input_file,action='read')
        read(1,nml=quasi_input)
        close(1)

        bad_input = .false.

        basis_dir = trim(basis_input_dir)
        if (.not.is_dir(basis_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The basis_input_dir you supplied: ", basis_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the basis_input_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        quasi_root_dir = trim(quasi_output_dir)
        if (.not.is_dir(quasi_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The quasi_output_dir you supplied: ", quasi_root_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the quasi_output_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(calc_param<0)) then
            write(stderr,*)
            write(stderr,*) "Error! The calc_param not set: ", calc_param
            write(stderr,*) "You need to supply an positive calc_param in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (other_param<0) then
            write(stderr,*)
            write(stderr,*) "Error! The other_param not set: ", other_param
            write(stderr,*) "You need to supply an positive other_param in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (direct_solver.eqv.block_precond) then
            write(stderr,*)
            write(stderr,*) "Error! The options direct_solver and block_precond must differ: ", direct_solver, block_precond
            write(stderr,*) "You need to supply correct values in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (reduce_basis.and.block_precond.and.(block_precond_type == "PQ")) then
            if (n_precond > n_relevant) then
                write(stderr,*)
                write(stderr,*) "Error! n_precond is greater than n_relevant: ", n_precond, n_relevant
                write(stderr,*) "You need to supply correct values in the input file"
                write(stderr,*)
                bad_input = .true.
            end if

            do i = 1,n_precond
                if (.not.any(precond_blocks_in(i)==relevant_blocks)) then
                    write(stderr,*)
                    write(stderr,*) "Error! precond_blocks_in(i) not present in relevant_blocks: ",&
                                    i, precond_blocks_in(i), relevant_blocks
                    write(stderr,*) "You need to supply correct values in the input file"
                    write(stderr,*)
                    bad_input = .true.
                end if
            end do
        end if

        if (reduce_basis) then
            relevant_blocks = relevant_blocks(1:n_relevant)
        else
            relevant_blocks = relevant_blocks(1:1)
        end if

        if ((block_precond).and.(block_precond_type == "PQ")) then
            precond_blocks_in = precond_blocks_in(1:n_precond)
        else
            precond_blocks_in = precond_blocks_in(1:1)
        end if

        write(stdout,*)
        write(stdout,*) "Quasienergy input:"
        write(stdout,nml=quasi_input)

        if (bad_input) call error_bad_input()
    end subroutine get_quasi_input

    subroutine get_diag_input(input_file)
        character(len=*), intent(in) :: input_file

        logical :: bad_input

        open(unit=1,file=input_file,action='read')
        read(1,nml=diag_input)
        close(1)

        bad_input = .false.

        basis_dir = trim(basis_input_dir)
        if (.not.is_dir(basis_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The basis_input_dir you supplied: ", basis_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the basis_input_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        diag_root_dir = trim(diag_output_dir)
        if (.not.is_dir(diag_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The diag_output_dir you supplied: ", diag_root_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the diag_output_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(calc_param<0)) then
            write(stderr,*)
            write(stderr,*) "Error! The calc_param not set: ", calc_param
            write(stderr,*) "You need to supply an positive calc_param in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (n_syms<1) then
            write(stderr,*)
            write(stderr,*) "Error! n_syms < 1: ", other_param
            write(stderr,*) "You need to supply an positive n_syms >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        blocks = blocks(1:n_syms)
        n_eigs = n_eigs(1:n_syms)
        target_energies = target_energies(1:n_syms)

        write(stdout,*)
        write(stdout,*) "Diagonalization input:"
        write(stdout,nml=diag_input)

        if (bad_input) call error_bad_input()
    end subroutine get_diag_input

    subroutine write_basis_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'basis_input.dat', newunit = unit, action = 'write')
        write(unit, nml = basis_input)
        close(unit)
    end subroutine write_basis_input

    subroutine write_quasi_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'quasi_input.dat', newunit = unit, action = 'write')
        write(unit, nml = quasi_input)
        close(unit)
    end subroutine write_quasi_input

    subroutine write_diag_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'diag_input.dat', newunit = unit, action = 'write')
        write(unit, nml = diag_input)
        close(unit)
    end subroutine write_diag_input

    subroutine set_basis_defaults()
        k = 6
        m = 3
        Z = 2
        h_max = 0.5_dp
        r_max = 15.0_dp
        r_2_max = -1_dp
        k_GL = k + 6
        CAP_order = 2
        CAP_r_0 = 10.0_dp
        CAP_eta = dcmplx(1.d-3,0.d0)
        max_L = 2
        max_l_1p = 5
        max_k = 4
        z_pol = .true.
        full = .true.
    end subroutine set_basis_defaults

    subroutine set_quasi_defaults()
        basis_input_dir = "-"
        quasi_output_dir = "-"
        calc_type = "-"
        calc_param = [-1.0_dp,-1.0_dp]
        n_calc = 1
        n_quasi = 1
        other_param = -1.0_dp
        n_blocks = [1,0]
        gs_energy = (0.0_dp,0.0_dp)
        store_vecs = .false.
        track_proj = .true.
        direct_solver = .true.
        block_precond = .false.
        block_precond_type = "Jacobi"
        couple_pq = .false.
        n_precond = 200
        allocate(precond_blocks_in(n_precond),source = -1)
        reduce_basis = .false.
        n_relevant = 200
        allocate(relevant_blocks(n_relevant),source = -1)
    end subroutine set_quasi_defaults

    subroutine set_diag_defaults()
        basis_input_dir = "-"
        diag_output_dir = "-"
        n_syms = 200
        allocate(blocks(n_syms),source = -1)
        allocate(n_eigs(n_syms),source = -1)
        allocate(target_energies(n_syms),source = (0.0_dp,0.0_dp))
    end subroutine set_diag_defaults

    subroutine error_bad_input()
        write(stderr,*)
        write(stderr,*) "Some of the input parameters were bad, see warnings above"
        stop
    end subroutine error_bad_input
end module input_tools