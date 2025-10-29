module input_tools
    use kind_tools
    use dir_tools, only: is_dir
    use iso_fortran_env, only: stdout => output_unit ,stderr => error_unit
    implicit none

    ! Basis input variables
    character(len=:), allocatable :: basis_root_dir
    character(len=64) :: basis_output_dir

    integer :: k,m,Z,max_k,k_GL,max_L,max_l_1p,max_l2
    double precision :: h_max,r_max,r_2_max,r_all_l

    integer :: CAP_order
    double precision :: CAP_r_0
    double complex :: CAP_eta

    character :: gauge
    logical :: z_pol
    logical :: full
    logical :: two_el

    namelist /basis_input/ &
    & basis_output_dir,&
    & k,&
    & m,&
    & Z,&
    & h_max,&
    & r_max,&
    & r_2_max,&
    & r_all_l,&
    & k_GL,&
    & CAP_order,&
    & CAP_r_0,&
    & CAP_eta,&
    & max_L,&
    & max_l_1p,&
    & max_l2,&
    & max_k,&
    & gauge,&
    & z_pol,&
    & full, &
    & two_el

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

    logical :: ILU_precond

    logical :: reduce_basis
    integer :: n_relevant
    integer, allocatable :: relevant_blocks(:)

    logical :: use_essential_states
    integer :: n_ess
    integer, allocatable :: target_blocks(:)
    integer, allocatable :: target_Floquet_blocks(:)
    complex(dp), allocatable :: targets(:)

    logical :: calc_H_eff

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
    & ILU_precond,&
    & reduce_basis,&
    & n_relevant,&
    & relevant_blocks,&
    & use_essential_states,&
    & n_ess,&
    & target_blocks,&
    & target_Floquet_blocks,&
    & targets,&
    & calc_H_eff

    ! Diagonalization input variables
    ! Some are shared with quasienergies so are not declared here
    character(len=:), allocatable :: diag_root_dir
    character(len=64) :: diag_output_dir

    integer :: n_syms
    integer, allocatable :: blocks(:)
    integer, allocatable :: n_eigs(:)
    complex(dp), allocatable :: target_energies(:)
    logical :: store_diag_vecs

    namelist /diag_input/ &
    & basis_input_dir,&
    & diag_output_dir,&
    & n_syms,&
    & blocks,&
    & n_eigs,&
    & target_energies,&
    & store_diag_vecs

    ! Function evaluation input variables
    character(len=:), allocatable :: eval_root_dir
    character(len=64) :: eval_dir

    integer(int32) :: N_r
    integer(int32) :: N_theta
    integer(int32) :: N_phi
    real(dp) :: r_limits(2)
    real(dp) :: theta_limits(2)
    real(dp) :: phi_limits(2)

    namelist /eval_input/ &
    & eval_dir,&
    & N_r,&
    & r_limits,&
    & N_theta,&
    & theta_limits,&
    & N_phi,&
    & phi_limits

    integer, private :: i

    ! Time propagation input variables (some have same name as quasienergies, so are declared there)
    character(len=:), allocatable :: time_prop_root_dir
    character(len=64) :: time_prop_output_dir

    real(dp) :: time_limits(2)
    real(dp) :: dt

    integer, parameter :: max_len_pulse_name = 64
    integer, parameter :: max_n_pulses = 20
    integer, parameter :: max_n_pulse_params = 20

    integer :: n_pulses
    character(len=max_len_pulse_name), allocatable :: env_names(:)
    character(len=max_len_pulse_name), allocatable :: carrier_names(:)
    real(dp), allocatable :: env_params(:,:)
    real(dp), allocatable :: carrier_params(:,:)
    integer, allocatable :: n_env_params(:)
    integer, allocatable :: n_carrier_params(:)
    real(dp), allocatable :: intensities(:)
    real(dp), allocatable :: t_0(:)
    character(len=1), allocatable :: pol(:)

    real(dp) :: t_vecs(2)

    namelist /time_prop_input/ &
    & basis_input_dir,&
    & time_prop_output_dir,&
    & time_limits,&
    & dt,&
    & n_pulses,&
    & env_names,&
    & n_env_params,&
    & env_params,&
    & carrier_names,&
    & n_carrier_params,&
    & carrier_params,&
    & intensities,&
    & t_0,&
    & pol,&
    & use_essential_states,&
    & n_ess,&
    & target_blocks,&
    & targets,&
    & reduce_basis,&
    & n_relevant,&
    & relevant_blocks,&
    & block_precond,&
    & block_precond_type,&
    & couple_pq,&
    & n_precond,&
    & precond_blocks_in,&
    & store_vecs,&
    & t_vecs

    ! RKH input variables
    character(len=:), allocatable :: RKH_root_dir
    character(len=64) :: RKH_output_dir

    real(dp) :: RKH_intensity
    real(dp) :: RKH_omega

    namelist /RKH_input/ &
    & basis_input_dir,&
    & RKH_output_dir,&
    & N_r,&
    & r_limits,&
    & RKH_intensity,&
    & RKH_omega,&
    & n_quasi,&
    & gs_energy

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

        if (r_2_max > r_max) then
            write(stderr,*)
            write(stderr,*) "Error! r_2_max > r_max: ", r_2_max, r_max
            write(stderr,*) "You need to supply r_2_max < r_max in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (r_all_l > r_max) then
            write(stderr,*)
            write(stderr,*) "Error! r_all_l > r_max: ", r_all_l, r_max
            write(stderr,*) "You need to supply r_all_l < r_max in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if ((r_2_max > 0).and.(r_2_max >= r_all_l)) then
            write(stderr,*)
            write(stderr,*) "Error! r_2_max >= r_all_l: ", r_2_max, r_all_l
            write(stderr,*) "You need to supply r_2_max < r_all_l in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if ((r_2_max <= 0).and.(r_all_l > 0)) then
            write(stderr,*)
            write(stderr,*) "Error! r_2_max <= 0 while r_all_l > 0: ", r_2_max, r_all_l
            write(stderr,*) "You need to supply r_2_max < 0 and r_all_l < 0, or 0 < r_2_max < r_all_l in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if ((max_l2 > max_l_1p)) then
            write(stderr,*)
            write(stderr,*) "Error! max_l2 > max_l_1p: ", max_l2, max_l_1p
            write(stderr,*) "You need to supply max_l2 < max_l_1p in the input file"
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
        integer :: solv_counter

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
            write(stderr,*) "You need to supply a positive calc_param in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (other_param<0) then
            write(stderr,*)
            write(stderr,*) "Error! The other_param not set: ", other_param
            write(stderr,*) "You need to supply a positive other_param in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        solv_counter = 0
        if (direct_solver) solv_counter = solv_counter + 1
        if (block_precond) solv_counter = solv_counter + 1
        if (ILU_precond) solv_counter = solv_counter + 1

        if (solv_counter /= 1) then
            write(stderr,*)
            write(stderr,*) "Error! One, and only one of the solver options must be set."
            write(stderr,*) "direct_solver: ", direct_solver
            write(stderr,*) "block_precond: ", block_precond
            write(stderr,*) "ILU_precond: ", ILU_precond
            write(stderr,*) "You need to supply correct values in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        call handle_basis_and_ess_states(bad_input)

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

    subroutine get_eval_input(input_file)
        character(len=*), intent(in) :: input_file

        logical :: bad_input

        open(unit=1,file=input_file,action='read')
        read(1,nml=eval_input)
        close(1)

        bad_input = .false.

        eval_root_dir = trim(eval_dir)//'/'
        if (.not.is_dir(eval_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The eval_dir you supplied: ", eval_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the eval_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(r_limits<0)) then
            write(stderr,*)
            write(stderr,*) "Error! Bad r_limits input: ", r_limits
            write(stderr,*) "You need to supply positive values for r_limits in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (N_r<1) then
            write(stderr,*)
            write(stderr,*) "Error! N_r < 1: ", N_r
            write(stderr,*) "You need to supply N_r >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        blocks = blocks(1:n_syms)
        n_eigs = n_eigs(1:n_syms)
        target_energies = target_energies(1:n_syms)

        write(stdout,*)
        write(stdout,*) "Function evaluation input:"
        write(stdout,nml=eval_input)

        if (bad_input) call error_bad_input()
    end subroutine get_eval_input

    subroutine get_time_prop_input(input_file)
        character(len=*), intent(in) :: input_file

        integer :: unit
        logical :: bad_input

        open(newunit=unit,file=input_file,action='read')
        read(unit,nml=time_prop_input)
        close(unit)

        bad_input = .false.

        basis_dir = trim(basis_input_dir)
        if (.not.is_dir(basis_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The basis_input_dir you supplied: ", basis_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the basis_input_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        time_prop_root_dir = trim(time_prop_output_dir)
        if (.not.is_dir(time_prop_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The time_prop_output_dir you supplied: ", time_prop_root_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the time_prop_output_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (dt <= 0) then
            write(stderr,*)
            write(stderr,*) "Error! dt <= 0: ", dt
            write(stderr,*) "You need to supply dt > 0 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (time_limits(1) >= time_limits(2)) then
            write(stderr,*)
            write(stderr,*) "Error! time_limits(1) >= time_limits(2): ", time_limits
            write(stderr,*) "You need to supply time_limits(1) < time_limits(2) in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (n_pulses < 0) then
            write(stderr,*)
            write(stderr,*) "Error! n_pulses < 0: ", n_pulses
            write(stderr,*) "You need to supply n_pulses >= 0 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        env_names = env_names(1:n_pulses)
        n_env_params = n_env_params(1:n_pulses)
        env_params = env_params(:,1:n_pulses)
        carrier_names = carrier_names(1:n_pulses)
        n_carrier_params = n_carrier_params(1:n_pulses)
        carrier_params = carrier_params(:,1:n_pulses)
        intensities = intensities(1:n_pulses)
        t_0 = t_0(1:n_pulses)
        pol = pol(1:n_pulses)

        if (any(n_env_params < 1)) then
            write(stderr,*)
            write(stderr,*) "Error! n_env_params < 1 found: ", n_env_params
            write(stderr,*) "You need to supply n_env_params >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(n_carrier_params < 1)) then
            write(stderr,*)
            write(stderr,*) "Error! n_carrier_params < 1 found: ", n_carrier_params
            write(stderr,*) "You need to supply n_carrier_params >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(intensities < 0)) then
            write(stderr,*)
            write(stderr,*) "Error! intensities < 0 found: ", intensities
            write(stderr,*) "You need to supply intensities >= 0 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (.not.((any(pol=='x')).or.(any(pol=='y')).or.(any(pol=='z')))) then
            write(stderr,*)
            write(stderr,*) "Error! pol /= x,y,z found: ", pol
            write(stderr,*) "The polarizaiton must be one of x,y,or z in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if ((n_ess < 1) .or. (.not.(use_essential_states))) then
            write(stderr,*)
            write(stderr,*) "Error! n_ess < 1 or use_essential_states == .false.", n_ess, use_essential_states
            write(stderr,*) "Must use one essential state as inital state"
            write(stderr,*)
            bad_input = .true.
        end if
        call handle_basis_and_ess_states(bad_input)

        if (store_vecs.and.(t_vecs(1) >= t_vecs(2))) then
            write(stderr,*)
            write(stderr,*) "Error! t_vecs(1) >= t_vecs(2)", t_vecs(1), t_vecs(2)
            write(stderr,*) "Must set t_vecs(1) < t_vecs(2) in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        write(stdout,*)
        write(stdout,*) "Time propagation input:"
        write(stdout,nml=time_prop_input)

        if (bad_input) call error_bad_input()
    end subroutine get_time_prop_input

    subroutine get_RKH_input(input_file)
        character(len=*), intent(in) :: input_file

        logical :: bad_input

        open(unit=1,file=input_file,action='read')
        read(1,nml=RKH_input)
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

        RKH_root_dir = trim(RKH_output_dir)//'/'
        if (.not.is_dir(RKH_root_dir)) then
            write(stderr,*)
            write(stderr,*) "Error! The RKH_output_dir you supplied: ", RKH_output_dir, " does not exist."
            write(stderr,*) "You need to supply an exisiting directory to the RKH_output_dir variable in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (any(r_limits<0)) then
            write(stderr,*)
            write(stderr,*) "Error! Bad r_limits input: ", r_limits
            write(stderr,*) "You need to supply positive values for r_limits in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (N_r<1) then
            write(stderr,*)
            write(stderr,*) "Error! N_r < 1: ", N_r
            write(stderr,*) "You need to supply N_r >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (RKH_intensity<0) then
            write(stderr,*)
            write(stderr,*) "Error! RKH_intensity < 0: ", RKH_intensity
            write(stderr,*) "You need to supply RKH_intensity < 0 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (RKH_omega<0) then
            write(stderr,*)
            write(stderr,*) "Error! RKH_omega < 0: ", RKH_omega
            write(stderr,*) "You need to supply RKH_omega < 0 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (N_quasi<1) then
            write(stderr,*)
            write(stderr,*) "Error! N_quasi < 1: ", N_quasi
            write(stderr,*) "You need to supply N_quasi >= 1 in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        write(stdout,*)
        write(stdout,*) "Function evaluation input:"
        write(stdout,nml=RKH_input)

        if (bad_input) call error_bad_input()
    end subroutine get_RKH_input

    subroutine handle_basis_and_ess_states(bad_input)
        logical, intent(inout) :: bad_input

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

        if (reduce_basis.and.use_essential_states) then
            do i=1,n_ess
                if (.not.any(target_blocks(i)==relevant_blocks)) then
                    write(stderr,*)
                    write(stderr,*) "Error! target_blocks(i) not present in relevant_blocks: ",&
                                    i, target_blocks(i), relevant_blocks
                    write(stderr,*) "You need to supply correct values in the input file"
                    write(stderr,*)
                    bad_input = .true.
                end if
            end do
        end if

        if (calc_H_eff.and.(.not.(use_essential_states))) then
            write(stderr,*)
            write(stderr,*) "Error! Must use_essential_states must be .true. if calc_H_eff is .true.: ", &
                                use_essential_states, calc_H_eff
            write(stderr,*) "You need to supply correct values in the input file"
            write(stderr,*)
            bad_input = .true.
        end if

        if (calc_H_eff.and.(n_ess /= n_quasi)) then
            write(stderr,*)
            write(stderr,*) "Error! n_ess and n_quasi must be equal if calc_H_eff is .true.: ", n_ess, n_quasi
            write(stderr,*) "(H_eff must be square!)"
            write(stderr,*) "You need to supply correct values in the input file"
            write(stderr,*)
            bad_input = .true.
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

        if (use_essential_states) then
            target_blocks = target_blocks(:n_ess)
            target_Floquet_blocks = target_Floquet_blocks(:n_ess)
            targets = targets(:n_ess)
        end if
    end subroutine handle_basis_and_ess_states

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

    subroutine write_eval_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'eval_input.dat', newunit = unit, action = 'write')
        write(unit, nml = eval_input)
        close(unit)
    end subroutine write_eval_input

    subroutine write_time_prop_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'time_prop_input.dat', newunit = unit, action = 'write')
        write(unit, nml = time_prop_input)
        close(unit)
    end subroutine write_time_prop_input

    subroutine write_RKH_input(res_dir)
        character(len=:), allocatable, intent(in) :: res_dir

        integer :: unit

        open(file = res_dir // 'RKH_input.dat', newunit = unit, action = 'write')
        write(unit, nml = RKH_input)
        close(unit)
    end subroutine write_RKH_input

    subroutine set_basis_defaults()
        k = 6
        m = 3
        Z = 2
        h_max = 0.5_dp
        r_max = 15.0_dp
        r_2_max = -1.0_dp
        r_all_l = -1.0_dp
        k_GL = k + 6
        CAP_order = 2
        CAP_r_0 = 10.0_dp
        CAP_eta = dcmplx(1.d-3,0.d0)
        max_L = 2
        max_l_1p = 5
        max_l2 = 5
        max_k = 4
        z_pol = .true.
        full = .true.
        two_el = .true.
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
        ILU_precond = .false.
        reduce_basis = .false.
        n_relevant = 200
        allocate(relevant_blocks(n_relevant),source = -1)
        use_essential_states = .false.
        n_ess = 200
        allocate(target_blocks(n_ess),source = -1)
        allocate(target_Floquet_blocks(n_ess),source = 0)
        allocate(targets(n_ess),source = (0.0_dp,0.0_dp))
        calc_H_eff = .false.
    end subroutine set_quasi_defaults

    subroutine set_diag_defaults()
        basis_input_dir = "-"
        diag_output_dir = "-"
        n_syms = 200
        allocate(blocks(n_syms),source = -1)
        allocate(n_eigs(n_syms),source = -1)
        allocate(target_energies(n_syms),source = (0.0_dp,0.0_dp))
        store_diag_vecs = .false.
    end subroutine set_diag_defaults

    subroutine set_eval_defaults()
        eval_dir = "-"
        N_r = 50
        r_limits = [0.01_dp,10.0_dp]
    end subroutine set_eval_defaults

    subroutine set_time_prop_defaults()
        basis_input_dir = "-"
        time_prop_output_dir = "-"
        time_limits = [0,100]
        dt = 0.05_dp
        n_pulses = -1
        allocate(env_names(max_n_pulses))
        allocate(n_env_params(max_n_pulses),source=-1)
        allocate(env_params(max_n_pulse_params,max_n_pulses))
        allocate(carrier_names(max_n_pulses))
        allocate(n_carrier_params(max_n_pulses),source=-1)
        allocate(carrier_params(max_n_pulse_params,max_n_pulses))
        allocate(intensities(max_n_pulses),source=-1.0_dp)
        allocate(t_0(max_n_pulses),source=0.0_dp)
        allocate(pol(max_n_pulses),source='-')
        use_essential_states = .true.
        n_ess = 1
        allocate(target_blocks(100),source=-1)
        target_blocks(1) = 1
        allocate(targets(100),source = (0.0_dp,0.0_dp))
        reduce_basis = .false.
        n_relevant = 200
        allocate(relevant_blocks(n_relevant), source = -1)
        block_precond = .true.
        block_precond_type = "Jacobi"
        couple_pq = .false.
        n_precond = 200
        allocate(precond_blocks_in(n_precond),source = -1)
        store_vecs = .false.
        t_vecs = [0.0d0,0.5d0]

        ! Just set this so that check passes
        allocate(target_Floquet_blocks(n_relevant),source = 0)
    end subroutine set_time_prop_defaults

    subroutine set_RKH_defaults()
        basis_input_dir = "-"
        RKH_output_dir = "-"
        N_r = -1
        r_limits = [-2.0_dp,-1.0_dp]
        RKH_intensity = -1.0_dp
        RKH_omega = -1.0_dp
        n_quasi = -1
        gs_energy = (0.0_dp,0.0_dp)
    end subroutine set_RKH_defaults

    subroutine error_bad_input()
        write(stderr,*)
        write(stderr,*) "Some of the input parameters were bad, see errors above"
        stop
    end subroutine error_bad_input
end module input_tools