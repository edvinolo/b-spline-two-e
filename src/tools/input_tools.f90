module input_tools
    use dir_tools, only: is_dir
    implicit none

    !Basis input variables
    character(len=:), allocatable :: basis_root_dir
    character(len=64) :: basis_output_dir

    integer :: k,m,Z,max_k,k_GL,max_L,max_l_1p
    double precision :: h_max,r_max

    integer :: CAP_order
    double precision :: CAP_r_0
    double complex :: CAP_eta

    character :: gauge
    logical :: z_pol

    namelist /basis_input/ &
    & basis_output_dir,&
    & k,&
    & m,&
    & Z,&
    & h_max,&
    & r_max,&
    & k_GL,&
    & CAP_order,&
    & CAP_r_0,&
    & CAP_eta,&
    & max_L,&
    & max_l_1p,&
    & max_k,&
    & gauge,&
    & z_pol
contains
    function get_input_file() result(res)
        character(len=:), allocatable :: res

        integer :: n_command
        character(len=256) :: input_file

        n_command = command_argument_count()
        if(n_command>1) then
            write(6,*)
            write(6,*) "More than 1 command line arguments supplied. Only the first one is used, the rest is ignored."
            write(6,*)
        else if (n_command == 0) then
            write(6,*)
            write(6,*) "No command line arguments supplied. You need to tell the program the name of your input file."
            write(6,*)
            stop
        end if

        call get_command_argument(1,input_file)
        res = input_file
    end function get_input_file

    subroutine get_basis_input(input_file)
        character(len=:), allocatable, intent(in) :: input_file

        open(unit=1,file=input_file,action='read')
        read(1,nml=basis_input)
        close(1)

        basis_root_dir = trim(basis_output_dir)
        if (.not.is_dir(basis_root_dir)) then
            write(6,*)
            write(6,*) "Error! The basis_output_dir you supplied: ", basis_root_dir, " does not exist."
            write(6,*) "You need to supply an exisiting directory to the basis_output_dir variable in the input file"
            write(6,*)
            stop
        end if

        if ((gauge /= 'l') .and. (gauge /= 'v')) then
            write(6,*)
            write(6,*) "Bad gauge input parameter: ", gauge
            write(6,*) "Must be l (length) or v (velocity)"
            write(6,*)
            stop
        end if

        write(6,nml=basis_input)
    end subroutine get_basis_input

    subroutine write_basis_input(basis_res_dir)
        character(len=:), allocatable, intent(in) :: basis_res_dir

        integer :: unit

        open(file = basis_res_dir // 'basis_input.dat', newunit = unit, action = 'write')
        write(unit, nml = basis_input)
        close(unit)
    end subroutine write_basis_input

    subroutine set_basis_defaults()
        k = 6
        m = 3
        Z = 2
        h_max = 0.5d0
        r_max = 15.d0
        k_GL = k + 6
        CAP_order = 2
        CAP_r_0 = 10.d0
        CAP_eta = dcmplx(1.d-3,0.d0)
        max_L = 2
        max_l_1p = 5
        max_k = 4
        z_pol = .true.
    end subroutine set_basis_defaults
end module input_tools