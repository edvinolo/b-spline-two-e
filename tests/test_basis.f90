program test_basis
    use orbital_tools
    implicit none

    integer :: max_L
    integer :: max_l_1p
    integer :: n_b
    logical :: z_pol
    double complex, dimension(:,:), allocatable :: eigs

    type(basis) :: bas_1,bas_2
    integer :: i,j
    logical :: configs_equal

    max_L = 2
    max_l_1p = 3
    n_b = 10
    z_pol = .false.
    allocate(eigs(n_b,n_b),source = dcmplx(0.d0,0.d0))

    call bas_1%init(max_L,max_l_1p,n_b,z_pol,eigs)

    do i = 1,bas_1%n_sym
        write(6,*) bas_1%syms(i)%l,bas_1%syms(i)%m,bas_1%syms(i)%pi,bas_1%syms(i)%n_config
    end do

    call bas_1%store('TESTOUTPUT/')
    call bas_2%load('TESTOUTPUT/')

    do i = 1,bas_1%n_sym
        configs_equal = .true.
        do j = 1,bas_1%syms(i)%n_config
            if (any(bas_1%syms(i)%configs(j)%n /= bas_2%syms(i)%configs(j)%n)) then
                configs_equal = .false.
                exit
            end if
            if (any(bas_1%syms(i)%configs(j)%l /= bas_2%syms(i)%configs(j)%l)) then
                configs_equal = .false.
                exit
            end if
            if (bas_1%syms(i)%configs(j)%eqv .neqv. bas_2%syms(i)%configs(j)%eqv) then
                configs_equal = .false.
                exit
            end if
        end do
        print *, bas_1%syms(i)%l == bas_2%syms(i)%l ,bas_1%syms(i)%m == bas_2%syms(i)%m ,bas_1%syms(i)%pi .eqv. bas_2%syms(i)%pi &
                ,bas_1%syms(i)%n_config == bas_2%syms(i)%n_config, configs_equal
    end do
end program test_basis