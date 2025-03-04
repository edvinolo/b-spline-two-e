program test_basis
    use orbital_tools
    implicit none

    integer :: max_L
    integer :: n_b
    type(basis) :: bas
    integer :: i

    max_L = 4
    n_b = 39

    call bas%init(max_L,n_b)

    do i = 1,bas%n_sym
        write(6,*) bas%syms(i)%l,bas%syms(i)%m,bas%syms(i)%pi,bas%syms(i)%n_config,&
                    log10(real(bas%syms(i)%n_config**2,kind=8)*sizeof(dcmplx(1.d0,1.d0)))
    end do

    write(6,*) sum(bas%syms%n_config),log10(real(sum(bas%syms%n_config)**2,kind=8)*sizeof(dcmplx(1.d0,1.d0)))
end program test_basis