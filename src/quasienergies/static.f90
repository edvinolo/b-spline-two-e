module static
    use kind_tools
    use constants_tools
    use block_tools, only: block_CS,block_diag_CS,APX,XPAY
    use orbital_tools, only: basis
    implicit none

contains
    subroutine setup_H_static(H_block, bas, H_0, S, D, shift, V_0)
        type(block_CS), intent(inout) :: H_block
        type(basis), intent(in) :: bas
        type(block_diag_CS), intent(in) ::  H_0
        type(block_diag_CS), intent(in) :: S
        type(block_CS), intent(in) :: D
        complex(dp), intent(in) :: shift
        real(dp), intent(in) :: V_0

        complex(dp) :: shifts(bas%n_sym)

        H_block = XPAY(H_0, D, cmplx(V_0,kind=dp)) !Could add an XPAY specific function that accepts real scalars. Why no templates in Fortran :(
        call H_block%compute_shape()

        shifts = -shift
        call APX(H_block, shifts, S)
    end subroutine setup_H_static

end module static