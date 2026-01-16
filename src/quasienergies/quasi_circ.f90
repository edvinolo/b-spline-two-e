module quasi_circ
    use kind_tools
    use constants_tools
    use block_tools, only: block_CS,block_diag_CS,AXPBY,APX,XPAY
    use orbital_tools, only: basis
    implicit none

contains
    subroutine setup_dip_circ(dip,D,gauge)
        type(block_CS), intent(inout) :: dip(-1:1)
        type(block_CS), intent(out) :: D
        character, intent(in) :: gauge

        integer :: q

        if (gauge == 'l') then
            ! Use X because it is complex symmetric in length gauge
            ! X = -/sqrt(2)*(T_1 - T_-1)
            D = AXPBY(dip(1), -(inv_sqrt_2,0.0_dp), dip(-1), (inv_sqrt_2,0.0_dp))
        else if (gauge == 'v') then
            ! Use Y because it is complex symetric in velocity gauge
            ! Y = i/sqrt(2)*(T_1 + T_-1)
            D = AXPBY(dip(1), i_*inv_sqrt_2, dip(-1), i_*inv_sqrt_2)
        end if

        do q = -1,1
            deallocate(dip(q)%blocks)
        end do
    end subroutine setup_dip_circ

    subroutine setup_H_circ(H_block, bas, H_0, S, D, omega, shift, V_0, B_subset)
        type(block_CS), intent(inout) :: H_block
        type(basis), intent(in) :: bas
        type(block_diag_CS), intent(in) ::  H_0
        type(block_diag_CS), intent(in) :: S
        type(block_CS), intent(in) :: D
        real(dp), intent(in) :: omega
        complex(dp), intent(in) :: shift
        real(dp), intent(in) :: V_0
        logical, intent(in) :: B_subset

        complex(dp) :: shifts(bas%n_sym)

        H_block = XPAY(H_0, D, cmplx(V_0,kind=dp)) !Could add an XPAY specific function that accepts real scalars. Why no templates in Fortran :(
        call H_block%compute_shape()

        shifts = -omega*bas%syms%m - shift
        call APX(H_block, shifts, S, B_subset)
    end subroutine setup_H_circ
end module quasi_circ