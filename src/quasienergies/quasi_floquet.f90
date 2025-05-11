module quasi_floquet
    use kind_tools
    use constants_tools
    use sparse_array_tools
    use block_tools
    use orbital_tools
    use input_tools
    implicit none

contains

    subroutine setup_dip_floquet(dip,D,bas)
        type(block_CS), intent(inout) :: dip
        type(block_CS), intent(out) :: D
        type(basis), intent(in) :: bas

    end subroutine setup_dip_floquet

    subroutine setup_diag_blocks_floquet(bas,H_0,S)
        type(basis), intent(in) :: bas
        type(block_diag_CS), intent(in) ::  H_0
        type(block_diag_CS), intent(in) :: S

    end subroutine setup_diag_blocks_floquet

    subroutine setup_H_0_floquet(bas, H_0, S, omega, shift)
        type(basis), intent(in) :: bas
        type(block_diag_CS), intent(in) ::  H_0
        type(block_diag_CS), intent(in) :: S
        real(dp), intent(in) :: omega
        complex(dp), intent(in) :: shift

        double complex, allocatable :: shifts(:)

        ! shifts = omega*bas%m + shift
        ! call APX(H_0, shifts, S)
    end subroutine setup_H_0_floquet

end module quasi_floquet