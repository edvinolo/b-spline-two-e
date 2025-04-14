module orbitals
    use orbital_tools
    use bspline_tools
    use potentials
    use CAP_tools
    use block_tools
    use eig_tools
    use mat_els
    implicit none
contains
    subroutine find_orbitals(b_splines,k_GL,pot,CAP_c,max_l,eigs,vecs)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: k_GL
        class(sph_pot), intent(in) :: pot
        class(CAP), intent(in) :: CAP_c
        integer, intent(in) :: max_l
        double complex, dimension(b_splines%n_b,0:max_l), intent(out) :: eigs
        type(block), dimension(0:max_l), intent(out) :: vecs

        integer :: l
        type(block), dimension(0:max_l) :: H_vec
        double complex, dimension(:,:), allocatable :: S

        allocate(S(b_splines%n_b,b_splines%n_b))
        S = 0
        call setup_S(b_splines,k_GL,S)

        write(6,*) ''
        write(6,*) 'Diagonalizing one particle Hamiltonians to find orbitals. Maximum l = ', max_l
        do l = 0,max_l
            write(6,*) 'l = ', l
            call H_vec(l)%init(b_splines%n_b,b_splines%n_b,.true.)
            call vecs(l)%init(b_splines%n_b,b_splines%n_b,.true.)
            call setup_H_one_particle(pot,CAP_c,l,b_splines,k_GL,H_vec(l)%data)
            call eig_general(H_vec(l)%data,S,eigs(:,l),vecs(l)%data)
            call sort_eig(b_splines%n_b,eigs(:,l),vecs(l)%data)
            write(6,*) eigs(1,l)
        end do
        write(6,*) 'Done!'

    end subroutine find_orbitals
end module orbitals