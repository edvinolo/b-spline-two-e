module hamiltonian
    use block_tools
    use sparse_array_tools
    use orbital_tools
    use mat_els
    use omp_lib, only: omp_get_wtime
    implicit none

contains
    !Need to find a way to check which orbitals should be included in each symmetry
    !Only certain l combinations are allowed, and to reduce the basis size I will restrict the allowed range
    !for the sum of one-particle energies. Another option is to inlcude all allowed orbitals in each symmetry and diagonalize, and then remove the states
    !that have high energies. The latter option is probably better from physical standpoint.
    !The question is then how to do with the dipoles, since they mix the symmetries. The most straight forward approach is probably to compute them
    ! before diagonalizing, and then transforming to the diagonalized basis. But maybe this is extra work since I any way need to diagonalize the
    ! combined light-matter hamiltonian in the end. If I don't bother diagonalizing, then I have to use the one-particle energies to do the basis reduction.
    !For the basis I also have to consider parity, since it matters for the dipole selection rules, 
    !and I think it is preserved by the Coulomb interaction as well.
    !So for the Hamiltonian that means one block for each LMP combination, where L is total ang. mom., M is total magnetic quatum number, 
    !and P is parity (+1,-1). 

    !So it is probably best to just step through each symmetry and then count the allowed combinations

    subroutine construct_block(H_block,term,eigs,vecs,max_k,Rk)
        type(block), intent(inout) :: H_block
        type(sym), intent(in) :: term
        double complex, dimension(:,:), allocatable, intent(in) :: eigs
        type(block), dimension(:), allocatable, intent(in) :: vecs
        integer, intent(in) :: max_k
        type(sparse_Slater), intent(in) :: Rk

        logical :: both
        integer :: i,j,n
        integer :: n_a,n_b,n_c,n_d
        integer :: l_a,l_b,l_c,l_d
        type(config), dimension(2) :: confs

        double precision :: t_1,t_2

        n = size(eigs)

        write(6,*) term%l,term%m,term%pi,term%n_config

        t_1 = omp_get_wtime()
        if (mod(term%l,2)/=0) then
            do i = 1,term%n_config
                n_a = term%configs(i)%n(1)
                n_b = term%configs(i)%n(2)
                l_a = term%configs(i)%l(1)
                l_b = term%configs(i)%l(2)
                confs(1) = term%configs(i)

                H_block%data(i,i) = H_block%data(i,i) + eigs(n_a,l_a) + eigs(n_b,l_b)

                do j = 1,term%n_config
                    n_c = term%configs(j)%n(1)
                    n_d = term%configs(j)%n(2)
                    l_c = term%configs(j)%l(1)
                    l_d = term%configs(j)%l(2)
                    confs(2) = term%configs(j)

                    H_block%data(j,i) = H_block%data(j,i) + c_mat_neq(n,term%l,confs,&
                    vecs(l_a)%data(:,n_a),vecs(l_b)%data(:,n_b),vecs(l_c)%data(:,n_c),vecs(l_d)%data(:,n_d),&
                    max_k,Rk)
                end do
            end do
        else
            !$omp parallel do schedule(dynamic) shared(H_block) private(j,confs,both,n_a,n_b,l_a,l_b,n_c,n_d,l_c,l_d)
            do i = 1,term%n_config
                n_a = term%configs(i)%n(1)
                n_b = term%configs(i)%n(2)
                l_a = term%configs(i)%l(1)
                l_b = term%configs(i)%l(2)
                confs(1) = term%configs(i)

                !H_block%data(i,i) = H_block%data(i,i) + eigs(n_a,l_a) + eigs(n_b,l_b)

                do j = 1,term%n_config
                    both = (term%configs(i)%eqv.and.term%configs(j)%eqv)
                    n_c = term%configs(j)%n(1)
                    n_d = term%configs(j)%n(2)
                    l_c = term%configs(j)%l(1)
                    l_d = term%configs(j)%l(2)
                    confs(2) = term%configs(j)

                    if (term%configs(i)%eqv.or.term%configs(j)%eqv) then
                        H_block%data(j,i) = H_block%data(j,i) + c_mat_eq(both,n,term%l,confs,&
                        vecs(l_a)%data(:,n_a),vecs(l_b)%data(:,n_b),vecs(l_c)%data(:,n_c),vecs(l_d)%data(:,n_d),&
                        max_k,Rk)
                    else
                        H_block%data(j,i) = H_block%data(j,i) + c_mat_neq(n,term%l,confs,&
                        vecs(l_a)%data(:,n_a),vecs(l_b)%data(:,n_b),vecs(l_c)%data(:,n_c),vecs(l_d)%data(:,n_d),&
                        max_k,Rk)
                    end if
                    !H_block%data(i,j) = H_block%data(j,i)
                end do
            end do
            !$omp end parallel do
        end if

        t_2 = omp_get_wtime()
        write(6,*) 'Time to construct H_block (s): ', t_2-t_1
    end subroutine construct_block

    subroutine construct_block_tensor(H_block,H,S_Block,S,b_splines,term,max_k,Rk,R_k)
        type(block), intent(inout) :: H_block
        type(block), dimension(:), allocatable, intent(in) :: H
        type(block), intent(inout) :: S_block
        double complex, dimension(:,:), intent(in) :: S
        type(b_spline), intent(in) :: b_splines
        type(sym), intent(in) :: term
        integer, intent(in) :: max_k
        type(sparse_Slater), intent(in) :: Rk
        double precision, dimension(:,:,:,:,:), allocatable,intent(in) :: R_k

        logical :: both
        integer :: i,j,ptr
        integer :: n_a,n_b,n_c,n_d
        integer :: l_a,l_b,l_c,l_d
        type(config), dimension(2) :: confs
        !double precision, dimension(:,:,:,:,:), allocatable :: R_k

        double precision :: t_1,t_2

        !allocate(R_k(0:max_k,b_splines%n_b,b_splines%n_b,b_splines%n_b,b_splines%n_b),source=0.d0)

        ! do j = 0,max_k
        !     do i = 1,Rk%nnz
        !         R_k(j,Rk%i(i),Rk%j(i),Rk%i_p(i),Rk%j_p(i)) = R_k(j,Rk%i(i),Rk%j(i),Rk%i_p(i),Rk%j_p(i)) + Rk%data(i,j)
        !     end do
        ! end do


        write(6,*) term%l,term%m,term%pi,term%n_config

        write(6,*) R_k(0,1,1,2,1),R_k(0,2,1,1,1)
        write(6,*) R_k(0,1,2,3,4),R_k(0,1,4,3,2),R_k(0,3,2,1,4),R_k(0,3,4,1,2)

        t_1 = omp_get_wtime()
        ptr = 1
        if (mod(term%l,2)/=0) then
            do i = 1,term%n_config
                n_a = term%configs(i)%n(1)
                n_b = term%configs(i)%n(2)
                l_a = term%configs(i)%l(1)
                l_b = term%configs(i)%l(2)
                confs(1) = term%configs(i)

                do j = 1,term%n_config
                    n_c = term%configs(j)%n(1)
                    n_d = term%configs(j)%n(2)
                    l_c = term%configs(j)%l(1)
                    l_d = term%configs(j)%l(2)
                    confs(2) = term%configs(j)

                    H_block%data(j,i) = H_block%data(j,i) + c_mat_neq_tens(term%l,confs,&
                    n_a,n_b,n_c,n_d,max_k,Rk,0,R_k)
                    H_block%data(j,i) = H_block%data(j,i) + H_1p_neq(confs,term%l,H,S)
                    S_block%data(j,i) = S_block%data(j,i) + S_mat_neq(confs,term%l,S)
                end do
            end do
        else
            !!$omp parallel do schedule(dynamic) shared(H_block,S_block,H,S) private(j,confs,both,n_a,n_b,l_a,l_b,n_c,n_d,l_c,l_d)
            do i = 1,term%n_config
                n_a = term%configs(i)%n(1)
                n_b = term%configs(i)%n(2)
                l_a = term%configs(i)%l(1)
                l_b = term%configs(i)%l(2)
                confs(1) = term%configs(i)

                do j = 1,term%n_config
                    n_c = term%configs(j)%n(1)
                    n_d = term%configs(j)%n(2)
                    l_c = term%configs(j)%l(1)
                    l_d = term%configs(j)%l(2)
                    !if (abs(n_a-n_c)>=b_splines%k) cycle
                    !if (abs(n_b-n_d)>=b_splines%k) cycle
                    confs(2) = term%configs(j)

                    !write(6,*) n_a,n_b,n_c,n_d,l_a,l_b,l_c,l_d
                    if (term%configs(i)%eqv.or.term%configs(j)%eqv) then
                        both = (term%configs(i)%eqv.and.term%configs(j)%eqv)
                        H_block%data(j,i) = H_block%data(j,i) + c_mat_neq_tens(term%l,confs,&
                        n_a,n_b,n_c,n_d,max_k,Rk,0,R_k)
                        H_block%data(j,i) = H_block%data(j,i) + H_1p_neq(confs,term%l,H,S)
                        S_block%data(j,i) = S_block%data(j,i) + S_mat_neq(confs,term%l,S)
                    else
                        H_block%data(j,i) = H_block%data(j,i) + c_mat_neq_tens(term%l,confs,&
                        n_a,n_b,n_c,n_d,max_k,Rk,0,R_k)
                        H_block%data(j,i) = H_block%data(j,i) + H_1p_neq(confs,term%l,H,S)
                        S_block%data(j,i) = S_block%data(j,i) + S_mat_neq(confs,term%l,S)
                    end if
                    !H_block%data(i,j) = H_block%data(j,i)
                    !S_block%data(i,j) = S_block%data(j,i)
                    !write(6,*) n_a,n_b,n_c,n_d,l_a,l_b,l_c,l_d, H_block%data(i,j)
                end do
            end do
            !!$omp end parallel do
        end if

        t_2 = omp_get_wtime()
        write(6,*) 'Time to construct H_block (s): ', t_2-t_1
    end subroutine construct_block_tensor
end module hamiltonian