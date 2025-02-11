program mat_els_test
    use grid_tools
    use bspline_tools
    use block_tools
    use eig_tools
    use potentials
    use CAP_tools
    use mat_els
    use orbitals
    use stdlib_linalg, only: eig,eigvals,linalg_state_type
    use stdlib_sorting, only: sort
    implicit none

    integer :: k,m,Z,i,j,l,max_k,k_GL
    double precision :: h_max,r_max
    double precision, dimension(:), allocatable :: grid
    double complex, dimension(:,:), allocatable :: H,S,vecs
    double complex, dimension(:), allocatable :: eigs
    double precision, dimension(:), allocatable :: eigs_real
    !double precision, dimension(:,:,:,:), allocatable :: r_k,r_m_k
    !double precision, dimension(:,:,:,:,:,:), allocatable :: r_d_k
    type(sparse_4d) :: r_k,r_m_k
    type(sparse_6d) :: r_d_k
    double precision :: diag_sum,off_diag_sum
    double complex :: test_int
    integer :: i_b,j_b,i_b_p,j_b_p,i_gs
    integer, dimension(1) :: index_gs
    type(linalg_state_type) :: err

    double complex, dimension(:,:), allocatable :: eigs_v
    type(block), dimension(:), allocatable :: vecs_v

    !for ploting orbitals
    double complex, dimension(:), allocatable :: vals
    double precision, dimension(:), allocatable :: r,c
    double precision :: dr

    type(b_spline) :: splines
    type(hydrogenic) :: pot

    type(CAP) :: CAP_c
    integer :: CAP_order = 2
    double precision :: CAP_r_0 = 60.d0
    double complex :: CAP_eta = dcmplx(1.d-3,0.d0)

    k = 8
    m = 3
    Z = 1
    h_max = 0.5d0
    r_max = 80.d0
    k_GL = splines%k + 6
    call generate_grid(k,m,Z,h_max,r_max,grid)
    call splines%init(k,grid)

    call CAP_c%init(CAP_order,CAP_r_0,CAP_eta)

    call pot%init(Z)
    l = 0

    allocate(H(splines%n-2,splines%n-2),S(splines%n-2,splines%n-2),vecs(splines%n-2,splines%n-2))
    H = 0.d0
    S = 0.d0
    vecs = 0.d0

    call setup_H_one_particle(pot,CAP_c,l,splines,k_GL,H)
    call setup_S(splines,k_GL,S)

    max_k = 0
    ! allocate(r_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b))
    ! allocate(r_m_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b))
    ! allocate(r_d_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b,splines%n_b,splines%n_b))
    !r_k = 0.d0
    !r_m_k = 0.d0
    !r_d_k = 0.d0
    call r_k%init(max_k,splines)
    call r_m_k%init(max_k,splines)
    call r_d_k%init(max_k,splines)
    call setup_Slater_integrals(splines,max_k,k_GL,r_k,r_m_k,r_d_k)

    allocate(eigs(splines%n_b),eigs_real(splines%n_b))

    ! open(unit=1,file='H_test.dat')
    ! do i = 1,splines%n-2
    !     write(1,*) real(H(i,:))
    ! end do
    ! close(1)
    ! open(unit=1,file='S_test.dat')
    ! do i = 1,splines%n-2
    !     write(1,*) real(S(i,:))
    ! end do
    ! close(1)

    ! eigs = eigvals(H,S)
    ! eigs_real = real(eigs)

    ! call sort(eigs_real)

    ! do i = 1,size(eigs_real)
    !     write(6,*) eigs_real(i)
    ! end do


    !call eig(H,S,eigs,right=vecs,err=err)
    call eig_general(H,S,eigs,vecs)
    call sort_eig(splines%n_b,eigs,vecs)

    open(unit=1,file="eigs.dat")
    do i = 1,splines%n_b
        eigs_real(i) = abs(eigs(i)+0.5d0)
        write(1,*) real(eigs(i)),aimag(eigs(i))
    end do
    close(1)

    index_gs = minloc(eigs_real)
    i_gs = index_gs(1)

    write(6,*) eigs(i_gs),i_gs


    test_int = 0.d0
    do i = 1,r_d_k%nnz
        test_int = test_int + vecs(r_d_k%i(i),i_gs)*vecs(r_d_k%j(i),i_gs)*&
                            vecs(r_d_k%i_p(i),i_gs)*vecs(r_d_k%j_p(i),i_gs)*r_d_k%data(i,0)

        ! test_int = test_int + vecs(r_d_k%j(i),i_gs)*vecs(r_d_k%i(i),i_gs)*&
        ! vecs(r_d_k%i_p(i),i_gs)*vecs(r_d_k%j_p(i),i_gs)*r_d_k%data(i)
    end do
    test_int = 2.d0*test_int

    do i = 1,r_k%nnz
        do j = 1,r_k%nnz
            off_diag_sum = 0.d0
            if (r_k%iv(i)<r_k%iv(j)) then
                off_diag_sum = r_k%data(i,0)*r_m_k%data(j,0)
            end if
            test_int = test_int + 2.d0*off_diag_sum*vecs(r_k%i(i),i_gs)*vecs(r_k%i(j),i_gs)*&
                        vecs(r_k%j(i),i_gs)*vecs(r_k%j(j),i_gs)
        end do
    end do
    write(6,*) test_int,test_int/0.625d0
    test_int = Slater_k(splines%n_b,vecs(:,i_gs),vecs(:,i_gs),vecs(:,i_gs),vecs(:,i_gs),0,r_d_k,r_k,r_m_k)
    write(6,*) test_int,test_int/0.625d0

    l = 1
    H = 0.d0
    S = 0.d0
    call setup_H_one_particle(pot,CAP_c,l,splines,k_GL,H)
    call setup_S(splines,k_GL,S)
    call eig_general(H,S,eigs,vecs)

    do i = 1,splines%n_b
        eigs_real(i) = abs(eigs(i)+0.375d0)
    end do

    index_gs = minloc(eigs_real)
    i_gs = index_gs(1)

    write(6,*) eigs(i_gs)

    !plot the 1s orbital
    allocate(r(1101),vals(1101),c(splines%n))
    c = 0.d0
    vals = 0.d0
    r = 0.d0
    dr = 0.01d0

    do i = 2,1101
        r(i) = r(i-1) + dr
    end do
    open(unit = 1, file = "1s.dat")
    do i = 1,1101
        do i_b = 1,splines%n_b
            c(i_b+1) = 1.d0
            vals(i) = vals(i) + vecs(i_b,i_gs)*splines%eval_d(r(i),c)
            c(i_b+1) = 0.d0
        end do
        write(1,*) r(i), -real(vals(i)), aimag(vals(i)), 2.d0*r(i)*exp(-r(i))
    end do
    close(1)

    allocate(eigs_v(splines%n_b,0:4),vecs_v(0:4))
    call find_orbitals(splines,k_GL,pot,CAP_c,4,eigs_v,vecs_v)

end program mat_els_test