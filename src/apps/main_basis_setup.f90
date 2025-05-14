program main_basis_setup
    use input_tools
    use dir_tools
    use grid_tools
    use bspline_tools
    use block_tools
    use orbital_tools
    use potentials
    use CAP_tools
    use mat_els
    use orbitals
    use hamiltonian
    use dipole
    use omp_lib, only: omp_get_wtime
    implicit none

    character(len=:), allocatable :: input_file
    character(len=:), allocatable :: basis_res_dir

    double precision, dimension(:), allocatable :: grid
    type(b_spline) :: splines
    integer :: max_n_b_2
    double complex, dimension(:,:), allocatable :: S
    double complex, dimension(:,:), allocatable :: radial_dip
    type(sparse_4d) :: r_k,r_m_k
    type(sparse_6d) :: r_d_k
    type(sparse_Slater) :: Rk
    double precision, dimension(:,:,:,:,:),allocatable :: R_p

    double complex, dimension(:,:), allocatable :: eigs_v
    type(block), dimension(:), allocatable :: vecs_v
    type(basis) :: bas
    type(block), dimension(:), allocatable :: H_vec
    type(block_diag_CS) :: H_diag,S_diag
    type(block_CS), dimension(-1:1) :: dipoles
    type(sym), dimension(2) :: syms

    type(hydrogenic) :: pot
    type(CAP) :: CAP_c

    integer :: i,j,q
    logical :: compute

    double precision :: t_start,t_end,t_prog_start,t_prog_end

    t_prog_start = omp_get_wtime()

    call set_basis_defaults()
    input_file = get_input_file()
    call get_basis_input(input_file)

    call generate_grid(k,m,Z,h_max,r_max,grid)
    call splines%init(k,grid)
    write(6,*) "n_b: ", splines%n_b
    write(6,*) splines%knots

    ! Limit the radial extent of the second electron if r_2_max > 0
    if (r_2_max > 0) then
        max_n_b_2 = splines%find_max_n_b(r_2_max)
    else
        max_n_b_2=splines%n_b
    end if

    call CAP_c%init(CAP_order,CAP_r_0,CAP_eta)

    call pot%init(Z)

    t_start = omp_get_wtime()
    call r_k%init(max_k,splines)
    call r_m_k%init(max_k,splines)
    call r_d_k%init(max_k,splines)
    call setup_Slater_integrals(splines,max_k,k_GL,r_k,r_m_k,r_d_k)
    t_end = omp_get_wtime()
    write(6,*) "Time computing primitive cell integrals (s): ", t_end -t_start
    call Rk%init(max_k,splines,r_d_k,r_k,r_m_k)

    allocate(R_p(0:max_k,splines%n_b,splines%n_b,splines%n_b,splines%n_b),source=0.d0)
    call compute_R_k(r_d_k,r_k,r_m_k,max_k,splines%n_b,R_p)

    allocate(eigs_v(splines%n_b,0:max_l_1p),vecs_v(0:max_l_1p))
    call find_orbitals(splines,k_GL,pot,CAP_c,max_l_1p,eigs_v,vecs_v)

    allocate(S(splines%n-2,splines%n-2))
    S = 0.d0
    call setup_S(splines,k_GL,S)

    allocate(H_vec(0:max_l_1p))
    do i =0,max_l_1p
        call H_vec(i)%init(splines%n_b,splines%n_b,.true.)
        call setup_H_one_particle(pot,CAP_c,i,splines,k_GL,H_vec(i)%data)
    end do

    call bas%init(max_L,max_l_1p,splines%n_b,max_n_b_2,z_pol,eigs_v)

    call H_diag%init(bas%n_sym)
    call S_diag%init(bas%n_sym)
    !$omp parallel do
    do i = 1,bas%n_sym
        call construct_block_tensor(H_vec,S,splines,bas%syms(i),max_k,Rk,R_p,H_diag%blocks(i),S_diag%blocks(i),full)
        write(6,'(a,i0,a,es0.4)') "Sparsity H(", i,"): ", real(H_diag%blocks(i)%nnz)/(H_diag%blocks(i)%shape(1)**2)
        write(6,'(a,i0,a,es0.4)') "Sparsity S(", i,"): ", real(S_diag%blocks(i)%nnz)/(S_diag%blocks(i)%shape(1)**2)
    end do
    !$omp end parallel do
    call H_diag%compute_shape()
    call S_diag%compute_shape()

    deallocate(R_p)

    ! Compute dipole matrix elements
    write(6,*)
    write(6,*) "Constructing dipoles..."
    allocate(radial_dip(splines%n_b,splines%n_b))
    call setup_radial_dip(splines,k_GL,radial_dip,gauge)
    t_start = omp_get_wtime()
    do q = -1,1
        call dipoles(q)%init([bas%n_sym,bas%n_sym])
        !$omp parallel do private(i,syms,compute) schedule(dynamic)
        do j = 1,bas%n_sym
            syms(2) = bas%syms(j)
            do i = 1,bas%n_sym
                if ((.not.full).and.(i>j)) then
                    compute = .false.
                else
                    compute = .true.
                end if
                syms(1) = bas%syms(i)
                call construct_dip_block_tensor(syms,q,splines,S,radial_dip,dipoles(q)%blocks(i,j),compute)
            end do
        end do
        !$omp end parallel do
    end do
    t_end = omp_get_wtime()
    write(6,*) "Done!"
    write(6,*) "Wall time for dipole construction (s): ", t_end-t_start
    write(6,*) "Wall time for dipole construction (m): ", (t_end-t_start)/60.d0
    write(6,*) "Wall time for dipole construction (h): ", (t_end-t_start)/3600.d0

    ! Store Matrix elements and basis information
    call make_res_dir(basis_root_dir,basis_res_dir)
    call H_diag%store(basis_res_dir//'H_diag.dat')
    call S_diag%store(basis_res_dir//'S_diag.dat')
    call dipoles(-1)%store(basis_res_dir//'D_-1.dat')
    call dipoles(0)%store(basis_res_dir//'D_0.dat')
    call dipoles(1)%store(basis_res_dir//'D_1.dat')
    call write_basis_input(basis_res_dir)
    call bas%store(basis_res_dir)

    t_prog_end = omp_get_wtime()
    write(6,*) "Total time for basis setup (s): ", t_prog_end - t_prog_start
    write(6,*) "Total time for basis setup (min): ", (t_prog_end - t_prog_start)/60.d0
    write(6,*) "Total time for basis setup (h): ", (t_prog_end - t_prog_start)/(60.d0*60.d0)
end program main_basis_setup