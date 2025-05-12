program mat_els_test
    use grid_tools
    use bspline_tools
    use block_tools
    use eig_tools
    use potentials
    use CAP_tools
    use mat_els
    use orbitals
    use hamiltonian
    use dipole
    use stdlib_linalg, only: eig,eigvals
    use stdlib_sorting, only: sort
    use omp_lib, only: omp_get_wtime
    implicit none

    integer :: k,m,Z,i,j,l,max_k,k_GL
    double precision :: h_max,r_max
    double precision, dimension(:), allocatable :: grid
    double complex, dimension(:,:), allocatable :: H,S,r_mat,dr_mat,vecs,ac,dip_block
    double complex, dimension(:), allocatable :: eigs,vec,vec_gs,vec_excited,vec_dip
    double precision, dimension(:), allocatable :: eigs_real
    type(sparse_4d) :: r_k,r_m_k
    type(sparse_6d) :: r_d_k
    type(sparse_Slater) :: Rk
    double precision :: off_diag_sum
    double complex :: test_int
    integer :: i_b,j_b,i_b_p,j_b_p,i_gs
    integer, dimension(1) :: index_gs

    double precision :: t_start,t_end,t_1_dip,t_2_dip

    double complex, dimension(:,:), allocatable :: eigs_v
    type(block), dimension(:), allocatable :: vecs_v
    type(block) :: H_block,S_block
    type(basis) :: bas
    type(block), dimension(:), allocatable :: H_vec
    type(config),dimension(2) :: confs
    type(sym), dimension(2) :: syms
    double precision, dimension(:,:,:,:,:),allocatable :: R_p
    double complex :: res
    integer, dimension(2) :: nnz
    type(CSC_matrix) H_spc,S_spc
    type(CSR_matrix) H_spr,S_spr,dip_block_spr
    type(block_diag_CS) :: H_diag,S_diag
    type(block_CS) :: dip,block_sum
    character :: gauge
    logical :: compute

    !for ploting orbitals
    double complex, dimension(:), allocatable :: vals
    double precision, dimension(:), allocatable :: r,c
    double precision :: dr

    type(b_spline) :: splines
    type(hydrogenic) :: pot

    type(CAP) :: CAP_c
    integer :: CAP_order = 2
    double precision :: CAP_r_0 = 20.d0
    double complex :: CAP_eta = dcmplx(1.d-3,0.d0)

    k = 4
    m = 3
    Z = 2
    h_max = 0.5d0
    r_max = 5.d0
    k_GL = splines%k + 6
    call generate_grid(k,m,Z,h_max,r_max,grid)
    call splines%init(k,grid)
    write(6,*) "n_b: ", splines%n_b
    write(6,*) splines%knots

    call CAP_c%init(CAP_order,CAP_r_0,CAP_eta)

    call pot%init(Z)
    l = 0

    allocate(H(splines%n-2,splines%n-2),S(splines%n-2,splines%n-2),vecs(splines%n-2,splines%n-2))
    allocate(r_mat(splines%n-2,splines%n-2),dr_mat(splines%n-2,splines%n-2))
    H = 0.d0
    S = 0.d0
    r_mat = 0.d0
    dr_mat = 0.d0
    vecs = 0.d0

    call setup_H_one_particle(pot,CAP_c,l,splines,k_GL,H)
    call setup_S(splines,k_GL,S)

    gauge = 'v'
    call setup_radial_dip(splines,k_GL,r_mat,gauge)

    max_k = 4
    ! allocate(r_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b))
    ! allocate(r_m_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b))
    ! allocate(r_d_k(0:max_k,size(splines%breakpoints)-1,splines%n_b,splines%n_b,splines%n_b,splines%n_b))
    !r_k = 0.d0
    !r_m_k = 0.d0
    !r_d_k = 0.d0
    t_start = omp_get_wtime()
    call r_k%init(max_k,splines)
    call r_m_k%init(max_k,splines)
    call r_d_k%init(max_k,splines)
    call setup_Slater_integrals(splines,max_k,k_GL,r_k,r_m_k,r_d_k)
    t_end = omp_get_wtime()
    write(6,*) "Time computing primitive cell integrals (s): ", t_end -t_start
    call Rk%init(max_k,splines,r_d_k,r_k,r_m_k)
    ! do i = 1,r_k%nnz
    !     do j = 1,r_k%nnz
    !         write(6,*) r_k%iv(i),r_k%i(i),r_k%j(i),r_k%iv(j),r_k%i(j),r_k%j(j)
    !     end do
    ! end do

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

    allocate(ac(splines%n_b,splines%n_b),source = dcmplx(0.d0))
    do i = 1,splines%n_b
        do j = 1,splines%n_b
            ac(j,i) = vecs(j,i_gs)*vecs(i,i_gs)
        end do
    end do

    t_start = omp_get_wtime()
    !test_int = Slater_k(splines%n_b,vecs(:,i_gs),vecs(:,i_gs),vecs(:,i_gs),vecs(:,i_gs),0,Rk)
    test_int = Slater_k(splines%n_b,ac,ac,0,Rk)
    t_end = omp_get_wtime()
    write(6,*) test_int,test_int/0.625d0
    write(6,*) "Slater integral time: ", t_end-t_start

    ! l = 1
    ! H = 0.d0
    ! S = 0.d0
    ! call setup_H_one_particle(pot,CAP_c,l,splines,k_GL,H)
    ! call setup_S(splines,k_GL,S)
    ! call eig_general(H,S,eigs,vecs)

    ! do i = 1,splines%n_b
    !     eigs_real(i) = abs(eigs(i)+0.375d0)
    ! end do

    ! index_gs = minloc(eigs_real)
    ! i_gs = index_gs(1)

    ! write(6,*) eigs(i_gs)

    ! !plot the 1s orbital
    ! allocate(r(1101),vals(1101),c(splines%n))
    ! c = 0.d0
    ! vals = 0.d0
    ! r = 0.d0
    ! dr = 0.01d0

    ! do i = 2,1101
    !     r(i) = r(i-1) + dr
    ! end do
    ! open(unit = 1, file = "1s.dat")
    ! do i = 1,1101
    !     do i_b = 1,splines%n_b
    !         c(i_b+1) = 1.d0
    !         vals(i) = vals(i) + vecs(i_b,i_gs)*splines%eval_d(r(i),c)
    !         c(i_b+1) = 0.d0
    !     end do
    !     write(1,*) r(i), -real(vals(i)), aimag(vals(i)), 2.d0*r(i)*exp(-r(i))
    ! end do
    ! close(1)

    allocate(eigs_v(splines%n_b,0:5),vecs_v(0:5))
    call find_orbitals(splines,k_GL,pot,CAP_c,5,eigs_v,vecs_v)

    ! Test single particle dipole
    allocate(vec(splines%n_b))
    vec = 0.d0
    res = 0.d0
    confs(1)%eqv = .false.
    confs(1)%l = 1
    confs(2)%l = 0

    do i = 1,splines%n_b
        confs(1)%n = i
        do j = 1,splines%n_b
            confs(2)%n = j
            vec(i) = vec(i) + vecs_v(0)%data(j,1)*dip_red_1p(confs,1,r_mat)
        end do
    end do
    vec = (-1)**1*three_j0(1,1,0)*vec
    res = dot_product(conjg(vecs_v(1)%data(:,1)),vec)
    write(6,*) "1s->2p z dipole: ", res

    S = 0.d0
    call setup_S(splines,k_GL,S)
    allocate(H_vec(0:3))
    do i =0,3
        call H_vec(i)%init(splines%n_b,splines%n_b,.true.)
        call setup_H_one_particle(pot,CAP_c,i,splines,k_GL,H_vec(i)%data)
    end do

    confs(1)%eqv = .false.
    confs(1)%l = [0,0]
    confs(1)%n = [splines%n_b,splines%n_b-2]
    confs(2)%eqv = .false.
    confs(2)%l = [0,0]
    confs(2)%n = [splines%n_b,splines%n_b-1]
    ! res = c_mat_eq(.true.,splines%n_b,0,confs,vecs_v(0)%data(:,1),vecs_v(0)%data(:,1),vecs_v(1)%data(:,1),vecs_v(1)%data(:,1),&
    !             max_k,Rk)
    ! write(6,*) res

    allocate(R_p(0:max_k,splines%n_b,splines%n_b,splines%n_b,splines%n_b),source=0.d0)

    ! do j = 0,max_k
    !     do i = 1,Rk%nnz
    !         R_p(j,Rk%i(i),Rk%j(i),Rk%i_p(i),Rk%j_p(i)) = R_p(j,Rk%i(i),Rk%j(i),Rk%i_p(i),Rk%j_p(i)) + Rk%data(i,j)
    !     end do
    ! end do

    call compute_R_k(r_d_k,r_k,r_m_k,max_k,splines%n_b,R_p)

    ! res = 0.d0
    ! do j_b_p = 1,splines%n_b
    !     do i_b_p = 1,splines%n_b
    !         do j = 1,splines%n_b
    !             do i = 1,splines%n_b
    !                 res = res + vecs_v(0)%data(i,1)*vecs_v(0)%data(j,1)*vecs_v(1)%data(i_b_p,1)&
    !                             *vecs_v(1)%data(j_b_p,1)*c_mat_eq_tens(.true.,0,confs,i,j,i_b_p,j_b_p,max_k,Rk,0,R_p)
    !             end do
    !         end do
    !     end do
    ! end do
    ! write(6,*) res

    res = c_mat_neq_tens(0,confs,splines%n_b,splines%n_b-2,splines%n_b,splines%n_b-1,max_k,Rk,0,R_p)
    write(6,*) res
    confs(1)%eqv = .false.
    confs(1)%l = [0,0]
    confs(1)%n = [splines%n_b,splines%n_b-1]
    confs(2)%eqv = .false.
    confs(2)%l = [0,0]
    confs(2)%n = [splines%n_b,splines%n_b-2]
    res = c_mat_neq_tens(0,confs,splines%n_b,splines%n_b-1,splines%n_b,splines%n_b-2,max_k,Rk,0,R_p)
    write(6,*) res
    !stop

    call bas%init(2,2,splines%n_b,.false.,eigs_v)
    ! call H_block%init(bas%syms(2)%n_config,bas%syms(2)%n_config,.true.)
    ! call S_block%init(bas%syms(2)%n_config,bas%syms(2)%n_config,.true.)
    write(6,*) bas%syms(1)%l,bas%syms(1)%m,bas%syms(1)%pi,bas%syms(1)%n_config
    nnz = count_nnz(splines,bas%syms(2),max_k,.true.)
    write(6,*) "nnz: ",nnz
    !call construct_block(H_block,bas%syms(1),eigs_v,vecs_v,max_k,Rk)
    ! call construct_block_tensor(H_block,H_vec,S_block,S,splines,bas%syms(2),max_k,Rk,R_p,H_spr,S_spr)
    !write(6,*) "Sparsity H: ", real(H_spr%nnz)/(H_spr%shape(1)**2)
    ! write(6,*) "Sparsity S: ", real(S_spr%nnz)/(S_spr%shape(1)**2)

    call H_diag%init(bas%n_sym)
    call S_diag%init(bas%n_sym)
    !$omp parallel do
    do i = 1,bas%n_sym
        call construct_block_tensor(H_vec,S,splines,bas%syms(i),max_k,Rk,R_p,H_diag%blocks(i),S_diag%blocks(i),.false.)
    end do
    !$omp end parallel do
    call H_diag%compute_shape()
    call S_diag%compute_shape()

    deallocate(R_p)

    ! write(6,*) H_spc%data(1),H_spr%data(1), "Hej!"
    ! write(6,*) H_block%data(1,1),S_block%data(1,1)
    ! write(6,*) H_vec(0)%data(1,1),S(1,1),2*H_vec(0)%data(1,1)*S(1,1),S(1,1)**2

    ! open(unit=1,file='H_test.dat')
    ! !H_block%data=abs(H_block%data-transpose(H_block%data))
    ! do i = 1,bas%syms(2)%n_config
    !     write(1,*) abs(H_block%data(i,:))
    ! end do
    ! close(1)
    ! open(unit=1,file='S_test.dat')
    ! do i = 1,bas%syms(2)%n_config
    !     write(1,*) abs(S_block%data(i,:))
    ! end do
    ! close(1)

    res = dcmplx(2.9d0,0.d0)
    ! call APX(H_diag,[res,res,res],S_diag)
    call FEAST(H_diag%blocks(1),S_diag%blocks(1),.true.,.false.,dcmplx(-2.5d0,0.0d0),.6d0,eigs,vecs,i)
    index_gs = minloc(real(eigs))
    i = index_gs(1)
    vec_gs = vecs(:,i)
    call FEAST(H_diag%blocks(3),S_diag%blocks(3),.true.,.false.,dcmplx(-2.5d0,0.0d0),.6d0,eigs,vecs,i)
    index_gs = minloc(real(eigs))
    i = index_gs(1)
    vec_excited = vecs(:,i)

    !call H_spr%store('H_csr.dat')
    !call H_spr%load('H_csr.dat')

    !write(6,*) "Check of symmetry, H: ", maxval(abs(H_block%data-transpose(H_block%data)))
    !write(6,*) "Check of symmetry, S: ", maxval(abs(S_block%data-transpose(S_block%data)))

    ! deallocate(H,S)
    ! call H_spr%get_dense(H)
    ! call S_spr%get_dense(S)
    ! write(6,*) maxval(abs(H-H_block%data)),maxval(abs(S-S_block%data)),H_spr%index_ptr(H_spr%shape(1)+1)
    call H_spr%deall()
    call S_spr%deall()
    call H_diag%store('H_diag.dat')
    call H_diag%deall()
    call H_diag%load('H_diag.dat')
    call H_diag%to_CS(H_spr,.false.)
    call S_diag%to_CS(S_spr,.false.)
    ! call FEAST(H_block%data,S_block%data,H_spr,S_spr,dcmplx(-2.5d0,0.0d0),0.6d0,eigs,vecs,i)

    ! deallocate(H)
    ! call H_spr%get_dense(H)
    ! open(unit=1,file='H_test.dat')
    ! !H_block%data=abs(H_block%data-transpose(H_block%data))
    ! do i = 1,H_spr%shape(1)
    !     write(1,*) abs(H(i,:))
    ! end do
    ! close(1)
    ! deallocate(H)
    ! call S_spr%get_dense(H)
    ! open(unit=1,file='S_test.dat')
    ! do i = 1,S_spr%shape(1)
    !     write(1,*) abs(H(i,:))
    ! end do
    ! close(1)

    !write(6,*) H_block%data(702,666),H_block%data(666,702)
    ! do i = 1,bas%syms(1)%n_config
    !     write(6,*) i,H_block%data(i,i)
    !     if (real(H_block%data(i,i))<0.d0) stop
    ! end do

    ! t_start = omp_get_wtime()
    ! eigs = eigvals(H_block%data)
    ! call sort_eigvals(size(eigs),eigs)
    ! t_end = omp_get_wtime()
    ! write(6,*) 'Time for diag and sort (s): ', t_end-t_start
    ! write(6,*) eigs(1),eigs(2)

    ! deallocate(eigs,vecs,H,S)
    ! allocate(eigs(bas%syms(1)%n_config),vecs(bas%syms(1)%n_config,bas%syms(1)%n_config))
    ! call H_diag%blocks(1)%get_dense(H)
    ! call S_diag%blocks(1)%get_dense(S)

    ! t_start = omp_get_wtime()
    ! call eig_general(H,S,eigs,vecs)
    ! call sort_eig(bas%syms(1)%n_config,eigs,vecs)
    ! t_end = omp_get_wtime()
    ! write(6,*) 'Time for diag and sort (s): ', t_end-t_start
    ! write(6,*) eigs(1),eigs(2)
    ! open(unit=1,file='eigs.dat')
    ! do i = 1,bas%syms(1)%n_config
    !     write(1,*) real(eigs(i)), aimag(eigs(i))
    ! end do
    ! close(1)
    ! deallocate(eigs,vecs)
    ! allocate(eigs(splines%n_b),vecs(splines%n_b,splines%n_b))
    ! call eig_general(H_vec(0)%data,S,eigs,vecs)
    ! call sort_eig(splines%n_b,eigs,vecs)
    ! write(6,*) eigs(1),eigs(2)

    allocate(dip_block(bas%syms(1)%n_config,bas%syms(3)%n_config))
    syms(1) = bas%syms(1)
    syms(2) = bas%syms(3)
    call construct_dip_block_tensor(syms,1,splines,S,r_mat,dip_block_spr,.true.)
    call construct_dip_block_dense(syms,1,S,r_mat,dip_block)

    allocate(vec_dip(size(vec_gs)))
    call zgemv('N',size(vec_gs),size(vec_excited),dcmplx(1.d0,0.d0),dip_block,size(vec_gs),vec_excited,1,&
                dcmplx(0.d0,0.d0),vec_dip,1)
    !call CSR_mv(S_diag%blocks(1),vec_excited,vec_dip)
    res = dot_product(conjg(vec_gs),vec_dip)
    !res = zdotu(size(vec_gs),vec_gs,size(vec_dip),vec_dip)
    write(6,*) "<1s^2 (0,0,e)|z|1s2p (1,0,o)>: ", res, abs(res),abs(res)/0.404d0

    call CSR_mv(dip_block_spr,vec_excited,vec_dip)
    !call CSR_mv(S_diag%blocks(1),vec_excited,vec_dip)
    res = dot_product(conjg(vec_gs),vec_dip)
    !res = zdotu(size(vec_gs),vec_gs,size(vec_dip),vec_dip)
    write(6,*) "<1s^2 (0,0,e)|z|1s2p (1,0,o)>: ", res, abs(res),abs(res)/0.404d0

    write(6,*) 'Constructing dipoles ...'
    t_1_dip = omp_get_wtime()
    call dip%init([bas%n_sym,bas%n_sym])
    !$omp parallel do private(i,syms,compute) schedule(dynamic)
    do j = 1,bas%n_sym
        syms(2) = bas%syms(j)
        do i = 1,bas%n_sym
            syms(1) = bas%syms(i)
            if (i>j) then
                compute = .false.
            else
                compute = .true.
            end if
            call construct_dip_block_tensor(syms,1,splines,S,r_mat,dip%blocks(i,j),compute)
        end do
    end do
    !$omp end parallel do
    call block_sum%init([bas%n_sym,bas%n_sym])
    !$omp parallel do private(i,syms,compute) schedule(dynamic)
    do j = 1,bas%n_sym
        syms(2) = bas%syms(j)
        do i = 1,bas%n_sym
            syms(1) = bas%syms(i)
            if (i>j) then
                compute = .false.
            else
                compute = .true.
            end if
            call construct_dip_block_tensor(syms,-1,splines,S,r_mat,block_sum%blocks(i,j),compute)
        end do
    end do
    !$omp end parallel do
    t_2_dip = omp_get_wtime()
    write(6,*) 'Done! Time for dipole construction (s): ', t_2_dip-t_1_dip

    dip = AXPBY(dip, (2.d0,0.d0), block_sum, (2.d0,0.d0))

    block_sum = XPAY(H_diag,dip,(1.d0,0.d0))

    ! call dip%store('dip_block.dat')
    ! call dip%deall()
    ! call dip%load('dip_block.dat')
    call dip_block_spr%deall()
    deallocate(dip_block)
    call block_sum%to_CS(dip_block_spr,.true.)
    call dip_block_spr%get_dense(dip_block)

    call FEAST(dip_block_spr,S_spr,.true.,.false.,dcmplx(-2.5d0,0.0d0),.6d0,eigs,vecs,i)

    ! res = maxval(abs(dip_block))
    ! dip_block=dip_block-transpose(dip_block)
    ! write(6,*) maxval(abs(dip_block)),maxval(abs(dip_block))/res
    ! Uncomment to write dipole block to file, useful for checking the sparsity pattern
    open(file='dip_test.dat',unit=1,action='write')
    do i = 1,ubound(dip_block,dim=1)
        write(1,*) abs(dip_block(i,:))
    end do
    close(1)
end program mat_els_test