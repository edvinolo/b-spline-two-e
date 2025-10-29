module RKH_hamiltonian
    use kind_tools
    use constants_tools
    use bspline_tools
    use quad_tools, only: gau_leg
    use CAP_tools, only: CAP
    use orbital_tools, only: sym,basis
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS, block_diag_CS, APX
    use PARDISO_tools, only: PARDISO_solver
    use eig_tools, only: drive_ARPACK_SI
    use adiabatic_potentials, only: r_k,r_0,max_k,ang_fact,write_channel_functions
    implicit none

    real(dp), allocatable :: r_k_int(:,:,:)
    complex(dp), allocatable :: H_T_and_CAP(:,:,:)
    type(gau_leg) :: g_l

    complex(dp), allocatable :: eigs(:)
    complex(dp), allocatable :: vecs(:,:)

contains
    subroutine init_g_l(b_splines, k_GL)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: k_GL

        call g_l%init(b_splines%breakpoints,k_GL)
    end subroutine init_g_l

    subroutine setup_H_T_and_CAP(bas,CAP_c,b_splines)
        type(basis), intent(in) :: bas
        class(CAP), intent(in) :: CAP_c
        type(b_spline), intent(in) :: b_splines

        integer :: l,i_b,j_b,i_r
        real(dp), allocatable :: c_i(:),c_j(:)

        allocate(H_T_and_CAP(b_splines%n_b,b_splines%n_b,0:bas%max_l_1p),source = (0.0_dp,0.0_dp))
        allocate(c_i(b_splines%n),c_j(b_splines%n))

        do l = 0,bas%max_l_1p
            c_j = 0.d0

            do j_b = 1, b_splines%n-2
                c_j(j_b+1) = 1.d0
                c_j(j_b) = 0.d0

                c_i = 0.d0
                do i_b = 1, b_splines%n-2
                    c_i(i_b+1) = 1.d0
                    c_i(i_b) = 0.d0

                    if (abs(j_b-i_b)>=b_splines%k) then
                        cycle
                    end if

                    do i_r = 1, size(b_splines%breakpoints)-1
                        if (b_splines%support(i_r,i_b+1).and.b_splines%support(i_r,j_b+1)) then
                            call compute_H_T_and_CAP(CAP_c,l,b_splines,i_b,j_b,c_i,c_j,i_r)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine setup_H_T_and_CAP

    subroutine compute_H_T_and_CAP(CAP_c,l,b_splines,i,j,c_i,c_j,i_r)
        class(CAP), intent(in) :: CAP_c
        integer, intent(in) :: l
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        real(dp), intent(in) :: c_i(b_splines%n)
        real(dp), intent(in) :: c_j(b_splines%n)
        integer, intent(in) :: i_r


        integer :: i_sum
        double precision :: B_i, B_j, D_B_j

        do i_sum = 1, g_l%N
            B_i = b_splines%eval_d(g_l%x(i_sum,i_r),c_i,i_r)
            B_j = b_splines%eval_d(g_l%x(i_sum,i_r),c_j,i_r)
            D_B_j = b_splines%d_eval_d(g_l%x(i_sum,i_r),c_j,2,i_r)
            H_T_and_CAP(i,j,l) = H_T_and_CAP(i,j,l) + g_l%w(i_sum,i_r)*(-0.5*B_i*D_B_j + (l*(l+1)/(2*g_l%x(i_sum,i_r)**2)&
                                + CAP_c%V(g_l%x(i_sum,i_r)))*B_i*B_j)
        end do

    end subroutine compute_H_T_and_CAP

    subroutine setup_r_k_int(b_splines)
        type(b_spline), intent(in) :: b_splines

        integer :: i_b,j_b,i_r
        real(dp), allocatable :: c_i(:),c_j(:)

        allocate(r_k_int(0:max_k,b_splines%n_b,b_splines%n_b),source = 0.0_dp)
        allocate(c_i(b_splines%n),c_j(b_splines%n))
        c_j = 0.d0

        do j_b = 1, b_splines%n-2
            c_j(j_b+1) = 1.d0
            c_j(j_b) = 0.d0

            c_i = 0.d0
            do i_b = 1, b_splines%n-2
                c_i(i_b+1) = 1.d0
                c_i(i_b) = 0.d0

                if (abs(j_b-i_b)>=b_splines%k) then
                    cycle
                end if

                do i_r = 1, size(b_splines%breakpoints)-1
                    if (b_splines%support(i_r,i_b+1).and.b_splines%support(i_r,j_b+1)) then
                        r_k_int(:,i_b,j_b)  = r_k_int(:,i_b,j_b) + &
                            compute_r_k_int(b_splines,c_i,c_j,i_r)
                    end if
                end do
            end do
        end do

    end subroutine setup_r_k_int

    function compute_r_k_int(b_splines,c_i,c_j,i_r) result(res)
        type(b_spline), intent(in) :: b_splines
        real(dp), intent(in) :: c_i(b_splines%n)
        real(dp), intent(in) :: c_j(b_splines%n)
        integer, intent(in) :: i_r
        real(dp) :: res(0:max_k)

        integer :: k

        res = 0.0_dp
        do k = 0,max_k
            res(k) = res(k) + r_k_integral(k,b_splines,c_i,c_j,i_r)
        end do
    end function compute_r_k_int

    function r_k_integral(k,b_splines,c_i,c_j,i_r) result(res)
        integer, intent(in) :: k
        type(b_spline), intent(in) :: b_splines
        real(dp), intent(in) :: c_i(b_splines%n)
        real(dp), intent(in) :: c_j(b_splines%n)
        integer, intent(in) :: i_r
        real(dp) :: res


        integer :: i_sum
        double precision :: B_i, B_j

        res = 0.0_dp
        do i_sum = 1, g_l%N
            B_i = b_splines%eval_d(g_l%x(i_sum,i_r),c_i,i_r)
            B_j = b_splines%eval_d(g_l%x(i_sum,i_r),c_j,i_r)
            res = res + g_l%w(i_sum,i_r)*r_k(g_l%x(i_sum,i_r),r_0,k)*B_i*B_j
        end do

    end function r_k_integral

    function shift_Coul_matel(i,j,syms,i_sym,j_sym) result(res)
        integer, intent(in) :: i
        integer, intent(in) :: j
        type(sym), intent(in) :: syms(2)
        integer, intent(in) :: i_sym
        integer, intent(in) :: j_sym
        real(dp) :: res

        integer :: k

        res = 0.0_dp

        associate(l => syms(1)%l,&
                l_p => syms(2)%l)

            do k = abs(l-l_p),l+l_p,2 ! Loop index based on selection rules
                res = res  + ang_fact(k,i_sym,j_sym)*r_k_int(k,i,j)
            end do

        end associate

        res = -res ! Negative sign for attractive potential

    end function shift_Coul_matel

    subroutine shift_Coul_Ham_off_diag_block(syms,i_sym,j_sym,b_splines,H_block)
        type(sym), intent(in) :: syms(2)
        integer, intent(in) :: i_sym
        integer, intent(in) :: j_sym
        type(b_spline), intent(in) :: b_splines
        type(CSR_matrix), intent(out) :: H_block

        integer :: nnz,ptr
        integer :: i,j,n_i,n_j

        nnz = 0

        if (i_sym > j_sym) then
            call H_block%init([syms(1)%n_config,syms(2)%n_config],nnz)
            return
        end if

        do i = 1,syms(1)%n_config
            n_i = syms(1)%configs(i)%n(1)
            do j = 1,syms(2)%n_config
                n_j = syms(2)%configs(j)%n(1)
                if (abs(n_j-n_i)<b_splines%k) then
                    nnz = nnz + 1
                end if
            end do
        end do

        call H_block%init([syms(1)%n_config,syms(2)%n_config],nnz)

        ptr = 1
        H_block%index_ptr(1) = ptr
        do i = 1, syms(1)%n_config
            n_i = syms(1)%configs(i)%n(1)
            do j = 1,syms(2)%n_config
                n_j = syms(2)%configs(j)%n(1)
                if (abs(n_j-n_i) < b_splines%k) then
                    H_block%indices(ptr) = j
                    H_block%data(ptr) = shift_Coul_matel(n_j,n_i,syms,j_sym,i_sym) ! The matrix element is symmetric under j-i exchange
                    ptr = ptr + 1
                end if
            end do
            H_block%index_ptr(i+1) = ptr
        end do

    end subroutine shift_Coul_Ham_off_diag_block

    subroutine shift_Coul_Ham_diag_block(sym_diag,i_sym,b_splines,H_block)
        type(sym), intent(in) :: sym_diag
        integer, intent(in) :: i_sym
        type(b_spline), intent(in) :: b_splines
        type(CSR_matrix), intent(out) :: H_block

        integer :: nnz,ptr
        integer :: i,j,n_i,n_j
        type(sym) :: syms(2)

        syms = [sym_diag,sym_diag]

        nnz = 0
        do i = 1,sym_diag%n_config
            n_i = sym_diag%configs(i)%n(1)
            do j = i,sym_diag%n_config ! Only upper triangle
                n_j = sym_diag%configs(j)%n(1)
                if (abs(n_i-n_j)<b_splines%k) then
                    nnz = nnz + 1
                end if
            end do
        end do

        call H_block%init([sym_diag%n_config,sym_diag%n_config],nnz)


        ptr = 1
        H_block%index_ptr(1) = ptr
        do i = 1, sym_diag%n_config
            n_i = sym_diag%configs(i)%n(1)
            do j = i,sym_diag%n_config ! Only upper triangle
                n_j = sym_diag%configs(j)%n(1)
                if (abs(n_i-n_j) < b_splines%k) then
                    H_block%indices(ptr) = j
                    H_block%data(ptr) = shift_Coul_matel(n_j,n_i,syms,i_sym,i_sym) ! The matrix element is symmetric under j-i exchange
                    H_block%data(ptr) = H_block%data(ptr) + H_T_and_CAP(n_j,n_i,sym_diag%l) ! The matrix element is symmetric under j-i exchange
                    ptr = ptr + 1
                end if
            end do
            H_block%index_ptr(i+1) = ptr
        end do

    end subroutine shift_Coul_Ham_diag_block

    subroutine setup_shift_Coul_Ham(H_block,S_block,omega,shift,bas,CAP_c,k_GL,b_splines)
        type(block_CS), intent(out) :: H_block
        type(block_diag_CS), intent(in) :: S_block
        real(dp), intent(in) :: omega
        complex(dp), intent(in) :: shift
        type(basis), intent(in) :: bas
        type(CAP), intent(in) :: CAP_c
        integer, intent(in) :: k_GL
        type(b_spline), intent(in) :: b_splines

        integer :: i,j
        type(sym) :: syms(2)
        complex(dp) :: shifts(bas%n_sym)

        call init_g_l(b_splines,k_GL)
        call setup_H_T_and_CAP(bas,CAP_c,b_splines)
        call setup_r_k_int(b_splines)

        call H_block%init([bas%n_sym,bas%n_sym])

        ! Set diagonal blocks
        do i = 1,bas%n_sym
            call shift_Coul_Ham_diag_block(bas%syms(i),i,b_splines,H_block%blocks(i,i))
        end do

        ! Set off diagonal blocks
        do i = 1,bas%n_sym
            syms(1) = bas%syms(i)
            do j = 1,bas%n_sym
                if (i==j) cycle
                syms(2) = bas%syms(j)
                call shift_Coul_Ham_off_diag_block(syms,i,j,b_splines,H_block%blocks(i,j))
            end do
        end do

        ! Shift diagonal blocks
        shifts = -bas%syms%m*omega - shift
        ! shifts = -shift
        call APX(H_block,shifts,S_block)

    end subroutine setup_shift_Coul_Ham

    subroutine compute_RKH_energies(H_block,S_block,shift,bas,n_eigs,b_splines,res_dir)
        type(block_CS), intent(inout) :: H_block
        type(block_diag_CS), intent(inout) :: S_block
        complex(dp), intent(in) :: shift
        type(basis), intent(in) :: bas
        integer, intent(in) :: n_eigs
        type(b_spline), intent(in) :: b_splines
        character(len=*), intent(in) :: res_dir

        type(CSR_matrix) :: H,S
        type(PARDISO_solver) :: solver
        logical, parameter :: full = .false.

        call H_block%to_CS(H,.true.)
        call S_block%to_CS(S,.false.)

        allocate(eigs(n_eigs),vecs(bas%n_states,n_eigs))

        call solver%setup(H%shape(1),H%nnz,H%data,H%index_ptr,H%indices,full)
        call solver%factor(H%data,H%index_ptr,H%indices)
        call drive_ARPACK_SI(H, solver, S, full, shift, n_eigs, eigs, vecs)

        call write_channel_functions(vecs,n_eigs,bas,b_splines,res_dir)

        call write_eigs(res_dir,n_eigs)

    end subroutine compute_RKH_energies

    subroutine write_eigs(res_dir,n_eigs)
        character(len=*), intent(in) :: res_dir
        integer, intent(in) :: n_eigs

        ! character(len=200) :: fmt
        integer :: i,unit
        open(file = res_dir//"energies.out", newunit = unit, action = 'write')

        write(unit,'(a)') "# Energies in [a.u]"
        do i = 1,n_eigs
            write(unit,'(a)', advance='no') '('

            if (real(eigs(i),kind = dp)<0) then
                write(unit,'(es25.17e3)', advance='no') real(eigs(i),kind = dp)
            else
                write(unit,'(es24.17e3)', advance='no') real(eigs(i),kind = dp)
            end if

            if (aimag(eigs(i))<0) then
                write(unit,'(es25.17e3)', advance='no') aimag(eigs(i))
            else
                write(unit,'(a)', advance='no') '+'
                write(unit,'(es24.17e3)', advance='no') aimag(eigs(i))
            end if

            write(unit,'(a)') 'j) '
        end do

        close(unit)
    end subroutine write_eigs
end module RKH_hamiltonian