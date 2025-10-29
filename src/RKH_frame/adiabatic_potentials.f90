module adiabatic_potentials
    use kind_tools
    use constants_tools
    use orbital_tools
    use wigner_tools
    use bspline_tools, only: b_spline
    use eig_tools, only: eig_real_sym
    use input_tools, only: N_r, r_limits
    use function_eval, only: array_func_1p, write_func_1p_from_vals,r_vec
    use stdlib_math, only: linspace
    implicit none

    integer :: n_sym
    integer :: max_k
    real(dp), allocatable :: U_r(:,:,:)
    real(dp), allocatable :: U_r_vecs(:,:,:)
    real(dp), allocatable :: U_r_vals(:,:)
    character(len = 1) :: uplo
    real(dp) :: r_0

    real(dp), allocatable :: ang_fact(:,:,:) ! Should a sparse data structure be used here?
contains

    function r_k(r,r_0_in,k) result(res)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: r_0_in
        integer, intent(in) :: k
        real(dp) :: res

        if (r>r_0) then
            res = r_0_in**k/r**(k+1)
        else
            res = r**k/r_0_in**(k+1)
        end if
    end function r_k

    subroutine setup_ang_fact(bas)
        type(basis), intent(in) :: bas

        integer i,j
        type(sym) :: syms(2)

        allocate(ang_fact(0:max_k,n_sym,n_sym), source = 0.0_dp)

        do j = 1,n_sym
            syms(2) = bas%syms(j)
            do i = 1,n_sym ! Store only upper triangle of matrix
                syms(1) = bas%syms(i)
                call compute_ang_fact(syms,i,j)
            end do
        end do
    end subroutine setup_ang_fact

    subroutine compute_ang_fact(syms,i,j)
        type(sym), intent(in) :: syms(2)
        integer, intent(in) :: i
        integer, intent(in) :: j

        integer :: k,q
        real(dp) :: k_fact,q_fact,Y_lm

        associate(l => syms(1)%l,&
                m => syms(1)%m,&
                l_p => syms(2)%l,&
                m_p => syms(2)%m)

            do k = abs(l-l_p),l+l_p,2 ! Triangle rule gives loop bounds, C-tensor l+k+l_p = even gives loop increment by two
                ! The conditions commented out below are equivalent to the above loop bounds and increment.
                ! if (mod(l+k+l_p,2) /= 0) cycle ! From C-tensor reduced matel
                ! if ((k<abs(l-l_p)).or.(k > l+l_p)) cycle ! Triangle rule
                k_fact = sqrt(4*pi/(2*k+1))*(-1)**k*C_red_mat(k,l,l_p)
                q_fact = (0.0_dp,0.0_dp)
                do q = -k,k,2 ! Only k + q = even are nonzero.
                    if (m /= q + m_p) cycle
                    Y_lm = real(conjg(sph_harm(k,q,pi/2,0.0_dp)),kind=dp) !The spherical harmonic is real along x-axis, so no actual need for conjg.
                    q_fact = q_fact + Y_lm*three_j(l,k,l_p,-m,q,m_p)*(-1)**(l-m)
                    ang_fact(k,i,j) = ang_fact(k,i,j) + q_fact
                end do 
                ang_fact(k,i,j) = ang_fact(k,i,j)*k_fact
            end do
        end associate

    end subroutine compute_ang_fact

    ! Matrix element of the operator u(r) from Tor Atomic Siegert states paper.
    ! The matrix element is evaluated between two Spherical harmonics, using the Laplace expansion for the shifted potential.
    function U_r_matel(i,j,syms,r,r_0_in,omega) result(res)
        integer, intent(in) :: i
        integer, intent(in) :: j
        type(sym), intent(in) :: syms(2)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: r_0_in
        real(dp), intent(in) :: omega
        real(dp) :: res

        integer :: k

        res = (0.0_dp,0.0_dp)

        associate(l => syms(1)%l,&
                m => syms(1)%m,&
                l_p => syms(2)%l,&
                m_p => syms(2)%m)

            if ((l == l_p).and.(m == m_p)) then
                res = l*(l+1)/(2*r**2) - m*omega
            end if

            do k = 0,l+l_p
                if (mod(l+k+l_p,2) /= 0) cycle ! From C-tensor reduced matel
                if ((k<abs(l-l_p)).or.(k > l+l_p)) cycle ! Triangle rule
                res = res - ang_fact(k,i,j)*r_k(r,r_0_in,k) ! Minus because the potential is attractive. 
            end do
        end associate
    end function U_r_matel

    subroutine setup_U_r(bas,intensity,omega)
        type(basis), intent(in) :: bas
        real(dp), intent(in) :: intensity
        real(dp), intent(in) :: omega

        type(sym) :: syms(2)
        integer :: i,j,i_r
        real(dp) :: E_0

        E_0 = E_0_au(intensity)/sqrt_2
        r_0 = E_0/(omega**2)
        print *,''
        print *, 'r_0: ', r_0
        print *,''

        if (.not.allocated(r_vec)) allocate(r_vec(N_r))
        r_vec = linspace(r_limits(1),r_limits(2),N_r)

        allocate(u_r(bas%n_sym,bas%n_sym,N_r),source = 0.0_dp)
        n_sym = bas%n_sym

        max_k = 2*bas%max_l_1p

        call setup_ang_fact(bas)

        uplo = 'U' ! Matrix is symmetric, build only upper triangular part
        do i_r = 1,N_r
            do j = 1,bas%n_sym
                syms(2) = bas%syms(j)
                do i = 1,j ! Matrix is symmetric, build only upper triangular part
                    syms(1) = bas%syms(i)
                    U_r(i,j,i_r) = U_r_matel(i,j,syms,r_vec(i_r),r_0,omega)
                end do
            end do
        end do
    end subroutine setup_U_r

    subroutine compute_adiabatic_potentials()
        integer :: i_r

        allocate(U_r_vecs(n_sym,n_sym,N_r),U_r_vals(n_sym,N_r))

        do i_r = 1,N_r
            call eig_real_sym(U_r(:,:,i_r),U_r_vals(:,i_r),U_r_vecs(:,:,i_r),uplo)
        end do
    end subroutine compute_adiabatic_potentials

    subroutine write_adiabatic_potentials(res_dir)
        character(len=*), intent(in) :: res_dir

        integer :: unit,i,j

        open(file = res_dir//'adiabatic_potentials.out',newunit = unit, action = 'write')
        write(unit,'(a)') '# Adiabatic potentials, first column is r, remaining are the ad. pots. in increasing order'
        do i = 1,N_r
            write(unit, '(es24.17e3)', advance = 'no') r_vec(i)
            do j = 1,n_sym
                write(unit, '(a,es25.17e3)', advance = 'no') ' ', U_r_vals(j,i)
            end do
            write(unit,*) ''
        end do
        close(unit)

    end subroutine write_adiabatic_potentials

    subroutine write_channel_functions(vecs,n_vecs,bas,b_splines,res_dir)
        complex(dp), intent(in) :: vecs(:,:)
        integer, intent(in) :: n_vecs
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: b_splines
        character(len=*), intent(in) :: res_dir

        complex(dp), allocatable :: vals_lm(:,:,:),vals_nu(:,:,:)
        integer :: i,j,k
        logical :: lm_basis
        character(len=:), allocatable :: res_file
        character(len=4) :: file

        allocate(vals_lm(N_r,n_sym,n_vecs),vals_nu(N_r,n_sym,n_vecs),source = (0.0_dp,0.0_dp))
        do i = 1,n_vecs
            do j = 1,n_sym
                call array_func_1p(vecs(:,i),bas,b_splines,vals_lm(:,:,i))
            end do
        end do

        do i = 1,n_vecs
            do j = 1,n_sym
                do k = 1,n_sym
                    vals_nu(:,j,i) = vals_nu(:,j,i) + U_r_vecs(k,j,:)*vals_lm(:,k,i)
                end do
            end do
        end do

        lm_basis = .true.
        do i = 1,n_vecs
            write(file,'(I4.4)') i
            res_file = res_dir//'rad_lm_'//file//'.out'
            call write_func_1p_from_vals(vals_lm(:,:,i),bas,res_file,lm_basis)
        end do

        lm_basis = .false.
        do i = 1,n_vecs
            write(file,'(I4.4)') i
            res_file = res_dir//'rad_nu_'//file//'.out'
            call write_func_1p_from_vals(vals_nu(:,:,i),bas,res_file,lm_basis)
        end do
    end subroutine write_channel_functions
end module adiabatic_potentials