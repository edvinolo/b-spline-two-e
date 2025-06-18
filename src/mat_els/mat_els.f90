module mat_els
    use quad_tools
    use bspline_tools
    use sparse_array_tools
    use block_tools
    use potentials
    use CAP_tools
    use wigner_tools
    use orbital_tools
    implicit none

    procedure(dip_red_1p_fun), pointer :: dip_red_1p

    type, public :: radial_dipole
        character :: gauge
        double complex, allocatable :: r_mat(:,:)
        double complex, allocatable :: r_inv_mat(:,:)
        double complex, allocatable :: dr_mat(:,:)
    end type radial_dipole

    abstract interface
        subroutine rad_dip(b_splines,i,j,c_i,c_j,i_r,k_GL,r,w,radial_dip)
            import :: b_spline
            import :: radial_dipole
            type(b_spline), intent(in) :: b_splines
            integer, intent(in) :: i
            integer, intent(in) :: j
            double precision, dimension(b_splines%n), intent(in) :: c_i
            double precision, dimension(b_splines%n), intent(in) :: c_j
            integer, intent(in) :: i_r
            integer, intent(in) :: k_GL
            double precision, dimension(k_GL), intent(in) :: r
            double precision, dimension(k_GL), intent(in) :: w
            type(radial_dipole), intent(inout) :: radial_dip
        end subroutine rad_dip

        pure function dip_red_1p_fun(confs,i,radial_dip) result(res)
            import :: config
            import :: radial_dipole
            type(config), dimension(2), intent(in) :: confs
            integer, intent(in) :: i
            type(radial_dipole), intent(in) :: radial_dip
            double complex :: res
        end function dip_red_1p_fun
    end interface
contains
    subroutine setup_H_one_particle(potential,CAP_c,l,b_splines,k_GL,H)
        class(sph_pot), intent(in) :: potential
        class(CAP), intent(in) :: CAP_c
        integer, intent(in) :: l
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: k_GL
        double complex, dimension(:,:), intent(out) :: H

        integer :: i_b,j_b,i_r
        double precision, dimension(:), allocatable :: c_i,c_j
        type(gau_leg) :: g_l

        call g_l%init(b_splines%breakpoints,k_GL)
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
                        call compute_H(potential,CAP_c,l,b_splines,i_b,j_b,c_i,c_j,i_r,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),H)
                    end if
                end do
            end do
        end do
    end subroutine setup_H_one_particle

    subroutine setup_S(b_splines,k_GL,S)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: k_GL
        double complex, dimension(:,:), intent(out) ::  S

        integer :: i_b,j_b,i_r
        double precision, dimension(:), allocatable :: c_i,c_j
        type(gau_leg) :: g_l

        call g_l%init(b_splines%breakpoints,k_GL)
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
                        call compute_S(b_splines,i_b,j_b,c_i,c_j,i_r,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),S)
                    end if
                end do
            end do
        end do
    end subroutine setup_S

    subroutine setup_radial_dip(b_splines,k_GL,radial_dip,gauge)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: k_GL
        type(radial_dipole), intent(out) :: radial_dip
        character, intent(in) :: gauge

        integer :: i_b,j_b,i_r
        double precision, dimension(:), allocatable :: c_i,c_j
        type(gau_leg) :: g_l

        procedure(rad_dip), pointer :: compute_radial_dip

        if (gauge == 'l') then
            compute_radial_dip => compute_radial_dip_len
            dip_red_1p => dip_red_1p_len
            allocate(radial_dip%r_mat(b_splines%n_b,b_splines%n_b),source=dcmplx(0.d0,0.d0))
        else if (gauge == 'v') then
            compute_radial_dip => compute_radial_dip_vel
            dip_red_1p => dip_red_1p_vel
            allocate(radial_dip%r_inv_mat(b_splines%n_b,b_splines%n_b),source=dcmplx(0.d0,0.d0))
            allocate(radial_dip%dr_mat(b_splines%n_b,b_splines%n_b),source=dcmplx(0.d0,0.d0))
        else
            write(6,*) "Error! Unrecognized gauge option in for radial dipole: ", gauge
            stop
        end if

        call g_l%init(b_splines%breakpoints,k_GL)
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
                        call compute_radial_dip(b_splines,i_b,j_b,c_i,c_j,i_r,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),radial_dip)
                    end if
                end do
            end do
        end do
    end subroutine

    subroutine setup_Slater_integrals(b_splines,max_k,k_GL,r_k,r_m_k,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        type(sparse_4d), intent(inout) :: r_k
        type(sparse_4d), intent(inout) :: r_m_k
        type(sparse_6d), intent(inout) :: r_d_k

        call setup_Slater_off_diag(b_splines,max_k,k_GL,r_k,r_m_k)
        call setup_Slater_diag(b_splines,max_k,k_GL,r_d_k)
    end subroutine setup_Slater_integrals

    subroutine setup_Slater_off_diag(b_splines,max_k,k_GL,r_k,r_m_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        type(sparse_4d), intent(inout) :: r_k
        type(sparse_4d), intent(inout) :: r_m_k

        integer :: i_b,j_b,i_r,k,ptr
        double precision, dimension(k_GL) :: x,w
        double precision, dimension(2) :: limits
        double precision, dimension(:), allocatable :: c_i,c_j

        call setup_GL(k_GL,-1.d0,1.d0,x,w)

        allocate(c_i(b_splines%n),c_j(b_splines%n))
        !$omp parallel do private(ptr,c_i,c_j,i_b,j_b,i_r,limits)
        do k = 0,max_k
            ptr = 0
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
                                ptr = ptr + 1
                                limits(1) = b_splines%breakpoints(i_r)
                                limits(2) = b_splines%breakpoints(i_r+1)
                                call compute_Slater_off_diag(b_splines,i_r,i_b,j_b,c_i,c_j,k_GL,x,w,limits,k,ptr,r_k,r_m_k)
                            end if
                        end do
                    end do
            end do
        end do
        !$omp end parallel do
    end subroutine setup_Slater_off_diag

    subroutine setup_Slater_diag(b_splines,max_k,k_GL,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        type(sparse_6d), intent(inout) :: r_d_k

        integer :: i_b,j_b,i_b_p,j_b_p,i_r,k,ptr
        integer, dimension(4) :: i
        double precision, dimension(k_GL) :: x,w
        double precision, dimension(2) :: limits
        double precision, dimension(:,:), allocatable :: c
        logical :: j_support, i_support

        call setup_GL(k_GL,-1.d0,1.d0,x,w)

        allocate(c(b_splines%n,4))
        i = 0
        !$omp parallel do private(ptr,c,i,i_b,j_b,i_b_p,j_b_p,i_r,i_support,j_support,limits)
        do k = 0,max_k
            ptr = 0
            c = 0.d0
            do j_b_p = 1, b_splines%n-2
                c(:,3) = 0.d0
                c(j_b_p+1,4) = 1.d0
                c(j_b_p,4) = 0.d0
                i(4) = j_b_p
                do j_b = 1, b_splines%n-2
                    if (abs(j_b-j_b_p)>=b_splines%k) then
                        cycle
                    end if
                    c(:,2) = 0.d0
                    c(j_b+1,3) = 1.d0
                    c(j_b,3) = 0.d0
                    i(3) = j_b
                    do i_b_p = 1, b_splines%n-2
                        c(i_b_p+1,2) = 1.d0
                        c(i_b_p,2) = 0.d0
                        i(2) = i_b_p
                        c(:,1) = 0.d0
                        do i_b = 1, b_splines%n-2
                            c(i_b+1,1) = 1.d0
                            c(i_b,1) = 0.d0
                            i(1) = i_b
                            if (abs(i_b-i_b_p)>=b_splines%k) then
                                cycle
                            end if
                                do i_r = 1, size(b_splines%breakpoints)-1
                                    i_support = (b_splines%support(i_r,i_b+1).and.b_splines%support(i_r,i_b_p+1))
                                    j_support = (b_splines%support(i_r,j_b+1).and.b_splines%support(i_r,j_b_p+1))
                                    if (i_support.and.j_support) then
                                        ptr = ptr + 1
                                        limits(1) = b_splines%breakpoints(i_r)
                                        limits(2) = b_splines%breakpoints(i_r+1)
                                        call compute_Slater_diag(b_splines,i_r,i,c,k_GL,x,w,limits,k,ptr,r_d_k)
                                    end if
                                end do
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine setup_Slater_diag

    subroutine compute_H(potential,CAP_c,l,b_splines,i,j,c_i,c_j,i_r,k_GL,r,w,H)
        ! This is way too many arguments. Need to think about how to do this smarter.
        class(sph_pot), intent(in) :: potential
        class(CAP), intent(in) :: CAP_c
        integer, intent(in) :: l
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: i_r
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        double complex, dimension(:,:), intent(out) :: H


        integer :: i_sum
        double precision :: B_i, B_j, D_B_j

        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i,i_r)
            B_j = b_splines%eval_d(r(i_sum),c_j,i_r)
            D_B_j = b_splines%d_eval_d(r(i_sum),c_j,2,i_r)
            H(i,j) = H(i,j) + w(i_sum)*(-0.5*B_i*D_B_j + (potential%V(r(i_sum),l) + CAP_c%V(r(i_sum)))*B_i*B_j)
        end do

    end subroutine compute_H

    subroutine compute_S(b_splines,i,j,c_i,c_j,i_r,k_GL,r,w,S)
        ! This is way too many arguments. Need to think about how to do this smarter.
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: i_r
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        double complex, dimension(:,:), intent(out) :: S


        integer :: i_sum
        double precision :: B_i, B_j

        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i,i_r)
            B_j = b_splines%eval_d(r(i_sum),c_j,i_r)
            S(i,j) = S(i,j) + w(i_sum)*B_i*B_j
        end do

    end subroutine compute_S

    subroutine compute_radial_dip_len(b_splines,i,j,c_i,c_j,i_r,k_GL,r,w,radial_dip)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: i_r
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        type(radial_dipole), intent(inout) :: radial_dip

        integer :: i_sum
        double precision :: B_i,B_j
        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i,i_r)
            B_j = b_splines%eval_d(r(i_sum),c_j,i_r)
            radial_dip%r_mat(i,j) = radial_dip%r_mat(i,j) + w(i_sum)*r(i_sum)*B_i*B_j
        end do
    end subroutine compute_radial_dip_len

    subroutine compute_radial_dip_vel(b_splines,i,j,c_i,c_j,i_r,k_GL,r,w,radial_dip)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: i_r
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        type(radial_dipole), intent(inout) :: radial_dip

        integer :: i_sum
        double precision :: B_i,B_j,D_B_j
        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i,i_r)
            B_j = b_splines%eval_d(r(i_sum),c_j,i_r)
            D_B_j = b_splines%d_eval_d(r(i_sum),c_j,1,i_r)
            radial_dip%dr_mat(i,j) = radial_dip%dr_mat(i,j) + dcmplx(0.d0,-1.d0)*w(i_sum)*B_i*D_B_j
            radial_dip%r_inv_mat(i,j) = radial_dip%r_inv_mat(i,j) + dcmplx(0.d0,-1.d0)*w(i_sum)*B_i*B_j/r(i_sum)
        end do
    end subroutine compute_radial_dip_vel

    subroutine compute_Slater_off_diag(b_splines,i_r,i,j,c_i,c_j,k_GL,x,w,limits,k,ptr,r_k,r_m_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i_r
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: x
        double precision, dimension(k_GL), intent(in) :: w
        double precision, dimension(2), intent(in) :: limits
        integer, intent(in) :: k
        integer, intent(in) :: ptr
        type(sparse_4d), intent(inout) :: r_k
        type(sparse_4d), intent(inout) :: r_m_k

        integer :: i_sum
        double precision :: scale,translate
        double precision, dimension(k_GL) :: r
        double precision, dimension(k_GL) :: B_i,B_j

        scale = 0.5d0*(limits(2)-limits(1))
        translate = 0.5d0*(limits(2)+limits(1))
        do i_sum = 1, k_GL
            r(i_sum) = scale*x(i_sum) + translate
            B_i(i_sum) = b_splines%eval_d(r(i_sum),c_i,i_r)
            B_j(i_sum) = b_splines%eval_d(r(i_sum),c_j,i_r)
        end do

        do i_sum = 1,k_GL
            r_k%data(ptr,k) = r_k%data(ptr,k) + w(i_sum)*B_i(i_sum)*B_j(i_sum)*r(i_sum)**k
            r_m_k%data(ptr,k) = r_m_k%data(ptr,k) + w(i_sum)*B_i(i_sum)*B_j(i_sum)/r(i_sum)**(k+1)
        end do

        !Apply jacobian
        r_k%data(ptr,k) = scale*r_k%data(ptr,k)
        r_m_k%data(ptr,k) = scale*r_m_k%data(ptr,k)

        if (k==0) then
            r_k%iv(ptr) = i_r
            r_k%i(ptr) = i
            r_k%j(ptr) = j

            r_m_k%iv(ptr) = i_r
            r_m_k%i(ptr) = i
            r_m_k%j(ptr) = j
        end if
    end subroutine compute_Slater_off_diag

    subroutine compute_Slater_diag(b_splines,i_r,i,c,k_GL,x,w,limits,k,ptr,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i_r
        integer, dimension(4), intent(in) :: i
        double precision, dimension(b_splines%n,4), intent(in) :: c
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: x
        double precision, dimension(k_GL), intent(in) :: w
        double precision, dimension(2), intent(in) :: limits
        integer, intent(in) :: k
        integer, intent(in) :: ptr
        type(sparse_6d), intent(inout) :: r_d_k

        integer :: i_sum,j_sum
        double precision :: scale_i,translate_i,scale_j,translate_j
        double precision, dimension(k_GL) :: r
        double precision, dimension(k_GL) :: B_i,B_i_p
        double precision :: int_j,B_j,B_j_p,r_j

        scale_i = 0.5d0*(limits(2)-limits(1))
        translate_i= 0.5d0*(limits(2)+limits(1))
        do i_sum = 1, k_GL
            r(i_sum) = scale_i*x(i_sum) + translate_i
            B_i(i_sum) = b_splines%eval_d(r(i_sum),c(:,1),i_r)
            B_i_p(i_sum) = b_splines%eval_d(r(i_sum),c(:,2),i_r)
        end do

        do i_sum = 1,k_GL
            scale_j = 0.5d0*(r(i_sum)-limits(1))
            translate_j = 0.5d0*(r(i_sum)+limits(1))
            int_j = 0.d0
            do j_sum = 1,k_GL
                r_j = scale_j*x(j_sum) + translate_j
                B_j = b_splines%eval_d(r_j,c(:,3),i_r)
                B_j_p = b_splines%eval_d(r_j,c(:,4),i_r)
                int_j = int_j + w(j_sum)*B_j*B_j_p*r_j**k
            end do
            int_j = scale_j*int_j
            r_d_k%data(ptr,k) = r_d_k%data(ptr,k) &
                                                + w(i_sum)*B_i(i_sum)*B_i_p(i_sum)*int_j/r(i_sum)**(k+1)
        end do

        r_d_k%data(ptr,k) = scale_i*r_d_k%data(ptr,k)
        if (k==0) then
            r_d_k%iv(ptr) = i_r
            r_d_k%i(ptr) = i(1)
            r_d_k%i_p(ptr) = i(2)
            r_d_k%j(ptr) = i(3)
            r_d_k%j_p(ptr) = i(4)
        end if
    end subroutine compute_Slater_diag

    function Slater_k(n_b,ac,bd,k,Rk) result(res)
        integer, intent(in) :: n_b
        !double complex, dimension(n_b,n_b), intent(in) :: a,b,c,d
        double complex, dimension(n_b,n_b), intent(in) :: ac,bd
        integer, intent(in) :: k
        type(sparse_Slater), intent(in) :: Rk
        double complex :: res

        double complex :: bcd,temp
        integer i,j,j_p

        res = 0.d0
        j = Rk%j(1)
        j_p = Rk%j_p(1)
        temp = 0.d0
        bcd = bd(Rk%j(1),Rk%j_p(1))
        do i = 1, Rk%nnz
            if (j/=Rk%j(i)) then
                res = res + bd(j,j_p)*temp
                temp = 0.d0
                j = Rk%j(i)
                j_p = Rk%j_p(i)
            end if
            temp = temp + ac(Rk%i(i),Rk%i_p(i))*(Rk%data(i,k))
        end do
        res = res + temp*bd(Rk%j(Rk%nnz),Rk%j_p(Rk%nnz))
    end function Slater_k

    function r_12(n_b,L,configs,a,b,c,d,max_k,Rk) result(res)
        integer, intent(in) :: n_b
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        double complex, dimension(n_b), intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(sparse_Slater), intent(in) :: Rk
        double complex :: res

        integer :: k,i,j
        double complex, dimension(:,:), allocatable :: ac,bd
        double precision :: ang

        allocate(ac(n_b,n_b),bd(n_b,n_b),source = dcmplx(0.d0,0.d0))

        do i = 1,n_b
            do j = 1,n_b
                ac(j,i) = a(j)*c(i)
                bd(j,i) = b(j)*d(i)
            end do
        end do


        res = 0.d0
        do k = 0, max_k
            ang = ang_k_LS(k,configs,L)
            if (abs(ang)<5.d-16) cycle
            res = res + Slater_k(n_b,ac,bd,k,Rk)*ang
        end do
    end function r_12

    function r_12_tens(L,configs,a,b,c,d,max_k,R_k) result(res)
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        integer, intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(Nd_DOK), intent(inout) :: R_k
        double precision :: res

        integer :: k
        double precision :: ang
        double precision :: vals(0:max_k)

        res = 0.d0
        call R_k%get_val([a,b,c,d],vals)
        do k = 0, max_k
            ang = ang_k_LS(k,configs,L)
            if (abs(ang)<5.d-16) cycle
            res = res + vals(k)*ang
        end do
    end function r_12_tens

    function c_mat_neq(n_b,L,configs,a,b,c,d,max_k,Rk) result(res)
        integer, intent(in) :: n_b
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        double complex, dimension(n_b), intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(sparse_Slater), intent(in) :: Rk
        double complex :: res

        type(config), dimension(2) :: configs_ex

        configs_ex(1) = configs(1)
        configs_ex(2)%n = configs(2)%n([2,1])
        configs_ex(2)%l = configs(2)%l([2,1])
        res = r_12(n_b,L,configs,a,b,c,d,max_k,Rk)
        res = res + (-1)**(configs(2)%l(1) + configs(2)%l(2) + L)&
                *r_12(n_b,L,configs_ex,a,b,d,c,max_k,Rk)

    end function c_mat_neq

    function c_mat_eq(both,n_b,L,configs,a,b,c,d,max_k,Rk) result(res)
        logical, intent(in) :: both
        integer, intent(in) :: n_b
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        double complex, dimension(n_b), intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(sparse_Slater), intent(in) :: Rk
        double complex :: res

        res = r_12(n_b,L,configs,a,b,c,d,max_k,Rk)

        if (.not.both) res = sqrt(2.d0)*res
    end function c_mat_eq

    function c_mat_neq_tens(L,configs,a,b,c,d,max_k,R_k,support,support_ex) result(res)
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        integer, intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(Nd_DOK), intent(inout) :: R_k
        logical, intent(in) :: support
        logical, intent(in) :: support_ex
        double precision :: res

        type(config), dimension(2) :: configs_ex

        configs_ex(1) = configs(1)
        configs_ex(2)%n = configs(2)%n([2,1])
        configs_ex(2)%l = configs(2)%l([2,1])
        res = 0.0d0

        if (support) then
            res = res + r_12_tens(L,configs,a,b,c,d,max_k,R_k)
        end if

        if (support_ex) then
            res = res + (-1)**(configs(2)%l(1) + configs(2)%l(2) + L)&
                    *r_12_tens(L,configs_ex,a,b,d,c,max_k,R_k)
        end if
    end function c_mat_neq_tens

    function c_mat_eq_tens(both,L,configs,a,b,c,d,max_k,R_k) result(res)
        logical, intent(in) :: both
        integer, intent(in) :: L
        type(config), dimension(2), intent(in) :: configs
        integer, intent(in) :: a,b,c,d
        integer, intent(in) :: max_k
        type(Nd_DOK), intent(inout) :: R_k
        double precision :: res

        res = r_12_tens(L,configs,a,b,c,d,max_k,R_k)

        if (.not.both) res = sqrt(2.d0)*res
    end function c_mat_eq_tens

    pure function S_mat_eq(both,confs,S) result(res)
        logical, intent(in) :: both
        type(config), dimension(2), intent(in) :: confs
        double complex, dimension(:,:), intent(in) :: S
        double complex :: res

        if ((confs(1)%l(1)/=confs(2)%l(1)).or.(confs(1)%l(2)/=confs(2)%l(2))) then
            res = 0.d0
            return
        end if

        res = S(confs(1)%n(1),confs(2)%n(1))*S(confs(1)%n(2),confs(2)%n(2))
        if (.not.both) res = sqrt(2.d0)*res
    end function S_mat_eq

    pure function S_mat_neq(confs,L,S) result(res)
        type(config), dimension(2), intent(in) :: confs
        integer, intent(in) :: L
        double complex, dimension(:,:), intent(in) :: S
        double complex :: res

        res = 0.d0
        if ((confs(1)%l(1)==confs(2)%l(1)).and.(confs(1)%l(2)==confs(2)%l(2))) then
            res = res + S(confs(1)%n(1),confs(2)%n(1))*S(confs(1)%n(2),confs(2)%n(2))
        end if

        if ((confs(1)%l(1)==confs(2)%l(2)).and.(confs(1)%l(2)==confs(2)%l(1))) then
            res = res + (-1)**(L+sum(confs(2)%l))*S(confs(1)%n(1),confs(2)%n(2))*S(confs(1)%n(2),confs(2)%n(1))
        end if
    end function S_mat_neq

    pure function H_1p_eq(both,confs,H,S) result(res)
        logical, intent(in) :: both
        type(config), dimension(2), intent(in) :: confs
        type(block), dimension(:), allocatable, intent(in) :: H
        double complex, dimension(:,:), intent(in) :: S
        double complex :: res

        if ((confs(1)%l(1)/=confs(2)%l(1)).or.(confs(1)%l(2)/=confs(2)%l(2))) then
            res = 0.d0
            return
        end if

        res = H(confs(1)%l(1))%data(confs(1)%n(1),confs(2)%n(1))*S(confs(1)%n(2),confs(2)%n(2)) &
              + H(confs(1)%l(2))%data(confs(1)%n(2),confs(2)%n(2))*S(confs(1)%n(1),confs(2)%n(1))
        if (.not.both) res = sqrt(2.d0)*res
    end function H_1p_eq

    pure function H_1p_neq(confs,L,H,S) result(res)
        type(config), dimension(2), intent(in) :: confs
        integer, intent(in) :: L
        type(block), dimension(:), allocatable, intent(in) :: H
        double complex, dimension(:,:), intent(in) :: S
        double complex :: res

        res = 0.d0
        if ((confs(1)%l(1)==confs(2)%l(1)).and.(confs(1)%l(2)==confs(2)%l(2))) then
            res = res + H(confs(1)%l(1))%data(confs(1)%n(1),confs(2)%n(1))*S(confs(1)%n(2),confs(2)%n(2)) &
            + H(confs(1)%l(2))%data(confs(1)%n(2),confs(2)%n(2))*S(confs(1)%n(1),confs(2)%n(1))
        end if

        if ((confs(1)%l(1)==confs(2)%l(2)).and.(confs(1)%l(2)==confs(2)%l(1))) then
            res = res + (H(confs(1)%l(1))%data(confs(1)%n(1),confs(2)%n(2))*S(confs(1)%n(2),confs(2)%n(1)) &
            + H(confs(1)%l(2))%data(confs(1)%n(2),confs(2)%n(1))*S(confs(1)%n(1),confs(2)%n(2))) &
            *(-1)**(L+sum(confs(2)%l))
        end if
    end function H_1p_neq

    pure function dip_mat_neq(syms,confs,S,radial_dip) result(res)
        type(sym), dimension(2), intent(in) :: syms
        type(config), dimension(2), intent(in) :: confs
        double complex, dimension(:,:), intent(in) :: S
        type(radial_dipole), intent(in) :: radial_dip
        double complex ::  res

        double complex :: red_mat, red_mat_ex
        type(config), dimension(2) :: confs_ex

        ! Compute the reduced matrix element part
        red_mat = dip_red_mat(syms(1)%l,syms(2)%l,confs,S,radial_dip)

        ! Compute the exchange part of the reduced matrix element
        confs_ex(1) = confs(1)
        confs_ex(2)%n = confs(2)%n([2,1])
        confs_ex(2)%l = confs(2)%l([2,1])
        red_mat_ex = (-1)**(syms(2)%l+sum(confs(2)%l))*dip_red_mat(syms(1)%l,syms(2)%l,confs_ex,S,radial_dip)

        res = red_mat + red_mat_ex
    end function dip_mat_neq

    pure function dip_red_mat(L_1,L_2,confs,S,radial_dip) result(res)
        integer, intent(in) :: L_1
        integer, intent(in) :: L_2
        type(config), dimension(2), intent(in) :: confs
        double complex, dimension(:,:), intent(in) :: S
        type(radial_dipole), intent(in) :: radial_dip
        double complex :: res

        double complex :: res_1, res_2

        res_1 = 0.d0
        res_2 = 0.d0

        ! These are derived from Lindgren and Morrison Eq. (4.23)
        if (confs(1)%l(2) == confs(2)%l(2)) then
            res_1 = (-1)**(confs(1)%l(1) + confs(1)%l(2) + L_2 + 1) &
                *sqrt(real((2*L_1+1)*(2*L_2+1),kind=8)) &
                *six_j(L_1,1,L_2,confs(2)%l(1),confs(1)%l(2),confs(1)%l(1)) &
                *dip_red_1p(confs,1,radial_dip) &
                *S(confs(1)%n(2),confs(2)%n(2))
        end if

        if (confs(1)%l(1) == confs(2)%l(1)) then
            res_2 = (-1)**(confs(1)%l(1) + confs(2)%l(2) + L_1 + 1) &
                *sqrt(real((2*L_1+1)*(2*L_2+1),kind=8)) &
                *six_j(L_1,1,L_2,confs(2)%l(2),confs(1)%l(1),confs(1)%l(2)) &
                *dip_red_1p(confs,2,radial_dip) &
                *S(confs(1)%n(1),confs(2)%n(1))
        end if

        res = res_1 + res_2
    end function dip_red_mat

    pure function dip_red_1p_len(confs,i,radial_dip) result(res)
        type(config), dimension(2), intent(in) :: confs
        integer, intent(in) :: i
        type(radial_dipole), intent(in) :: radial_dip
        double complex :: res

        integer :: n,n_p,l,l_p

        n = confs(1)%n(i)
        n_p = confs(2)%n(i)
        l = confs(1)%l(i)
        l_p = confs(2)%l(i)

        res = radial_dip%r_mat(n,n_p)*C_red_mat(1,l,l_p)
    end function dip_red_1p_len

    pure function dip_red_1p_vel(confs,i,radial_dip) result(res)
        type(config), dimension(2), intent(in) :: confs
        integer, intent(in) :: i
        type(radial_dipole), intent(in) :: radial_dip
        double complex :: res

        integer :: n,n_p,l,l_p

        n = confs(1)%n(i)
        n_p = confs(2)%n(i)
        l = confs(1)%l(i)
        l_p = confs(2)%l(i)

        ! Quantum theory of angular momentum, Varshalovich et al p.486
        if (abs(l-l_p)/=1) then
            res = 0.d0
        else if (l>l_p) then
            res = sqrt(real(l,kind=8))*(radial_dip%dr_mat(n,n_p)-(l_p+1)*radial_dip%r_inv_mat(n,n_p))
        else !l<l_p
            res = -sqrt(real(l_p,kind=8))*(radial_dip%dr_mat(n,n_p)+l_p*radial_dip%r_inv_mat(n,n_p))
        end if
    end function dip_red_1p_vel

    pure function ang_dip_red(L_1,L_2,confs) result(res)
        integer, intent(in) :: L_1
        integer, intent(in) :: L_2
        type(config), dimension(2), intent(in) :: confs
        double precision :: res

        res = 0.d0

        if (confs(1)%l(2) == confs(2)%l(2)) then
            res = res + abs(six_j(L_1,1,L_2,confs(2)%l(1),confs(1)%l(2),confs(1)%l(1)) &
            *C_red_mat(1,confs(1)%l(1),confs(2)%l(1)))
        end if

        if (confs(1)%l(1) == confs(2)%l(1)) then
            res = res + abs(six_j(L_1,1,L_2,confs(2)%l(2),confs(1)%l(1),confs(1)%l(2)) &
            *C_red_mat(1,confs(1)%l(2),confs(2)%l(2)))
        end if
    end function ang_dip_red
end module mat_els