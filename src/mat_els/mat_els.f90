module mat_els
    use quad_tools
    use bspline_tools
    use grid_tools
    use sparse_array_tools
    use potentials
    use CAP_tools
    implicit none

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
                        call compute_H(potential,CAP_c,l,b_splines,i_b,j_b,c_i,c_j,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),H)
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
                        call compute_S(b_splines,i_b,j_b,c_i,c_j,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),S)
                    end if
                end do
            end do
        end do
        
    end subroutine setup_S

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
        c_j = 0.d0
        do k = 0,max_k
            ptr = 0
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
        c = 0.d0
        i = 0
        do k = 0,max_k
            ptr = 0
            do j_b_p = 1, b_splines%n-2
                c(:,3) = 0.d0
                c(j_b_p+1,4) = 1.d0
                c(j_b_p,4) = 0.d0
                i(4) = j_b_p
                do i_b_p = 1, b_splines%n-2
                    c(:,2) = 0.d0
                    c(i_b_p+1,3) = 1.d0
                    c(i_b_p,3) = 0.d0
                    i(3) = i_b_p
                    do j_b = 1, b_splines%n-2
                        c(j_b+1,2) = 1.d0
                        c(j_b,2) = 0.d0
                        i(2) = j_b

                        if (abs(j_b-j_b_p)>=b_splines%k) then
                            cycle
                        end if

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
    end subroutine setup_Slater_diag

    subroutine compute_H(potential,CAP_c,l,b_splines,i,j,c_i,c_j,k_GL,r,w,H)
        ! This is way too many arguments. Need to think about how to do this smarter.
        class(sph_pot), intent(in) :: potential
        class(CAP), intent(in) :: CAP_c
        integer, intent(in) :: l
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        double complex, dimension(:,:), intent(out) :: H


        integer :: i_sum
        double precision :: B_i, B_j, D_B_j

        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i)
            B_j = b_splines%eval_d(r(i_sum),c_j)
            D_B_j = b_splines%d_eval_d(r(i_sum),c_j,2)
            H(i,j) = H(i,j) + w(i_sum)*(-0.5*B_i*D_B_j + (potential%V(r(i_sum),l) + CAP_c%V(r(i_sum)))*B_i*B_j)
        end do

    end subroutine compute_H

    subroutine compute_S(b_splines,i,j,c_i,c_j,k_GL,r,w,S)
        ! This is way too many arguments. Need to think about how to do this smarter.
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        double complex, dimension(:,:), intent(out) :: S


        integer :: i_sum
        double precision :: B_i, B_j

        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i)
            B_j = b_splines%eval_d(r(i_sum),c_j)
            S(i,j) = S(i,j) + w(i_sum)*B_i*B_j
        end do

    end subroutine compute_S

    subroutine compute_radial_dip(b_splines,i,j,c_i,c_j,k_GL,r,w,r_mat,dr_mat)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i
        integer, intent(in) :: j
        double precision, dimension(b_splines%n), intent(in) :: c_i
        double precision, dimension(b_splines%n), intent(in) :: c_j
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: r
        double precision, dimension(k_GL), intent(in) :: w
        double complex, dimension(:,:), intent(out) :: r_mat
        double complex, dimension(:,:), intent(out) :: dr_mat

        integer :: i_sum
        double precision :: B_i,B_j,D_B_j
        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i)
            B_j = b_splines%eval_d(r(i_sum),c_j)
            D_B_j = b_splines%d_eval_d(r(i_sum),c_j,1)
            r_mat(i,j) = r_mat(i,j) + w(i_sum)*r(i_sum)*B_i*B_j
            dr_mat(i,j) = dr_mat(i,j) + dcmplx(0.d0,-1.d0)*w(i_sum)*B_i*D_B_j
        end do
    end subroutine compute_radial_dip

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
            B_i(i_sum) = b_splines%eval_d(r(i_sum),c_i)
            B_j(i_sum) = b_splines%eval_d(r(i_sum),c_j)
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
            B_i(i_sum) = b_splines%eval_d(r(i_sum),c(:,1))
            B_i_p(i_sum) = b_splines%eval_d(r(i_sum),c(:,3))
        end do

        do i_sum = 1,k_GL
            scale_j = 0.5d0*(r(i_sum)-limits(1))
            translate_j = 0.5d0*(r(i_sum)+limits(1))
            int_j = 0.d0
            do j_sum = 1,k_GL
                r_j = scale_j*x(j_sum) + translate_j
                B_j = b_splines%eval_d(r_j,c(:,2))
                B_j_p = b_splines%eval_d(r_j,c(:,4))
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
            r_d_k%j(ptr) = i(2)
            r_d_k%i_p(ptr) = i(3)
            r_d_k%j_p(ptr) = i(4)
        end if
    end subroutine compute_Slater_diag

    function Slater_k(n_b,a,b,c,d,k,r_d_k,r_k,r_m_k) result(res)
        integer, intent(in) :: n_b
        double complex, dimension(n_b), intent(in) :: a,b,c,d
        integer, intent(in) :: k
        type(sparse_6d), intent(in) :: r_d_k
        type(sparse_4d), intent(in) :: r_k
        type(sparse_4d), intent(in) :: r_m_k

        double complex :: res
        double complex :: coeffs
        integer i,j


        res = 0.d0
        do i = 1, r_d_k%nnz
            coeffs = a(r_d_k%i(i))*b(r_d_k%j(i))*c(r_d_k%i_p(i))*d(r_d_k%j_p(i)) &
                    + b(r_d_k%i(i))*a(r_d_k%j(i))*d(r_d_k%i_p(i))*c(r_d_k%j_p(i))
            res = res + coeffs*r_d_k%data(i,k)
        end do

        do i = 1,r_k%nnz
            do j = 1,r_k%nnz
                if (r_k%iv(i) < r_k%iv(j)) then
                    coeffs = a(r_k%i(i))*b(r_k%i(j))*c(r_k%j(i))*d(r_k%j(j)) &
                    + b(r_k%i(i))*a(r_k%i(j))*d(r_k%j(i))*c(r_k%j(j))
                    res = res + coeffs*r_k%data(i,k)*r_m_k%data(j,k)
                end if
            end do
        end do

    end function Slater_k
end module mat_els