module mat_els
    use quad_tools
    use bspline_tools
    use grid_tools
    use potentials
    implicit none

contains
    subroutine setup_H_one_particle(potential,l,b_splines,H,S)
        class(sph_pot), intent(in) :: potential
        integer, intent(in) :: l
        type(b_spline), intent(in) :: b_splines
        double complex, dimension(:,:), intent(out) :: H
        double complex, dimension(:,:), intent(out) :: S

        integer :: i_b,j_b,i_r
        integer :: k_GL
        double precision, dimension(:), allocatable :: c_i,c_j
        type(gau_leg) :: g_l

        ! k_GL = 2*b_splines%k-2 ! exact for overlap
        k_GL = b_splines%k + 3
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
                        call compute_H_and_S(potential,l,b_splines,i_b,j_b,c_i,c_j,k_GL,g_l%x(:,i_r),g_l%w(:,i_r),H,S)
                    end if
                end do
            end do
        end do
    end subroutine setup_H_one_particle

    subroutine setup_Slater_integrals(b_splines,max_k,k_GL,r_k,r_m_k,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        double precision, dimension(:,:,:,:), intent(out) :: r_k
        double precision, dimension(:,:,:,:), intent(out) :: r_m_k
        double precision, dimension(:,:,:,:,:,:), intent(out) :: r_d_k

        call setup_Slater_off_diag(b_splines,max_k,k_GL,r_k,r_m_k)
        call setup_Slater_diag(b_splines,max_k,k_GL,r_d_k)
    end subroutine setup_Slater_integrals

    subroutine setup_Slater_off_diag(b_splines,max_k,k_GL,r_k,r_m_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        double precision, dimension(:,:,:,:), intent(out) :: r_k
        double precision, dimension(:,:,:,:), intent(out) :: r_m_k

        integer :: i_b,j_b,i_r
        double precision, dimension(k_GL) :: x,w
        double precision, dimension(2) :: limits
        double precision, dimension(:), allocatable :: c_i,c_j

        call setup_GL(k_GL,-1.d0,1.d0,x,w)

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
                        limits(1) = b_splines%breakpoints(i_r)
                        limits(2) = b_splines%breakpoints(i_r+1)
                        call compute_Slater_off_diag(b_splines,i_r,i_b,j_b,c_i,c_j,k_GL,x,w,limits,max_k,r_k,r_m_k)
                    end if
                end do
            end do
        end do
    end subroutine setup_Slater_off_diag

    subroutine setup_Slater_diag(b_splines,max_k,k_GL,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        integer, intent(in) :: k_GL
        double precision, dimension(:,:,:,:,:,:), intent(out) :: r_d_k

        integer :: i_b,j_b,i_b_p,j_b_p,i_r
        integer, dimension(4) :: i
        double precision, dimension(k_GL) :: x,w
        double precision, dimension(2) :: limits
        double precision, dimension(:,:), allocatable :: c
        logical :: j_support, i_support

        call setup_GL(k_GL,-1.d0,1.d0,x,w)

        allocate(c(b_splines%n,4))
        c = 0.d0
        i = 0
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
                                limits(1) = b_splines%breakpoints(i_r)
                                limits(2) = b_splines%breakpoints(i_r+1)
                                call compute_Slater_diag(b_splines,i_r,i,c,k_GL,x,w,limits,max_k,r_d_k)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine setup_Slater_diag

    subroutine compute_H_and_S(potential,l,b_splines,i,j,c_i,c_j,k_GL,r,w,H,S)
        ! This is way too many arguments. Need to think about how to do this smarter.
        class(sph_pot), intent(in) :: potential
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
        double complex, dimension(:,:), intent(out) :: S


        integer :: i_sum
        double precision :: B_i, B_j, D_B_j

        do i_sum = 1, k_GL
            B_i = b_splines%eval_d(r(i_sum),c_i)
            B_j = b_splines%eval_d(r(i_sum),c_j)
            D_B_j = b_splines%d_eval_d(r(i_sum),c_j,2)
            H(i,j) = H(i,j) + w(i_sum)*(-0.5*B_i*D_B_j + potential%V(r(i_sum),l)*B_i*B_j)
            S(i,j) = S(i,j) + w(i_sum)*B_i*B_j
        end do

    end subroutine compute_H_and_S

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

    subroutine compute_Slater_off_diag(b_splines,i_r,i,j,c_i,c_j,k_GL,x,w,limits,max_k,r_k,r_m_k)
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
        integer, intent(in) :: max_k
        double precision, dimension(:,:,:,:), intent(out) :: r_k
        double precision, dimension(:,:,:,:), intent(out) :: r_m_k

        integer :: k,i_sum
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

        do k = 0,max_k
            do i_sum = 1,k_GL
                r_k(k+1,i_r,i,j) = r_k(k+1,i_r,i,j) + w(i_sum)*B_i(i_sum)*B_j(i_sum)*r(i_sum)**k
                r_m_k(k+1,i_r,i,j) = r_m_k(k+1,i_r,i,j) + w(i_sum)*B_i(i_sum)*B_j(i_sum)/r(i_sum)**(k+1)
            end do
            !Apply jacobian
            r_k(k+1,i_r,i,j) = scale*r_k(k+1,i_r,i,j)
            r_m_k(k+1,i_r,i,j) = scale*r_m_k(k+1,i_r,i,j)
        end do
    end subroutine compute_Slater_off_diag

    subroutine compute_Slater_diag(b_splines,i_r,i,c,k_GL,x,w,limits,max_k,r_d_k)
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: i_r
        integer, dimension(4), intent(in) :: i
        double precision, dimension(b_splines%n,4), intent(in) :: c
        integer, intent(in) :: k_GL
        double precision, dimension(k_GL), intent(in) :: x
        double precision, dimension(k_GL), intent(in) :: w
        double precision, dimension(2), intent(in) :: limits
        integer, intent(in) :: max_k
        double precision, dimension(:,:,:,:,:,:), intent(out) :: r_d_k

        integer :: k,i_sum,j_sum
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

        do k = 0,max_k
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
                r_d_k(k+1,i_r,i(1),i(2),i(3),i(4)) = r_d_k(k+1,i_r,i(1),i(2),i(3),i(4)) &
                                                    + w(i_sum)*B_i(i_sum)*B_i_p(i_sum)*int_j/r(i_sum)**(k+1)
            end do
            r_d_k(k+1,i_r,i(1),i(2),i(3),i(4)) = scale_i*r_d_k(k+1,i_r,i(1),i(2),i(3),i(4))
        end do
    end subroutine compute_Slater_diag
end module mat_els