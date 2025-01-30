module bspline_tools
    implicit none

    type, public :: b_spline
        integer :: k
        integer :: n
        integer :: n_b
        double precision, dimension(:), allocatable :: knots
        double precision, dimension(:), allocatable :: breakpoints
        logical, dimension(:,:), allocatable :: support

    contains
        procedure :: init
        procedure :: init_support
        procedure :: eval_d
        procedure :: d_eval_d
        procedure :: eval_z
        procedure :: d_eval_z

    end type b_spline

contains

    subroutine init(this,k,t)
        class(b_spline), intent(inout) :: this
        integer, intent(in) ::  k
        double precision, dimension(:), intent(in) :: t

        this%k = k
        this%n = size(t)-k
        this%n_b = this%n-2
        this%knots = t
        this%breakpoints = t(k:size(t) - k + 1)

        call this%init_support()

    end subroutine init

    subroutine init_support(this)
        class(b_spline), intent(inout) :: this

        integer :: i,j
        allocate(this%support(size(this%breakpoints)-1,this%n))

        this%support = .false.

        do i = 1, size(this%breakpoints)-1
            do j = i,i+this%k-1
                this%support(i,j) = .true.
            end do
        end do

    end subroutine init_support

    function eval_d(this,x,c) result(retval)
        class(b_spline), intent(in) :: this
        double precision, intent(in) :: x
        double precision, dimension(this%n), intent(in) :: c

        double precision :: retval

        retval = BVALUE_D(this%knots,c,this%n,this%k,x,0)

    end function eval_d

    function d_eval_d(this,x,c,n_deriv) result(retval)
        class(b_spline), intent(in) :: this
        double precision, intent(in) :: x
        double precision, dimension(this%n), intent(in) :: c
        integer, intent(in) :: n_deriv

        double precision :: retval

        retval = BVALUE_D(this%knots,c,this%n,this%k,x,n_deriv)

    end function d_eval_d

    function eval_z(this,x,c) result(retval)
        class(b_spline), intent(in) :: this
        double precision, intent(in) :: x
        double complex, dimension(this%n), intent(in) :: c

        double complex :: retval

        retval = BVALUE_Z(this%knots,c,this%n,this%k,x,0)

    end function eval_z

    function d_eval_z(this,x,c,n_deriv) result(retval)
        class(b_spline), intent(in) :: this
        double precision, intent(in) :: x
        double complex, dimension(this%n), intent(in) :: c
        integer, intent(in) :: n_deriv

        double complex :: retval

        retval = BVALUE_Z(this%knots,c,this%n,this%k,x,n_deriv)

    end function d_eval_z

    ! subroutine BSPLDR(t,a,n,k,ADIF,NDERIV)
    !     double precision, dimension(:), intent(in) :: t
    !     double complex, dimension(n),intent(in) :: a
    !     integer,intent(in) :: n
    !     integer,intent(in) :: k
    !     double precision ,dimension(n,NDERIV),intent(inout) :: ADIF
    !     integer, intent(in) :: NDERIV
    ! end subroutine BSPLDR

    subroutine BSPLEV(t,a,n,k,x,s_value,n_deriv)
        double precision, dimension(:), intent(in) :: t
        double complex, dimension(n),intent(in) :: a
        integer,intent(in) :: n
        integer,intent(in) :: k
        double precision, intent(in) :: x
        double complex, dimension(n_deriv), intent(out) :: s_value
        integer, intent(in) :: n_deriv

        integer :: i,id,m_flag

        id = max(min(n_deriv,k),1)
        s_value(1:id) = dcmplx(0.d0,0.d0)

        call interv(t,n+1,x,i,m_flag)

    end subroutine BSPLEV

    function BVALUE_D(t,a,n,k,x,i_deriv) result(val)
        double precision, dimension(n+k), intent(in) :: t
        double precision, dimension(n),intent(in) :: a
        integer,intent(in) :: n
        integer,intent(in) :: k
        double precision, intent(in) :: x
        integer, intent(in) :: i_deriv

        double precision :: val
        integer :: j,jj,i,m_flag
        integer :: i_p_1,i_p_j,i_p_1_m_j,i_m_k,i_m_k_p_j,k_m_i_der,k_m_1,k_m_j,i_h_i,i_h_m_k_m_j
        double precision, dimension(20) :: a_j,d_p,d_m

        val = dcmplx(0.d0,0.d0)
        k_m_i_der = k-i_deriv

        if (k_m_i_der <= 0) then
            return
        end if

        k_m_1 = k - 1
        call INTERV(t,n+1,x,i,m_flag)

        if (m_flag /= 0) then
            if (x < t(k)) return
            if (x > t(i)) return

            do while(x == t(i))
                i = i - 1
                if (i == k) return
            end do
        end if

        i_m_k = i - k
        do j = 1,k
            i_m_k_p_j = i_m_k + j
            a_j(j) = a(i_m_k_p_j)
        end do

        do j = 1,i_deriv
            k_m_j = k - j
            do jj = 1,k_m_j
                i_h_i = i + jj
                i_h_m_k_m_j = i_h_i - k_m_j
                a_j(jj) = (a_j(jj+1)-a_j(jj))/(t(i_h_i)-t(i_h_m_k_m_j))*k_m_j
            end do
        end do

        i_p_1 = i + 1
        do j = 1,k_m_i_der
            i_p_j = i + j
            d_p(j) = t(i_p_j) - x
            i_p_1_m_j = i_p_1 - j
            d_m(j) = x - t(i_p_1_m_j)
        end do

        do j = i_deriv + 1, k_m_1
            k_m_j = k - j
            i = k_m_j
            do jj = 1, k_m_j
                a_j(jj) = (a_j(jj+1)*d_m(i) + a_j(jj)*d_p(jj))/(d_m(i)+d_p(jj))
                i = i-1
            end do
        end do
        val = a_j(1) !I think this is correct, but not sure
        return

    end function BVALUE_D

    function BVALUE_Z(t,a,n,k,x,i_deriv) result(val)
        double precision, dimension(n+k), intent(in) :: t
        double complex, dimension(n),intent(in) :: a
        integer,intent(in) :: n
        integer,intent(in) :: k
        double precision, intent(in) :: x
        integer, intent(in) :: i_deriv

        double complex :: val
        integer :: j,jj,i,m_flag
        integer :: i_p_1,i_p_j,i_p_1_m_j,i_m_k,i_m_k_p_j,k_m_i_der,k_m_1,k_m_j,i_h_i,i_h_m_k_m_j
        double complex, dimension(20) :: a_j,d_p,d_m

        val = dcmplx(0.d0,0.d0)
        k_m_i_der = k-i_deriv

        if (k_m_i_der <= 0) then
            return
        end if

        k_m_1 = k - 1
        call INTERV(t,n+1,x,i,m_flag)

        if (m_flag /= 0) then
            if (x < t(k)) return
            if (x > t(i)) return

            do while(x == t(i))
                i = i - 1
                if (i == k) return
            end do
        end if

        i_m_k = i - k
        do j = 1,k
            i_m_k_p_j = i_m_k + j
            a_j(j) = a(i_m_k_p_j)
        end do

        do j = 1,i_deriv
            k_m_j = k - j
            do jj = 1,k_m_j
                i_h_i = i + jj
                i_h_m_k_m_j = i_h_i - k_m_j
                a_j(jj) = (a_j(jj+1)-a_j(jj))/(t(i_h_i)-t(i_h_m_k_m_j))*k_m_j
            end do
        end do

        i_p_1 = i + 1
        do j = 1,k_m_i_der
            i_p_j = i + j
            d_p(j) = t(i_p_j) - x
            i_p_1_m_j = i_p_1 - j
            d_m(j) = x - t(i_p_1_m_j)
        end do

        do j = i_deriv + 1, k_m_1
            k_m_j = k - j
            i = k_m_j
            do jj = 1, k_m_j
                a_j(jj) = (a_j(jj+1)*d_m(i) + a_j(jj)*d_p(jj))/(d_m(i)+d_p(jj))
                i = i-1
            end do
        end do
        val = a_j(1) !I think this is correct, but not sure
        return

    end function BVALUE_Z

    subroutine INTERV(xt,lxt,x,i_left,m_flag)
        double precision, dimension(lxt), intent(in) :: xt
        integer, intent(in) :: lxt
        double precision, intent(in) :: x
        integer, intent(out) :: i_left
        integer, intent(out) :: m_flag

        integer, save :: i_lo = 1
        integer :: i_hi


        if (x < xt(1)) then
            i_left = 1
            m_flag = -1
        else if (xt(lxt) <= x) then
            i_left = lxt
            m_flag = 1
        else
            if (i_lo .lt. lxt) then
                i_hi = i_lo + 1

            else
                i_lo = lxt -1
                i_hi = lxt
            end if

            if ((xt(i_lo) <= x) .and. (x < xt(i_hi))) then
                i_left = i_lo
                m_flag = 0
            else
                if ( xt(i_lo) > x ) then
                    do while (xt(i_lo) > x)
                        i_lo = i_lo - 1
                    end do
                    i_left = i_lo
                    m_flag = 0
                else
                    do while(xt(i_hi) <= x)
                        i_hi = i_hi + 1
                    end do
                    i_lo = i_hi-1
                    i_left = i_lo
                    m_flag = 0
                end if

            end if
        end if
    end subroutine INTERV

end module bspline_tools