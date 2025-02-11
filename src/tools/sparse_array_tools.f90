module sparse_array_tools
    use bspline_tools
    implicit none

    type, public :: sparse_4d
    integer :: nnz
    double precision, dimension(:,:), allocatable :: data
    integer, dimension(:), allocatable :: iv
    integer, dimension(:), allocatable :: i
    integer, dimension(:), allocatable :: j
    contains
        procedure :: init => init_4d
        procedure :: count_nnz => count_nnz_4d
    end type sparse_4d

    type, public :: sparse_6d
        integer :: nnz
        double precision, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: iv
        integer, dimension(:), allocatable :: i
        integer, dimension(:), allocatable :: j
        integer, dimension(:), allocatable :: i_p
        integer, dimension(:), allocatable :: j_p
    contains
        procedure :: init => init_6d
        procedure :: count_nnz => count_nnz_6d
    end type sparse_6d
contains

    subroutine init_4d(this,max_k,b_splines)
        class(sparse_4d), intent(inout) :: this
        integer, intent(in) ::  max_k
        type(b_spline), intent(in) :: b_splines

        call count_nnz_4d(this,b_splines)
        allocate(this%data(this%nnz,0:max_k),source=0.d0)
        allocate(this%iv(this%nnz),this%i(this%nnz),this%j(this%nnz),source = 0)
    end subroutine init_4d

    subroutine init_6d(this,max_k,b_splines)
        class(sparse_6d), intent(inout) :: this
        integer, intent(in) ::  max_k
        type(b_spline), intent(in) :: b_splines

        call count_nnz_6d(this,b_splines)
        allocate(this%data(this%nnz,0:max_k),source=0.d0)
        allocate(this%iv(this%nnz),this%i(this%nnz),&
                this%j(this%nnz),this%i_p(this%nnz),this%j_p(this%nnz),source = 0)
    end subroutine init_6d

    subroutine count_nnz_4d(this, b_splines)
        class(sparse_4d), intent(inout) :: this
        type(b_spline), intent(in) :: b_splines

        integer :: i,j,iv
        logical :: support

        this%nnz = 0

        do j = 1, b_splines%n_b
            do i = 1,b_splines%n_b
                if (abs(i-j) >= b_splines%k) cycle
                do iv = 1,size(b_splines%breakpoints)-1
                    support = (b_splines%support(iv,i+1).and.b_splines%support(iv,j+1))
                    if (support) then
                        this%nnz = this%nnz + 1
                    end if
                end do
            end do
        end do
    end subroutine count_nnz_4d

    subroutine count_nnz_6d(this, b_splines)
        class(sparse_6d), intent(inout) :: this
        type(b_spline), intent(in) :: b_splines

        integer :: i,j,i_p,j_p,iv
        logical :: support_i, support_j

        this%nnz = 0

        do j_p = 1, b_splines%n_b
            do i_p = 1, b_splines%n_b
                do j = 1, b_splines%n_b
                    if (abs(j-j_p)>=b_splines%k) cycle
                    do i = 1,b_splines%n_b
                        if (abs(i-i_p) >= b_splines%k) cycle
                        do iv = 1,size(b_splines%breakpoints)-1
                            support_i = (b_splines%support(iv,i+1).and.b_splines%support(iv,i_p+1))
                            support_j = (b_splines%support(iv,j+1).and.b_splines%support(iv,j_p+1))
                            if (support_i.and.support_j) then
                                this%nnz = this%nnz + 1
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine count_nnz_6d
end module sparse_array_tools