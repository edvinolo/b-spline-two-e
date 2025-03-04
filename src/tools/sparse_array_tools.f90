module sparse_array_tools
    use bspline_tools
    use stdlib_sorting, only: sort_index
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

    type, public :: sparse_Slater
        integer :: nnz
        double precision, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: i
        integer, dimension(:), allocatable :: j
        integer, dimension(:), allocatable :: i_p
        integer, dimension(:), allocatable :: j_p
    contains
        procedure :: init => init_sparse_Slater
        procedure :: count_nnz => count_nnz_Slater
    end type sparse_Slater
contains

    pure subroutine init_4d(this,max_k,b_splines)
        class(sparse_4d), intent(inout) :: this
        integer, intent(in) ::  max_k
        type(b_spline), intent(in) :: b_splines

        call count_nnz_4d(this,b_splines)
        allocate(this%data(this%nnz,0:max_k),source=0.d0)
        allocate(this%iv(this%nnz),this%i(this%nnz),this%j(this%nnz),source = 0)
    end subroutine init_4d

    pure subroutine init_6d(this,max_k,b_splines)
        class(sparse_6d), intent(inout) :: this
        integer, intent(in) ::  max_k
        type(b_spline), intent(in) :: b_splines

        call count_nnz_6d(this,b_splines)
        allocate(this%data(this%nnz,0:max_k),source=0.d0)
        allocate(this%iv(this%nnz),this%i(this%nnz),&
                this%j(this%nnz),this%i_p(this%nnz),this%j_p(this%nnz),source = 0)
    end subroutine init_6d

    subroutine init_sparse_Slater(this,max_k,b_splines,r_d_k,r_k,r_m_k)
        class(sparse_Slater), intent(inout) :: this
        integer, intent(in) :: max_k
        type(b_spline), intent(in) :: b_splines
        type(sparse_6d), intent(in) :: r_d_k
        type(sparse_4d), intent(in) :: r_k
        type(sparse_4d), intent(in) :: r_m_k

        integer :: k,ptr,ptr_d,ptr_j,ptr_i,i,j,i_p,j_p,iv,jv
        integer :: ptr_t
        logical :: i_support,j_support
        type(sparse_6d) :: r_d_k_P

        call count_nnz_Slater(this,b_splines)
        allocate(this%data(this%nnz,0:max_k),source=0.d0)
        allocate(this%i(this%nnz),this%j(this%nnz),this%i_p(this%nnz),this%j_p(this%nnz),source = 0)

        r_d_k_P = permute_6d(r_d_k)
        write(6,*) r_k%nnz,r_d_k%nnz,this%nnz
        do k = 0,max_k
            ptr = 1
            ptr_d = 1
            do j_p = 1,b_splines%n_b
                do j = 1,b_splines%n_b
                    if (abs(j-j_p)>=b_splines%k) cycle
                    do i_p = 1,b_splines%n_b
                        do i = 1,b_splines%n_b
                            if (abs(i-i_p)>=b_splines%k) cycle
                            this%i(ptr) = i
                            this%i_p(ptr) = i_p
                            this%j(ptr) = j
                            this%j_p(ptr) = j_p
                            do iv = 1,size(b_splines%breakpoints)-1
                                j_support = b_splines%support(iv,j+1).and.b_splines%support(iv,j_p+1)
                                i_support = b_splines%support(iv,i+1).and.b_splines%support(iv,i_p+1)
                                if (i_support.and.j_support) then
                                    this%data(ptr,k) = this%data(ptr,k) &
                                           + r_d_k%data(ptr_d,k)
                                           !+ r_d_k_P%data(ptr_d,k)
                                    if (this%j_p(ptr)/=r_d_k%j_p(ptr_d)) then
                                       write(6,*) 'Different indices', j,r_d_k%j(ptr_d),ptr,ptr_d
                                       stop
                                    end if
                                    ptr_d = ptr_d + 1
                                end if
                            end do
                            ptr = ptr + 1
                        end do
                    end do
                end do
            end do
        end do

        ! do k = 0,max_k
        !     ptr = 1
        !     ptr_t = ptr
        !     ptr_d = 1
        !     ptr_j = 1
        !     do j_p = 1,b_splines%n_b
        !         do j = 1,b_splines%n_b
        !             if (abs(j-j_p)>=b_splines%k) cycle
        !             do jv = 1,size(b_splines%breakpoints)-1
        !                 j_support = b_splines%support(jv,j+1).and.b_splines%support(jv,j_p+1)
        !                 if (.not.j_support) cycle
        !                 ptr_i = 1
        !                 ptr_t = ptr
        !                 do i_p = 1,b_splines%n_b
        !                     do i = 1,b_splines%n_b
        !                         if (abs(i-i_p)>=b_splines%k) cycle
        !                         do iv = 1,size(b_splines%breakpoints)-1
        !                             i_support = b_splines%support(iv,i+1).and.b_splines%support(iv,i_p+1)
        !                             if (.not.i_support) cycle
        !                             if (iv<jv) then
        !                                 this%data(ptr_t,k) = this%data(ptr_t,k)&
        !                                         + r_k%data(ptr_i,k)*r_m_k%data(ptr_j,k)
        !                             else if (jv<iv) then
        !                                 this%data(ptr_t,k) = this%data(ptr_t,k) &
        !                                         + r_k%data(ptr_j,k)*r_m_k%data(ptr_i,k)
        !                             end if
        !                             ptr_i = ptr_i + 1
        !                         end do
        !                         ! this%i(ptr_t) = i
        !                         ! this%i_p(ptr_t) = i_p
        !                         ! this%j(ptr_t) = j
        !                         ! this%j_p(ptr_t) = j_p
        !                         ptr_t = ptr_t + 1
        !                     end do
        !                 end do
        !                 ptr_j = ptr_j + 1
        !             end do
        !             ptr = ptr_t
        !         end do
        !     end do
        ! end do
    end subroutine init_sparse_Slater

    pure subroutine count_nnz_4d(this, b_splines)
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

    pure subroutine count_nnz_6d(this, b_splines)
        class(sparse_6d), intent(inout) :: this
        type(b_spline), intent(in) :: b_splines

        integer :: i,j,i_p,j_p,iv
        logical :: support_i, support_j

        this%nnz = 0

        do j_p = 1, b_splines%n_b
            do j = 1, b_splines%n_b
                if (abs(j-j_p)>=b_splines%k) cycle
                do i_p = 1, b_splines%n_b
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

    pure subroutine count_nnz_Slater(this,b_splines)
        class(sparse_Slater), intent(inout) :: this
        type(b_spline), intent(in) :: b_splines

        integer :: i,j,i_p,j_p

        this%nnz = 0

        do j_p = 1, b_splines%n_b
            do j = 1, b_splines%n_b
                if (abs(j-j_p)>=b_splines%k) cycle
                do i_p = 1, b_splines%n_b
                    do i = 1,b_splines%n_b
                        if (abs(i-i_p) >= b_splines%k) cycle
                        this%nnz = this%nnz + 1
                    end do
                end do
            end do
        end do
    end subroutine count_nnz_Slater

    function permute_6d(A) result(res)
        type(sparse_6d), intent(in) :: A
        type(sparse_6d) :: res

        integer, dimension(:), allocatable :: index

        res%nnz = A%nnz
        if (allocated(A%data)) then
            res%data = A%data
        else
            write(6,*) "Data array not allocated in sparse Slate transpose!"
            stop
        end if
        res%iv = A%iv
        res%i = A%j
        res%j = A%i
        res%i_p = A%j_p
        res%j_p = A%i_p

        allocate(index(res%nnz))
        call sort_index(res%j_p,index)

        res%data = res%data(index,:)
        res%iv = res%iv(index)
        res%i = res%i(index)
        res%j = res%j(index)
        res%i_p = res%i_p(index)
    end function permute_6d

    subroutine compute_R_k(r_d_k,r_k,r_m_k,max_k,n_b,R)
        type(sparse_6d), intent(in) :: r_d_k
        type(sparse_4d), intent(in) :: r_k
        type(sparse_4d), intent(in) :: r_m_k
        integer, intent(in) :: max_k
        integer, intent(in) :: n_b
        double precision, dimension(0:max_k,n_b,n_b,n_b,n_b), intent(out) :: R

        integer :: n,m,k

        R = 0.d0
        do k = 0,max_k
            do n=1,r_k%nnz
                do m=1,r_k%nnz
                    if (r_k%iv(n)<r_k%iv(m)) then
                        R(k,r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)) =  R(k,r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)) &
                            + r_k%data(n,k)*r_m_k%data(m,k)
                    else if (r_k%iv(n)>r_k%iv(m)) then
                        R(k,r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)) = R(k,r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)) &
                            + r_k%data(m,k)*r_m_k%data(n,k)
                    end if
                end do
            end do

            do n=1,r_d_k%nnz
                R(k,r_d_k%i(n),r_d_k%j(n),r_d_k%i_p(n),r_d_k%j_p(n)) =  R(k,r_d_k%i(n),r_d_k%j(n),r_d_k%i_p(n),r_d_k%j_p(n))&
                    + r_d_k%data(n,k)

                R(k,r_d_k%j(n),r_d_k%i(n),r_d_k%j_p(n),r_d_k%i_p(n)) =  R(k,r_d_k%j(n),r_d_k%i(n),r_d_k%j_p(n),r_d_k%i_p(n))&
                    + r_d_k%data(n,k)
            end do
        end do
    end subroutine compute_R_k
end module sparse_array_tools