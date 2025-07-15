module block_tools
    use sparse_array_tools, only: CSR_matrix
    implicit none

    type, public :: block
        integer :: n
        integer :: m
        double complex, dimension(:,:), allocatable :: data
        logical :: nz
    contains
        procedure :: init => init_block
    end type block

    type, public :: block_mat
        integer :: n
        integer :: m
        type(block), dimension(:,:), allocatable :: blocks
    contains
        procedure :: init => init_block_mat
    end type block_mat

    type, public :: block_CS
        integer, dimension(2) :: block_shape
        integer, dimension(2) :: shape
        type(CSR_matrix), dimension(:,:), allocatable :: blocks
    contains
        procedure :: init => init_block_CS
        procedure :: deall => deallocate_block_CS
        procedure :: compute_shape => compute_shape_block_CS
        procedure :: to_CS => CS_from_block
        procedure :: store => CS_block_store
        procedure :: load => CS_block_load
        procedure :: scale => CS_block_scale
    end type block_CS

    type, public :: block_diag_CS
        integer, dimension(2) :: block_shape
        integer, dimension(2) :: shape
        type(CSR_matrix), dimension(:), allocatable :: blocks
    contains
        procedure :: init => init_block_diag_CS
        procedure :: deall => deallocate_block_diag_CS
        procedure :: compute_shape => compute_shape_block_diag_CS
        procedure :: to_CS => CS_from_block_diag
        procedure :: store => CS_block_diag_store
        procedure :: load => CS_block_diag_load
        procedure :: scale => CS_block_diag_scale
    end type  block_diag_CS

    interface XPAY
        module procedure :: CS_diag_X_P_A_block_Y
    end interface

    interface AXPBY
        module procedure :: CS_A_block_X_P_B_block_Y
    end interface

    interface APX
        module procedure :: CS_A_P_diag_X, CS_A_P_block_X, CS_A_diag_B_P_diag_X, CS_A_diag_B_P_block_X, &
                         CS_A_block_B_P_block_X
    end interface
contains
    subroutine init_block(this,n,m,nz)
        class(block), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: m
        logical, intent(in) :: nz

        this%n = n
        this%m = m
        this%nz = nz

        if (this%nz) then
            allocate(this%data(n,m),source=dcmplx(0.d0,0.d0))
        end if
    end subroutine init_block

    subroutine init_block_mat(this,n,m)
        class(block_mat), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: m

        this%n = n
        this%m = m

        allocate(this%blocks(n,m))
    end subroutine init_block_mat

    subroutine init_block_CS(this,shape_in)
        class(block_CS), intent(inout) :: this
        integer, dimension(2), intent(in) :: shape_in

        this%block_shape = shape_in
        allocate(this%blocks(shape_in(1),shape_in(2)))
        this%blocks%nnz = 0
    end subroutine init_block_CS

    subroutine init_block_diag_CS(this,n)
        class(block_diag_CS), intent(inout) :: this
        integer, intent(in) :: n

        this%block_shape = [n,n]
        allocate(this%blocks(n))
        this%blocks%nnz = 0
    end subroutine init_block_diag_CS

    subroutine deallocate_block_CS(this)
        class(block_CS), intent(inout) :: this

        if (allocated(this%blocks)) deallocate(this%blocks)
    end subroutine deallocate_block_CS

    subroutine deallocate_block_diag_CS(this)
        class(block_diag_CS), intent(inout) :: this

        if (allocated(this%blocks)) deallocate(this%blocks)
    end subroutine deallocate_block_diag_CS

    subroutine compute_shape_block_CS(this)
        class(block_CS), intent(inout) :: this

        integer :: i

        this%shape = 0
        do i = 1, this%block_shape(1)
            this%shape(1) = this%shape(1) + this%blocks(i,1)%shape(1)
        end do

        do i = 1, this%block_shape(2)
            this%shape(2) = this%shape(2) + this%blocks(1,i)%shape(2)
        end do
    end subroutine compute_shape_block_CS

    subroutine compute_shape_block_diag_CS(this)
        class(block_diag_CS), intent(inout) :: this

        integer :: i

        this%shape = 0
        do i = 1, this%block_shape(1)
            this%shape = this%shape + this%blocks(i)%shape
        end do
    end subroutine compute_shape_block_diag_CS

    subroutine CS_from_block(blocks,CS,deall_block)
        class(block_CS), intent(inout) :: blocks
        type(CSR_matrix), intent(inout) :: CS
        logical, intent(in) :: deall_block

        integer :: i, j, nnz

        nnz = 0
        do j = 1,blocks%block_shape(2)
            do i = 1,blocks%block_shape(1)
                if (blocks%blocks(i,j)%arrays_allocated()) then
                    nnz = nnz + blocks%blocks(i,j)%nnz
                end if
            end do
        end do

        call blocks%compute_shape()
        if ((CS%nnz /= nnz).or.(any(CS%shape /= blocks%shape))) then
            if (CS%arrays_allocated()) call CS%deall()
            call CS%init(blocks%shape,nnz)
        end if

        if (nnz /= 0) then
            call CSR_from_block(blocks,CS)
        end if

        if (deall_block) then
            do j = 1,blocks%block_shape(2)
                do i = 1,blocks%block_shape(1)
                    call blocks%blocks(i,j)%deall()
                end do
            end do
        end if
    end subroutine CS_from_block

    subroutine CSR_from_block(blocks,CS)
        class(block_CS), intent(inout) :: blocks
        type(CSR_matrix), intent(inout) :: CS

        integer :: i, j, row_i
        integer :: ptr, row, col, start, end, CS_start, CS_end

        ptr = 1
        row = 0

        do i = 1,blocks%block_shape(1)
            do row_i = 1, blocks%blocks(i,1)%shape(1)
                col = 0
                row = row + 1
                CS%index_ptr(row) = ptr
                do j = 1,blocks%block_shape(2)
                    if (blocks%blocks(i,j)%arrays_allocated()) then
                        start = blocks%blocks(i,j)%index_ptr(row_i)
                        end  = blocks%blocks(i,j)%index_ptr(row_i+1)-1
                        CS_start = ptr
                        CS_end = ptr + (end-start)
                        CS%indices(CS_start:CS_end) = col + blocks%blocks(i,j)%indices(start:end)
                        CS%data(CS_start:CS_end) = blocks%blocks(i,j)%data(start:end)
                        ptr = CS_end + 1
                    end if
                    col = col + blocks%blocks(i,j)%shape(2)
                end do
            end do
        end do
        CS%index_ptr(row+1) = ptr

        if (CS%index_ptr(CS%shape(1)+1) /= CS%nnz+1) then
            write(6,*) "Error in CSR from block assembly, index_ptr(N_rows+1) /= nnz + 1"
            write(6,*) CS%index_ptr(CS%shape(1)+1), CS%nnz+1
            stop
        end if

    end subroutine CSR_from_block

    ! subroutine CSC_from_block(blocks,CS)
    !     class(block_CS), intent(inout) :: blocks
    !     type(CSC_matrix), intent(inout) :: CS

    !     integer :: i, j, col_i
    !     integer :: ptr, row, col, start, end, CS_start, CS_end

    !     ptr = 1
    !     col = 0

    !     do i = 1,blocks%block_shape(2)
    !         do col_i = 1, blocks%blocks(1,i)%shape(2)
    !             col = col +1
    !             row = 0
    !             CS%index_ptr(col) = ptr
    !             do j = 1,blocks%block_shape(1)
    !                 if (blocks%blocks(j,i)%arrays_allocated()) then
    !                     start = blocks%blocks(j,i)%index_ptr(col_i)
    !                     end  = blocks%blocks(j,i)%index_ptr(col_i+1)-1
    !                     CS_start = ptr
    !                     CS_end = ptr + (end-start)
    !                     CS%indices(CS_start:CS_end) = col + blocks%blocks(j,i)%indices(start:end)
    !                     CS%data(CS_start:CS_end) = blocks%blocks(j,i)%data(start:end)
    !                     ptr = CS_end + 1
    !                 end if
    !                 row = row + blocks%blocks(j,i)%shape(1)
    !             end do
    !         end do
    !     end do
    !     CS%index_ptr(col+1) = ptr

    !     if (CS%index_ptr(CS%shape(2)+1) /= CS%nnz+1) then
    !         write(6,*) "Error in CSC from block assembly, index_ptr(N_cols+1) /= nnz + 1"
    !         write(6,*) CS%index_ptr(CS%shape(2)+1), CS%nnz+1
    !         stop
    !     end if
    ! end subroutine CSC_from_block

    subroutine CS_from_block_diag(blocks,CS,deall_block)
        class(block_diag_CS), intent(inout) :: blocks
        type(CSR_matrix), intent(inout) :: CS
        logical, intent(in) :: deall_block

        integer :: i, nnz

        nnz = 0
        do i = 1,blocks%block_shape(1)
            if (blocks%blocks(i)%arrays_allocated()) then
                nnz = nnz + blocks%blocks(i)%nnz
            end if
        end do

        call blocks%compute_shape()
        if ((CS%nnz /= nnz).or.(any(CS%shape /= blocks%shape))) then
            if (CS%arrays_allocated()) call CS%deall()
            call CS%init(blocks%shape,nnz)
        end if

        if (nnz /= 0) then
            call CSR_from_block_diag(blocks,CS)
        end if


        if (deall_block) then
            do i = 1,blocks%block_shape(1)
                call blocks%blocks(i)%deall()
            end do
        end if
    end subroutine CS_from_block_diag

    subroutine CSR_from_block_diag(blocks,CS)
        type(block_diag_CS), intent(inout) :: blocks
        type(CSR_matrix), intent(inout) :: CS

        integer :: i, row_i
        integer :: ptr, row, col, start, end, CS_start, CS_end

        ptr = 1
        row = 0

        col = 0
        do i = 1,blocks%block_shape(1)
            do row_i = 1, blocks%blocks(i)%shape(1)
                row = row + 1
                CS%index_ptr(row) = ptr

                if (blocks%blocks(i)%arrays_allocated()) then
                    start = blocks%blocks(i)%index_ptr(row_i)
                    end  = blocks%blocks(i)%index_ptr(row_i+1)-1
                    CS_start = ptr
                    CS_end = ptr + (end-start)
                    CS%indices(CS_start:CS_end) = col + blocks%blocks(i)%indices(start:end)
                    CS%data(CS_start:CS_end) = blocks%blocks(i)%data(start:end)
                    ptr = CS_end + 1
                end if
            end do
            col = col + blocks%blocks(i)%shape(2)
        end do
        CS%index_ptr(row+1) = ptr

        if (CS%index_ptr(CS%shape(1)+1) /= CS%nnz+1) then
            write(6,*) "Error in CSR from block assembly, index_ptr(N_rows+1) /= nnz + 1"
            stop
        end if
    end subroutine CSR_from_block_diag

    ! subroutine CSC_from_block_diag(blocks,CS)
    !     type(block_diag_CS), intent(inout) :: blocks
    !     type(CSC_matrix), intent(inout) :: CS

    !     integer :: i, col_i
    !     integer :: ptr, row, col, start, end, CS_start, CS_end

    !     ptr = 1
    !     col = 0

    !     row = 0
    !     do i = 1,blocks%block_shape(1)
    !         do col_i = 1, blocks%blocks(i)%shape(1)
    !             col = col + 1
    !             CS%index_ptr(col) = ptr

    !             if (blocks%blocks(i)%arrays_allocated()) then
    !                 start = blocks%blocks(i)%index_ptr(col_i)
    !                 end  = blocks%blocks(i)%index_ptr(col_i+1)-1
    !                 CS_start = ptr
    !                 CS_end = ptr + (end-start)
    !                 CS%indices(CS_start:CS_end) = col + blocks%blocks(i)%indices(start:end)
    !                 CS%data(CS_start:CS_end) = blocks%blocks(i)%data(start:end)
    !                 ptr = CS_end + 1
    !             end if
    !         end do
    !         row = row + blocks%blocks(i)%shape(2)
    !     end do
    !     CS%index_ptr(col+1) = ptr

    !     if (CS%index_ptr(CS%shape(2)+1) /= CS%nnz+1) then
    !         write(6,*) "Error in CSC from block assembly, index_ptr(N_cols+1) /= nnz + 1"
    !         stop
    !     end if
    ! end subroutine CSC_from_block_diag

    subroutine CS_block_store(this,loc)
        class(block_CS), intent(inout) :: this
        character(len=*), intent(in) :: loc

        integer :: i,j

        call this%compute_shape()

        open(unit=1,file=loc,action='write',form='unformatted')

        write(1) 'CSR'
        write(1) this%block_shape
        write(1) this%shape

        do j = 1, this%block_shape(2)
            do i = 1, this%block_shape(1)
                write(1) this%blocks(i,j)%shape
                write(1) this%blocks(i,j)%nnz

                if (this%blocks(i,j)%nnz > 0) then
                    write(1) this%blocks(i,j)%index_ptr
                    write(1) this%blocks(i,j)%indices
                    write(1) this%blocks(i,j)%data
                end if
            end do
        end do

        close(1)

    end subroutine CS_block_store

    subroutine CS_block_load(this, loc)
        class(block_CS), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        integer :: i,j
        integer :: nnz
        character(len=3) :: type
        integer, dimension(2) :: shape, block_shape

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        read(1) block_shape
        read(1) shape

        if (type == 'CSR') then
            call this%init(block_shape)
        else
            write(6,*) "Error! Only CSR blocks are implemented."
            write(6,*) type
            stop
        end if
        this%shape = shape

        do j = 1, this%block_shape(2)
            do i = 1, this%block_shape(1)
                read(1) shape
                read(1) nnz

                call this%blocks(i,j)%init(shape,nnz)

                if (this%blocks(i,j)%nnz > 0) then
                    read(1) this%blocks(i,j)%index_ptr
                    read(1) this%blocks(i,j)%indices
                    read(1) this%blocks(i,j)%data
                end if
            end do
        end do

        close(1)
    end subroutine CS_block_load

    subroutine CS_block_diag_store(this,loc)
        class(block_diag_CS), intent(inout) :: this
        character(len=*), intent(in) :: loc

        integer :: i

        call this%compute_shape()

        open(unit=1,file=loc,action='write',form='unformatted')

        write(1) 'CSR'
        write(1) this%block_shape
        write(1) this%shape

        do i = 1, this%block_shape(1)
            write(1) this%blocks(i)%shape
            write(1) this%blocks(i)%nnz

            if (this%blocks(i)%nnz > 0) then
                write(1) this%blocks(i)%index_ptr
                write(1) this%blocks(i)%indices
                write(1) this%blocks(i)%data
            end if
        end do

        close(1)

    end subroutine CS_block_diag_store

    subroutine CS_block_diag_load(this, loc)
        class(block_diag_CS), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        integer :: i
        integer :: nnz
        character(len=3) :: type
        integer, dimension(2) :: shape, block_shape

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        read(1) block_shape
        read(1) shape

        if (type == 'CSR') then
            call this%init(block_shape(1))
        else
            write(6,*) "Error! Blocks must be CSR!"
            write(6,*) type
            stop
        end if
        this%shape = shape

        do i = 1, this%block_shape(1)
            read(1) shape
            read(1) nnz

            call this%blocks(i)%init(shape,nnz)

            if (this%blocks(i)%nnz > 0) then
                read(1) this%blocks(i)%index_ptr
                read(1) this%blocks(i)%indices
                read(1) this%blocks(i)%data
            end if
        end do

        close(1)
    end subroutine CS_block_diag_load

    subroutine CS_block_scale(this, alpha)
        class(block_CS), intent(inout) :: this
        double complex, intent(in) ::  alpha

        integer i,j

        do j = 1,this%block_shape(2)
            do i = 1,this%block_shape(1)
                call this%blocks(i,j)%scale(alpha)
            end do
        end do
    end subroutine CS_block_scale

    subroutine CS_block_diag_scale(this, alpha)
        class(block_diag_CS), intent(inout) :: this
        double complex, intent(in) ::  alpha

        integer i

        do i = 1,this%block_shape(1)
            call this%blocks(i)%scale(alpha)
        end do
    end subroutine CS_block_diag_scale

    ! Function that returns the result of adding a block diagonal sparse matrix X to a block sparse matrix Y.
    ! The second block matrix is assumed to have empty blocks on the diagonal.
    ! The second matrix is scaled by alpha
    function CS_diag_X_P_A_block_Y(X,Y,alpha) result(res)
        type(block_diag_CS), intent(in) :: X
        type(block_CS), intent(in) :: Y
        double complex, intent(in) :: alpha
        type(block_CS) :: res

        integer :: i,j

        if (any(X%block_shape /= Y%block_shape)) then
            write(6,*) "Error in CS_diag_X_P_A_block_Y! Incompatible block shapes"
            stop
        end if

        call res%init(Y%block_shape)

        !$omp parallel do private(i) schedule(dynamic)
        do j = 1,Y%block_shape(2)
            res%blocks(j,j) = X%blocks(j)
            do i = 1,Y%block_shape(1)
                if (i == j) cycle
                res%blocks(i,j) = alpha*Y%blocks(i,j)
            end do
        end do
        !$omp end parallel do
    end function CS_diag_X_P_A_block_Y

    ! Function that returns the result of adding a block sparse matrix X to a block sparse matrix Y.
    ! It it assumed that X and Y do not have overlapping nonzero blocks
    ! The first matrix is scaled by alpha and the second by beta
    function CS_A_block_X_P_B_block_Y(X,alpha,Y,beta) result(res)
        type(block_CS), intent(in) :: X
        double complex, intent(in) :: alpha
        type(block_CS), intent(in) :: Y
        double complex, intent(in) :: beta
        type(block_CS) :: res

        integer :: i,j

        if (any(X%block_shape /= Y%block_shape)) then
            write(6,*) "Error in CS_A_block_X_P_B_block_Y! Incompatible block shapes"
            stop
        end if

        call res%init(Y%block_shape)

        !$omp parallel do private(i) schedule(dynamic)
        do j = 1,Y%block_shape(2)
            do i = 1,Y%block_shape(1)
                if (X%blocks(i,j)%nnz > 0) then
                    res%blocks(i,j) = alpha*X%blocks(i,j)
                else if (Y%blocks(i,j)%nnz > 0) then
                    res%blocks(i,j) = beta*Y%blocks(i,j)
                else
                    ! Block is zero so just do normal assignment
                    res%blocks(i,j) = X%blocks(i,j)
                end if
            end do
        end do
        !$omp end parallel do
    end function CS_A_block_X_P_B_block_Y

    subroutine CS_A_P_diag_X(X,alpha)
        type(block_diag_CS), intent(inout) :: X
        double complex, dimension(:), intent(in) :: alpha

        integer :: i

        if (X%block_shape(1) /= size(alpha)) then
            write(6,*) "Error in CS_A_P_diag_X! Incompatible block and alpha shapes"
            stop
        end if

        !$omp parallel do
        do i = 1,X%block_shape(1)
            call X%blocks(i)%shift(alpha(i))
        end do
        !$omp end parallel do
    end subroutine CS_A_P_diag_X

    subroutine CS_A_P_block_X(X,alpha)
        type(block_CS), intent(inout) :: X
        double complex, dimension(:), intent(in) :: alpha

        integer :: i

        if (X%block_shape(1) /= size(alpha)) then
            write(6,*) "Error in CS_A_P_block_X! Incompatible block and alpha shapes"
            stop
        end if

        !$omp parallel do
        do i = 1,X%block_shape(1)
            call X%blocks(i,i)%shift(alpha(i))
        end do
        !$omp end parallel do
    end subroutine CS_A_P_block_X

    subroutine CS_A_diag_B_P_diag_X(X,alpha,B)
        type(block_diag_CS), intent(inout) :: X
        double complex, dimension(:), intent(in) :: alpha
        type(block_diag_CS), intent(in) :: B

        integer :: i

        if (X%block_shape(1) /= size(alpha)) then
            write(6,*) "Error in CS_A_diag_B_P_diag_X! Incompatible block and alpha shapes"
            stop
        end if

        !!$omp parallel do
        do i = 1,X%block_shape(1)
            call X%blocks(i)%shift_B(alpha(i),B%blocks(i))
        end do
        !!$omp end parallel do
    end subroutine CS_A_diag_B_P_diag_X

    subroutine CS_A_diag_B_P_block_X(X,alpha,B)
        type(block_CS), intent(inout) :: X
        double complex, dimension(:), intent(in) :: alpha
        type(block_diag_CS), intent(in) :: B

        integer :: i

        if (X%block_shape(1) /= size(alpha)) then
            write(6,*) "Error in CS_A_diag_B_P_block_X! Incompatible block and alpha shapes"
            stop
        end if

        !!$omp parallel do
        do i = 1,X%block_shape(1)
            call X%blocks(i,i)%shift_B(alpha(i),B%blocks(i))
        end do
        !!$omp end parallel do
    end subroutine CS_A_diag_B_P_block_X

    subroutine CS_A_block_B_P_block_X(X,alpha,B)
        type(block_CS), intent(inout) :: X
        double complex, dimension(:,:), intent(in) :: alpha
        type(block_CS), intent(in) :: B

        integer :: i,j

        if (any(X%block_shape /= shape(alpha))) then
            write(6,*) "Error in CS_A_block_B_P_block_X! Incompatible block and alpha shapes"
            stop
        end if

        if (any(X%block_shape /= B%block_shape)) then
            write(6,*) "Error in CS_A_block_B_P_block_X! Incompatible block and alpha shapes"
            stop
        end if

        !!$omp parallel do private(i)
        do j = 1,X%block_shape(2)
            do i = 1,X%block_shape(1)
                call X%blocks(i,j)%shift_B(alpha(i,j),B%blocks(i,j))
            end do
        end do
        !!$omp end parallel do
    end subroutine CS_A_block_B_P_block_X
end module block_tools