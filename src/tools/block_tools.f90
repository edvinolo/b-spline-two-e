module block_tools
    use sparse_array_tools, only: CS_matrix,CSR_matrix,CSC_matrix
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
        class(CS_matrix), dimension(:,:), allocatable :: blocks
    contains
        procedure :: init => init_block_CS
        procedure :: deall => deallocate_block_CS
        procedure :: compute_shape => compute_shape_block_CS
        procedure :: to_CS => CS_from_block
        procedure :: store => CS_block_store
        procedure :: load => CS_block_load
    end type block_CS

    type, public :: block_diag_CS
        integer, dimension(2) :: block_shape
        integer, dimension(2) :: shape
        class(CS_matrix), dimension(:), allocatable :: blocks
    contains
        procedure :: init => init_block_diag_CS
        procedure :: deall => deallocate_block_diag_CS
        procedure :: compute_shape => compute_shape_block_diag_CS
        procedure :: to_CS => CS_from_block_diag
        procedure :: store => CS_block_diag_store
        procedure :: load => CS_block_diag_load
    end type  block_diag_CS

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

    subroutine init_block_CS(this,shape_in,block_type)
        class(block_CS), intent(inout) :: this
        integer, dimension(2), intent(in) :: shape_in
        class(CS_matrix), intent(in) :: block_type

        this%block_shape = shape_in

        select type(block_type)
        type is (CSR_matrix)
            allocate(CSR_matrix :: this%blocks(shape_in(1),shape_in(2)))
        type is (CSC_matrix)
            allocate(CSC_matrix :: this%blocks(shape_in(1),shape_in(2)))
        class default
            write(6,*) "block_type must be CSR_matrix or CSC_matrix"
        end select
    end subroutine init_block_CS

    subroutine init_block_diag_CS(this,n,block_type)
        class(block_diag_CS), intent(inout) :: this
        integer, intent(in) :: n
        class(CS_matrix), intent(in) :: block_type

        this%block_shape = [n,n]

        select type(block_type)
        type is (CSR_matrix)
            allocate(CSR_matrix :: this%blocks(n))
        type is (CSC_matrix)
            allocate(CSC_matrix :: this%blocks(n))
        class default
            write(6,*) "block_type must be CSR_matrix or CSC_matrix"
        end select
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
            this%shape = this%shape + this%blocks(i,i)%shape
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
        class(CS_matrix), intent(inout) :: CS
        logical, intent(in) :: deall_block

        integer :: i, j, nnz

        ! check that the individual blocks and CS are of the same type
        if (.not.same_type_as(CS,blocks%blocks(1,1))) then
            write(6,*) "CS and blocks must be of the same type"
            stop
        end if

        nnz = 0
        do j = 1,blocks%block_shape(2)
            do i = 1,blocks%block_shape(1)
                if (blocks%blocks(i,j)%arrays_allocated()) then
                    nnz = nnz + blocks%blocks(i,j)%nnz
                end if
            end do
        end do

        call blocks%compute_shape()
        call CS%init(blocks%shape,nnz)

        select type(CS)
        type is(CSR_matrix)
            call CSR_from_block(blocks,CS)
        type is(CSC_matrix)
            call CSC_from_block(blocks,CS)
        class default
            write(6,*) "blocks must be CSR_matrix or CSC_matrix"
            stop
        end select

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
            stop
        end if

    end subroutine CSR_from_block

    subroutine CSC_from_block(blocks,CS)
        class(block_CS), intent(inout) :: blocks
        type(CSC_matrix), intent(inout) :: CS

        integer :: i, j, col_i
        integer :: ptr, row, col, start, end, CS_start, CS_end

        ptr = 1
        col = 0

        do i = 1,blocks%block_shape(2)
            do col_i = 1, blocks%blocks(1,i)%shape(2)
                col = col +1
                row = 0
                CS%index_ptr(col) = ptr
                do j = 1,blocks%block_shape(1)
                    if (blocks%blocks(j,i)%arrays_allocated()) then
                        start = blocks%blocks(j,i)%index_ptr(col_i)
                        end  = blocks%blocks(j,i)%index_ptr(col_i+1)-1
                        CS_start = ptr
                        CS_end = ptr + (end-start)
                        CS%indices(CS_start:CS_end) = col + blocks%blocks(j,i)%indices(start:end)
                        CS%data(CS_start:CS_end) = blocks%blocks(j,i)%data(start:end)
                        ptr = CS_end + 1
                    end if
                    row = row + blocks%blocks(j,i)%shape(1)
                end do
            end do
        end do
        CS%index_ptr(col+1) = ptr

        if (CS%index_ptr(CS%shape(1)+1) /= CS%nnz+1) then
            write(6,*) "Error in CSC from block assembly, index_ptr(N_cols+1) /= nnz + 1"
            stop
        end if
    end subroutine CSC_from_block

    subroutine CS_from_block_diag(blocks,CS,deall_block)
        class(block_diag_CS), intent(inout) :: blocks
        class(CS_matrix), intent(inout) :: CS
        logical, intent(in) :: deall_block

        integer :: i, j, nnz, row_i
        integer :: ptr, row, col, start, end, CS_start, CS_end

        ! check that the individual blocks and CS are of the same type
        if (.not.same_type_as(CS,blocks%blocks(1))) then
            write(6,*) "CS and blocks must be of the same type"
            stop
        end if

        nnz = 0
        do i = 1,blocks%block_shape(1)
            if (blocks%blocks(i)%arrays_allocated()) then
                nnz = nnz + blocks%blocks(i)%nnz
            end if
        end do

        call blocks%compute_shape()
        call CS%init(blocks%shape,nnz)

        select type(CS)
        type is(CSR_matrix)
            call CSR_from_block_diag(blocks,CS)
        type is (CSC_matrix)
            call CSC_from_block_diag(blocks,CS)
        class default
            write(6,*) "blocks must be CSR_matrix or CSC_matrix"
            stop
        end select


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

        do i = 1,blocks%block_shape(1)
            do row_i = 1, blocks%blocks(i)%shape(1)
                col = 0
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

                col = col + blocks%blocks(i)%shape(2)
            end do
        end do
        CS%index_ptr(row+1) = ptr

        if (CS%index_ptr(CS%shape(1)+1) /= CS%nnz+1) then
            write(6,*) "Error in CSR from block assembly, index_ptr(N_rows+1) /= nnz + 1"
            stop
        end if
    end subroutine CSR_from_block_diag

    subroutine CSC_from_block_diag(blocks,CS)
        type(block_diag_CS), intent(inout) :: blocks
        type(CSC_matrix), intent(inout) :: CS

        integer :: i, col_i
        integer :: ptr, row, col, start, end, CS_start, CS_end

        ptr = 1
        col = 0

        do i = 1,blocks%block_shape(1)
            do col_i = 1, blocks%blocks(i)%shape(1)
                row = 0
                col = col + 1
                CS%index_ptr(col) = ptr

                if (blocks%blocks(i)%arrays_allocated()) then
                    start = blocks%blocks(i)%index_ptr(col_i)
                    end  = blocks%blocks(i)%index_ptr(col_i+1)-1
                    CS_start = ptr
                    CS_end = ptr + (end-start)
                    CS%indices(CS_start:CS_end) = col + blocks%blocks(i)%indices(start:end)
                    CS%data(CS_start:CS_end) = blocks%blocks(i)%data(start:end)
                    ptr = CS_end + 1
                end if

                row = row + blocks%blocks(i)%shape(2)
            end do
        end do
        CS%index_ptr(col+1) = ptr

        if (CS%index_ptr(CS%shape(2)+1) /= CS%nnz+1) then
            write(6,*) "Error in CSC from block assembly, index_ptr(N_cols+1) /= nnz + 1"
            stop
        end if
    end subroutine CSC_from_block_diag

    subroutine CS_block_store(this,loc)
        class(block_CS), intent(inout) :: this
        character(len=*), intent(in) :: loc

        integer :: i,j

        call this%compute_shape()

        open(unit=1,file=loc,action='write',form='unformatted')

        associate (element => this%blocks(1,1))
            select type (element)
            type is (CSR_matrix)
                write(1) 'CSR'
            type is (CSC_matrix)
                write(1) 'CSC'
            class default
                write(6,*) 'Blocks must be CSR_matrix or CSC_matrix'
                close(1)
                stop
            end select
        end associate

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
        type(CSR_matrix) :: CSR
        type(CSC_matrix) :: CSC

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        read(1) block_shape
        read(1) shape

        if (type == 'CSR') call this%init(block_shape,CSR)
        if (type == 'CSC') call this%init(block_shape,CSC)
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

        associate (element => this%blocks(1))
            select type (element)
            type is (CSR_matrix)
                write(1) 'CSR'
            type is (CSC_matrix)
                write(1) 'CSC'
            class default
                write(6,*) 'Blocks must be CSR_matrix or CSC_matrix'
                close(1)
                stop
            end select
        end associate

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
        type(CSR_matrix) :: CSR
        type(CSC_matrix) :: CSC

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        read(1) block_shape
        read(1) shape

        if (type == 'CSR') call this%init(block_shape(1),CSR)
        if (type == 'CSC') call this%init(block_shape(1),CSC)
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
end module block_tools