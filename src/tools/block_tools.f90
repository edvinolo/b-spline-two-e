module block_tools
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
end module block_tools