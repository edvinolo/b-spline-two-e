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

    type, abstract, public :: CS_matrix
        integer, dimension(2) :: shape
        integer :: nnz
        integer, dimension(:), allocatable :: indices
        integer, dimension(:), allocatable :: index_ptr
        double complex, dimension(:), allocatable :: data
    contains
        procedure :: init => init_CS
        procedure :: enter => enter_CS
        procedure :: convert => convert_CS
        procedure :: arrays_allocated
        procedure :: deall => deall_CS
        procedure :: store_CS
        procedure :: load_CS
        procedure(store), deferred :: store
        procedure(load), deferred :: load
        procedure(ptr_size), deferred :: ptr_size
        procedure(get_dense), deferred :: get_dense
    end type CS_matrix

    !Deferred procedures need to provide an interface
    abstract interface
        subroutine get_dense(this,matrix)
            import :: CS_matrix
            class(CS_matrix), intent(in) :: this
            double complex, dimension(:,:), allocatable, intent(inout) :: matrix
        end subroutine get_Dense

        function ptr_size (this) result(out)
            import :: CS_matrix
            class(CS_matrix), intent(in) :: this
            integer :: out
        end function ptr_size

        subroutine store(this,loc)
            import :: CS_matrix
            class(CS_matrix), intent(in) :: this
            character(len=*), intent(in) ::  loc
        end subroutine store

        subroutine load(this,loc)
            import :: CS_matrix
            class(CS_matrix), intent(inout) :: this
            character(len=*), intent(in) ::  loc
        end subroutine load
    end interface

    !The child types of CS_matrix. They must overwrite the deferred procedures
    !in the parent type.
    type, extends(CS_matrix), public :: CSR_matrix
    !Stores elements row by row, with pointers to beginning of each row.
    contains
        procedure :: store => CSR_store
        procedure :: load => CSR_load
        procedure :: ptr_size => CSR_ptr_size
        procedure :: get_dense => CSR_get_dense
    end type CSR_matrix

    type, extends(CS_matrix), public :: CSC_matrix
    !Stores elements column by column, with pointers to beginning of each column.
    contains
        procedure :: store => CSC_store
        procedure :: load => CSC_load
        procedure :: ptr_size => CSC_ptr_size
        procedure :: get_dense => CSC_get_dense
    end type CSC_matrix
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

    subroutine init_CS(this,shape,nnz)
        class(CS_matrix), intent(inout) :: this
        integer, dimension(2), intent(in) ::  shape
        integer, intent(in) :: nnz

        this%shape = shape
        this%nnz = nnz

        if (nnz>0) then
            allocate(this%indices(nnz),this%index_ptr(this%ptr_size()),source = 0)
            allocate(this%data(nnz),source = dcmplx(0.d0,0.d0))
        end if
    end subroutine init_CS

    !> Subroutine that sets the data array components
    !! using arrays provided by the user.
    !!
    !! @param indices
    !! @param index_ptr
    !! @param elements
    subroutine enter_CS(this,indices,index_ptr,elements)
        class(CS_matrix), intent(inout) :: this
        integer, dimension(:), intent(in) :: indices
        integer, dimension(:), intent(in) :: index_ptr
        double complex, dimension(:), intent(in) :: elements

        if (.not.(this%arrays_allocated())) then
            write(6,*) 'Arrays of CS matrix have not been allocated, please setup before entering elements'
            stop
        end if

        if ((lbound(this%indices,dim = 1).ne.lbound(indices,dim = 1)).or.&
            (ubound(this%indices,dim = 1).ne.ubound(indices,dim = 1))) then
            write(6,*) 'Error: indices array has different bounds than that in CS object'
            stop
        end if

        if ((lbound(this%index_ptr, dim = 1).ne.lbound(index_ptr, dim = 1)).or.&
            (ubound(this%index_ptr, dim = 1).ne.ubound(index_ptr, dim = 1))) then
            write(6,*) 'Error: index_ptr array has different bounds than that in CS object'
            stop
        end if

        if ((lbound(this%data, dim = 1).ne.lbound(elements, dim  = 1)).or.&
            (ubound(this%data, dim = 1).ne.ubound(elements, dim = 1))) then
            write(6,*) 'Error: data array has different bounds than that in CS object'
            stop
        end if

        this%indices = indices
        this%index_ptr = index_ptr
        this%data = elements
    end subroutine enter_CS

    !> Subroutine for converting between the two CS types
    !!
    !! @param A: the new representation
    subroutine convert_CS(this,A)
        class(CS_matrix), intent(in) :: this
        class(CS_matrix), intent(out) :: A

        integer :: i,j,k,ptr

        !The way things are setup A is canonical (sorted) by default
        call A%init(this%shape,this%nnz)

        !Count nnz for each index along the new dimension
        !put it in the next element up
        A%index_ptr = 0
        do i = lbound(this%indices, dim  = 1), ubound(this%indices, dim = 1)
            A%index_ptr(this%indices(i)+1) = A%index_ptr(this%indices(i)+1) + 1
        end do

        !Compute the new index pointers by adding the length of previous row/col
        A%index_ptr(1) = 1
        do i = lbound(A%index_ptr, dim = 1), ubound(A%index_ptr, dim = 1) - 1
            A%index_ptr(i+1) = A%index_ptr(i) + A%index_ptr(i+1)
        end do

        !Fill indices and elements of new CS_matrix.
        !index_ptr is updated to the next element in the row/col.
        !This loop automatically puts A in canonical order, since i is increasing.
        do i  = lbound(this%index_ptr, dim = 1), ubound(this%index_ptr,dim = 1) - 1
            do k = this%index_ptr(i),this%index_ptr(i+1)-1
                j = this%indices(k)
                ptr = A%index_ptr(j)
                A%data(ptr) = this%data(k)
                A%indices(ptr) = i
                A%index_ptr(j) = ptr + 1
            end do
        end do

        !Reset pointer after the filling
        do i = ubound(A%index_ptr,dim = 1) - 1, lbound(A%index_ptr,dim = 1),-1
            A%index_ptr(i+1) = A%index_ptr(i)
        end do
        A%index_ptr(1) = 1
    end subroutine convert_CS

    function arrays_allocated(this) result(retval)
        class(CS_matrix), intent(in) :: this
        logical :: retval

        if (allocated(this%index_ptr).and.allocated(this%indices).and.allocated(this%data)) then
            retval = .true.
        else
            retval = .false.
        end if
    end function arrays_allocated

    subroutine deall_CS(this)
        class(CS_matrix), intent(inout) :: this

        if (allocated(this%index_ptr)) deallocate(this%index_ptr)
        if (allocated(this%indices)) deallocate(this%indices)
        if (allocated(this%data)) deallocate(this%data)
    end subroutine deall_CS

    subroutine store_CS(this,loc,type)
        class(CS_matrix), intent(in) :: this
        character(len=*), intent(in) :: loc
        character(len=3), intent(in) :: type

        open(unit=1,file=loc,action='write',form='unformatted')
        write(1) type
        write(1) this%shape
        write(1) this%nnz
        write(1) this%index_ptr
        write(1) this%indices
        write(1) this%data
        close(1)
    end subroutine store_CS

    subroutine CSR_store(this,loc)
        class(CSR_matrix), intent(in) :: this
        character(len=*), intent(in) :: loc

        call this%store_CS(loc,'CSR')
    end subroutine CSR_store

    subroutine CSC_store(this,loc)
        class(CSC_matrix), intent(in) :: this
        character(len=*), intent(in) :: loc

        call this%store_CS(loc,'CSC')
    end subroutine CSC_store

    subroutine load_CS(this,loc)
        class(CS_matrix), intent(inout) :: this
        character(len=*), intent(in) :: loc

        character(len=3) :: type

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        read(1) this%shape
        read(1) this%nnz

        call this%deall()
        allocate(this%index_ptr(this%ptr_size()),this%indices(this%nnz),this%data(this%nnz))
        read(1) this%index_ptr
        read(1) this%indices
        read(1) this%data

        close(1)
    end subroutine load_CS

    subroutine CSR_load(this,loc)
        class(CSR_matrix), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        character(len=3) :: type

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        close(1)

        if (type/='CSR') then
            write(6,*) "Wrong matrix type: ", type, " when trying to load CSR."
            stop
        end if

        call this%load_CS(loc)
    end subroutine CSR_load

    subroutine CSC_load(this,loc)
        class(CSC_matrix), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        character(len=3) :: type

        open(unit=1,file=loc,action='read',form='unformatted')
        read(1) type
        close(1)

        if (type/='CSC') then
            write(6,*) "Wrong matrix type: ", type, " when trying to load CSC."
            stop
        end if

        call this%load_CS(loc)
    end subroutine CSC_load

    function CSR_ptr_size(this) result(res)
        class(CSR_matrix), intent(in) :: this
        integer :: res

        res = this%shape(1) + 1
    end function CSR_ptr_size

    function CSC_ptr_size(this) result(res)
        class(CSC_matrix), intent(in) :: this
        integer :: res

        res = this%shape(2) + 1
    end function CSC_ptr_size

    !> Subroutines that generates the dense representation for the underlying matrix.
    !!
    !! @param matrix: array that holds the dense matrix representation,
    !! should be unallocated on entry.
    subroutine CSR_get_dense(this,matrix)
        class(CSR_matrix), intent(in) :: this
        double complex, dimension(:,:), allocatable, intent(inout) :: matrix

        integer :: i,row_start,row_end

        allocate(matrix(this%shape(1),this%shape(2)),source = dcmplx(0.d0,0.d0))

        do i  = lbound(this%index_ptr,dim = 1), ubound(this%index_ptr, dim = 1) - 1
            row_start = this%index_ptr(i)
            row_end = this%index_ptr(i+1)-1
            if (row_start <= row_end) then !True if row has elements
                matrix(i,this%indices(row_start:row_end)) = this%data(row_start:row_end)
            end if
        end do
    end subroutine CSR_get_dense

    subroutine CSC_get_dense(this,matrix)
        class(CSC_matrix), intent(in) :: this
        double complex, dimension(:,:), allocatable, intent(inout) :: matrix

        integer :: i,col_start,col_end

        allocate(matrix(this%shape(1),this%shape(2)),source = dcmplx(0.d0,0.d0))

        do i  = lbound(this%index_ptr,dim = 1), ubound(this%index_ptr, dim = 1) - 1
            col_start = this%index_ptr(i)
            col_end = this%index_ptr(i+1)-1
            if (col_start <= col_end) then !True if col has elements
                matrix(this%indices(col_start:col_end),i) = this%data(col_start:col_end)
            end if
        end do
    end subroutine CSC_get_dense
end module sparse_array_tools