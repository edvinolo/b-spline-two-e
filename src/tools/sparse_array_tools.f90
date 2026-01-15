module sparse_array_tools
    use kind_tools
    use bspline_tools
    use stdlib_sorting, only: sort_index
    use stdlib_hashmaps, only: open_hashmap_type
    use omp_lib, only: omp_get_wtime
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

    type, public :: ptr_type
        integer :: ptr
    end type ptr_type

    type, public :: Nd_DOK
        integer :: N
        integer :: N_val
        integer :: N_elements
        type(ptr_type) :: entry
        type(open_hashmap_type) :: map
        real(dp), allocatable :: data(:)
        integer(int8), allocatable :: indices(:)
    contains
        procedure :: init => init_Nd_DOK
        procedure :: set_val => set_val_Nd_DOK
        procedure :: get_val => get_val_Nd_DOK
    end type Nd_DOK

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
        procedure :: transp => transpose_CS
        procedure :: arrays_allocated
        procedure :: deall => deall_CS
        procedure :: store_CS
        procedure :: load_CS
        procedure :: shift => shift_CS
        procedure :: shift_B => shift_B_CS
        procedure :: scale => scale_CS
        procedure, pass(B) :: assign_CS
        generic :: assignment(=) => assign_CS
        generic :: operator(*) => scalar_mult
        procedure(scalar_mult), pass(X), deferred :: scalar_mult
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

        function scalar_mult(alpha,X) result(res)
            import :: CS_matrix
            double complex, intent(in) :: alpha
            class(CS_matrix), intent(in) :: X
            class(CS_matrix), allocatable :: res
        end function scalar_mult
    end interface

    !The child types of CS_matrix. They must overwrite the deferred procedures
    !in the parent type.
    type, extends(CS_matrix), public :: CSR_matrix
    !Stores elements row by row, with pointers to beginning of each row.
    contains
        procedure, pass(X) :: scalar_mult => scalar_mult_CSR
        procedure :: store => CSR_store
        procedure :: load => CSR_load
        procedure :: ptr_size => CSR_ptr_size
        procedure :: get_dense => CSR_get_dense
    end type CSR_matrix

    type, extends(CS_matrix), public :: CSC_matrix
    !Stores elements column by column, with pointers to beginning of each column.
    contains
        procedure, pass(X) :: scalar_mult => scalar_mult_CSC
        procedure :: store => CSC_store
        procedure :: load => CSC_load
        procedure :: ptr_size => CSC_ptr_size
        procedure :: get_dense => CSC_get_dense
    end type CSC_matrix

    ! Abstract interface for matvecs
    abstract interface
        subroutine mv(A,x,y)
            import :: dp
            import :: CSR_matrix
            type(CSR_matrix), intent(in) :: a
            complex(dp), intent(in) :: x(:)
            complex(dp), intent(out) :: y(:)
        end subroutine mv
    end interface
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

    pure function count_nnz_R_k(b_splines) result(res)
        type(b_spline), intent(in) :: b_splines
        integer :: res

        integer :: i,j,i_p,j_p

        res = 0

        do j_p = 1, b_splines%n_b
            do j = 1, b_splines%n_b
                if (abs(j-j_p)>=b_splines%k) cycle
                do i_p = 1, b_splines%n_b
                    do i = 1,b_splines%n_b
                        if (abs(i-i_p) >= b_splines%k) cycle
                        res = res + 1
                    end do
                end do
            end do
        end do
    end function count_nnz_R_k

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
        ! integer :: i,j,nnz

        write(6,*)
        write(6,*) 'Constructing R_K array...'
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
        write(6,*) 'Done!'
        ! write(6,*) R(:,8,9,10,11)

        ! nnz = 0
        ! do i = 1,n_b
        !     do j = 1,n_b
        !         do n = 1,n_b
        !             do m = 1,n_b
        !                 do k = 0,max_k
        !                     if (abs(R(k,m,n,j,i)) > 1.0d-12) nnz = nnz +1
        !                 end do
        !             end do
        !         end do
        !     end do
        ! end do

        ! write(6,*) 'R_k density: ', real(nnz,kind=8)/(n_b**4)

    end subroutine compute_R_k

    subroutine compute_R_K_map(r_d_k,r_k,r_m_k,b_splines,max_k,R)
        type(sparse_6d), intent(in) :: r_d_k
        type(sparse_4d), intent(in) :: r_k
        type(sparse_4d), intent(in) :: r_m_k
        type(b_spline), intent(in) :: b_splines
        integer, intent(in) :: max_k
        type(Nd_DOK), intent(inout) :: R

        integer :: n,m,nnz
        double precision, allocatable :: val(:)
        double precision :: t1,t2

        write(6,*)
        write(6,*) 'Constructing R_K DOK array...'
        t1 = omp_get_wtime()

        allocate(val(0:max_k))
        nnz = count_nnz_R_k(b_splines)
        call R%init(4,max_k+1,nnz) !Four indices, and max_k+1 values

        do n=1,r_k%nnz
            do m=1,r_k%nnz
                if (r_k%iv(n)<r_k%iv(m)) then
                    val = r_k%data(n,:)*r_m_k%data(m,:)
                    call R%set_val([r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)],val)
                else if (r_k%iv(n)>r_k%iv(m)) then
                    val = r_k%data(m,:)*r_m_k%data(n,:)
                    call R%set_val([r_k%i(n),r_k%i(m),r_k%j(n),r_k%j(m)],val)
                end if
            end do
        end do

        do n=1,r_d_k%nnz
            val = r_d_k%data(n,:)
            call R%set_val([r_d_k%i(n),r_d_k%j(n),r_d_k%i_p(n),r_d_k%j_p(n)],val)
            call R%set_val([r_d_k%j(n),r_d_k%i(n),r_d_k%j_p(n),r_d_k%i_p(n)],val)
        end do
        t2 = omp_get_wtime()
        write(6,*) 'Done!, Time for construction (s): ', t2-t1
        ! call R%get_val([8,9,10,11],val)
        ! write(6,*) val
    end subroutine compute_R_K_map

    subroutine init_Nd_DOK(this,N,N_val,nnz)
        class(Nd_DOK), intent(out) :: this
        integer, intent(in) :: N
        integer, intent(in) :: N_val
        integer, intent(in) :: nnz

        call this%map%init()
        this%N = N
        this%N_val = N_val
        this%N_elements = nnz
        this%entry%ptr = 1
        allocate(this%data(N_val*this%N_elements))
    end subroutine init_Nd_DOK

    subroutine set_val_Nd_DOK(this,indices,val)
        class(Nd_DOK), intent(inout) :: this
        integer, intent(in) :: indices(this%N)
        real(dp), intent(in) :: val(:)

        logical(int32) :: exists
        class(*), allocatable :: temp

        this%indices = transfer(indices,[0_int8])

        call this%map%key_test(this%indices,exists)

        if (exists) then
            call this%map%get_other_data(this%indices,temp)
            select type (temp)
            type is (ptr_type)
                this%data(temp%ptr:temp%ptr+this%N_val-1) = this%data(temp%ptr:temp%ptr+this%N_val-1) + val
            class default
                error stop "Wrong dynamic type of ptr!"
            end select
        else
            this%data(this%entry%ptr:this%entry%ptr+this%N_val-1) = val
            call this%map%map_entry(this%indices,this%entry)
            this%entry%ptr = this%entry%ptr + this%N_val
        end if
    end subroutine set_val_Nd_DOK

    subroutine get_val_Nd_DOK(this,indices,values)
        class(Nd_DOK), intent(inout) :: this
        integer, intent(in) :: indices(this%N)
        real(dp), intent(out) :: values(:)

        integer(int8), allocatable :: key(:)
        class(*), allocatable :: temp

        key = transfer(indices,[0_int8])

        !$omp critical(hashmap)
        call this%map%get_other_data(key,temp)
        !$omp end critical(hashmap)
        select type (temp)
        type is (ptr_type)
            values = this%data(temp%ptr:temp%ptr+this%N_val-1)
        class default
            error stop 'Wrong dynamic type of ptr'
        end select
    end subroutine get_val_Nd_DOK

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
    subroutine convert_CS(this,A,transp)
        class(CS_matrix), intent(in) :: this
        class(CS_matrix), intent(out) :: A
        logical, optional, intent(in) :: transp

        integer :: i,j,k,ptr

        !The way things are setup A is canonical (sorted) by default
        if (present(transp).and.transp) then
            call A%init(this%shape([2,1]),this%nnz)
        else
            call A%init(this%shape,this%nnz)
        end if

        if (this%nnz == 0) return !Don't do anything if nnz = 0

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

    subroutine transpose_CS(this,A)
        class(CS_matrix), intent(in) :: this
        class(CS_matrix), intent(out) :: A

        if (.not.same_type_as(this,A)) then
            write(6,*) "Error in transpose_CS! this and A should have same dynamic type!"
        end if

        call this%convert(A,transp=.true.)

    end subroutine transpose_CS

    pure function arrays_allocated(this) result(retval)
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
        this%nnz = 0
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

    ! Shift the diagonal of the matrix in place. Right now only works if all diagonal elements are nonzero.
    subroutine shift_CS(this,shift)
        class(CS_matrix), intent(inout) :: this
        double complex, intent(in) :: shift

        integer, dimension(:), allocatable :: diag
        integer :: i,ptr,n_diag, diag_count

        n_diag = min(this%shape(1),this%shape(2))

        allocate(diag(n_diag))

        diag_count = 0
        do i = 1,n_diag
            do ptr = this%index_ptr(i),this%index_ptr(i+1)-1
                if (i == this%indices(ptr)) then
                    diag(i) = ptr
                    diag_count = diag_count + 1
                    exit
                end if
            end do
        end do

        if (diag_count /= n_diag) then
            write(6,*) "Error in shift_CS diag_count /= n_diag! Changes in sparsity structure is not yet implemented."
            stop
        end if

        this%data(diag) = this%data(diag) + shift
    end subroutine shift_CS

    ! Perform A = A + shift*B, where B's sparsity pattern is a subset of A's.
    subroutine shift_B_CS(this,shift,B)
        class(CS_matrix), intent(inout) :: this
        double complex, intent(in) :: shift
        class(CS_matrix), intent(in) :: B

        integer :: i, ptr_a, ptr_b

        if (.not.same_type_as(this,B)) then
            write(6,*) "Error in shift_B_CS! B must have same type as this"
            stop
        end if

        if (any(this%shape /= B%shape)) then
            write(6,*) "Error in shift_B_CS! this and B must have same shape"
            stop
        end if

        if (this%nnz < B%nnz) then
            write(6,*) "Error in shift_B_CS! this%nnz must be >= B%nnz"
            stop
        end if

        !!$omp parallel do private(ptr_a,ptr_b)
        do i = 1,this%ptr_size()-1
            ptr_b = B%index_ptr(i)
            do ptr_a = this%index_ptr(i),this%index_ptr(i+1)-1
                if (this%indices(ptr_a) == B%indices(ptr_b)) then
                    this%data(ptr_a) = this%data(ptr_a) + shift*B%data(ptr_b)
                    ptr_b = ptr_b + 1
                end if
            end do
        end do
        !!$omp end parallel
    end subroutine shift_B_CS

    subroutine scale_CS(this, alpha)
        class(CS_matrix), intent(inout) :: this
        double complex, intent(in) ::  alpha

        integer :: i

        if (this%nnz>0) then
            !$omp parallel do
            do i = 1,this%nnz
                this%data(i) = alpha*this%data(i)
            end do
            !$omp end parallel do
        end if
    end subroutine scale_CS

    elemental subroutine assign_CS(A,B)
        class(CS_matrix), intent(inout) :: A
        class(CS_matrix), intent(in) :: B

        if (same_type_as(A,B)) then
            A%shape = B%shape
            A%nnz = B%nnz

            if (B%arrays_allocated()) then
                A%index_ptr = B%index_ptr
                A%indices = B%indices
                A%data = B%data
            end if
        ! else
            ! error stop
        end if
    end subroutine assign_CS

    function scalar_mult_CSR(alpha,X) result(res)
        double complex, intent(in) :: alpha
        class(CSR_matrix), intent(in) :: X
        class(CS_matrix), allocatable :: res

        allocate(CSR_matrix :: res)

        res%shape = X%shape
        res%nnz = X%nnz

        if (X%arrays_allocated()) then
            res%index_ptr = X%index_ptr
            res%indices = X%indices
            res%data = alpha*X%data
        end if
    end function scalar_mult_CSR

    function scalar_mult_CSC(alpha,X) result(res)
        double complex, intent(in) :: alpha
        class(CSC_matrix), intent(in) :: X
        class(CS_matrix), allocatable :: res

        allocate(CSC_matrix :: res)

        res%shape = X%shape
        res%nnz = X%nnz

        if (X%arrays_allocated()) then
            res%index_ptr = X%index_ptr
            res%indices = X%indices
            res%data = alpha*X%data
        end if
    end function scalar_mult_CSC

    subroutine CSR_mv(A,x,y)
        type(CSR_matrix), intent(in) :: A
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)

        integer :: i,j

        if ((A%shape(2)/=size(x)).or.(A%shape(1)/=size(y))) then
            write(6,*) "Error! Incompatible array sizes in CSR_mv!"
            write(6,*) "A shape: ", A%shape
            write(6,*) "x shape: ", size(x)
            write(6,*) "y shape: ", size(y)
            stop
        end if

        if (A%nnz == 0) then
            y = 0_dp
            return
        end if

        !$omp parallel do private(j)
        do i = 1,A%shape(1)
            y(i) = 0_dp
            do j = A%index_ptr(i),A%index_ptr(i+1)-1
                y(i) = y(i) + A%data(j)*x(A%indices(j))
            end do
        end do
        !$omp end parallel do
    end subroutine CSR_mv

    subroutine CSR_T_mv(A,x,y)
        type(CSR_matrix), intent(in) :: A
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)

        integer :: i,j

        if ((A%shape(1)/=size(x)).or.(A%shape(2)/=size(y))) then
            write(6,*) "Error! Incompatible array sizes in CSR_T_mv!"
            write(6,*) "A shape: ", A%shape
            write(6,*) "x shape: ", size(x)
            write(6,*) "y shape: ", size(y)
            stop
        end if

        if (A%nnz == 0) then
            y = 0_dp
            return
        end if

        !$omp parallel private(j)
        !$omp do
        do i = 1,A%shape(1)
            y(i) = 0_dp
        end do
        !$omp end do

        !$omp do reduction(+:y)
        do i = 1,A%shape(1)
            do j = A%index_ptr(i),A%index_ptr(i+1)-1
                y(A%indices(j)) = y(A%indices(j)) + A%data(j)*x(i)
            end do
        end do
        !$omp end do
        !$omp end parallel
    end subroutine CSR_T_mv

    subroutine CSR_mv_sym(A,x,y)
        type(CSR_matrix), intent(in) :: A
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)

        integer :: i,j

        if ((A%shape(2)/=size(x)).or.(A%shape(1)/=size(y))) then
            write(6,*) "Error! Incompatible array sizes in CSR_mv!"
            write(6,*) "A shape: ", A%shape
            write(6,*) "x shape: ", size(x)
            write(6,*) "y shape: ", size(y)
            stop
        end if

        if (A%nnz == 0) then
            y = 0.0_dp
            return
        end if

        call mkl_zcsrsymv('U',A%shape(1),A%data,A%index_ptr,A%indices,x,y)

        ! This is correct, but for large problems the omp reduction clause runs into stack problems
        ! !$omp parallel private(j)
        ! !$omp do
        ! do i = 1,A%shape(1)
        !     y(i) = 0.0_dp
        !     do j = A%index_ptr(i),A%index_ptr(i+1)-1
        !         y(i) = y(i) + A%data(j)*x(A%indices(j))
        !     end do
        ! end do
        ! !$omp end do

        ! !$omp do reduction(+:y)
        ! do i = 1,A%shape(1)
        !     do j = A%index_ptr(i)+1,A%index_ptr(i+1)-1
        !         y(A%indices(j)) = y(A%indices(j)) + A%data(j)*x(i)
        !     end do
        ! end do
        ! !$omp end do
        ! !$omp end parallel
    end subroutine CSR_mv_sym

    subroutine CSR_dsymv(n,a,ia,ja,alpha,x,beta,y)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(:)
        integer, intent(in) :: ia(:)
        integer, intent(in) :: ja(:)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: beta
        real(dp), intent(inout) :: y(:)

        integer :: i,j

        ! if ((A%shape(2)/=size(x)).or.(A%shape(1)/=size(y))) then
        !     write(6,*) "Error! Incompatible array sizes in CSR_mv!"
        !     write(6,*) "A shape: ", A%shape
        !     write(6,*) "x shape: ", size(x)
        !     write(6,*) "y shape: ", size(y)
        !     stop
        ! end if

        ! y=beta*y
        !$omp parallel private(j)
        !$omp do
        do i = 1,n
            y(i) = beta*y(i)
            do j = ia(i),ia(i+1)-1
                y(i) = y(i) + alpha*a(j)*x(ja(j))
            end do
        end do
        !$omp end do

        !$omp do reduction(+:y)
        do i = 1,n
            do j = ia(i)+1,ia(i+1)-1
                y(ja(j)) = y(ja(j)) + alpha*a(j)*x(i)
            end do
        end do
        !$omp end do
        !$omp end parallel

    end subroutine CSR_dsymv

    subroutine check_diag_dom(A)
        class(CS_matrix), intent(in) :: A

        integer :: i,j
        real(dp) :: sum, diag
        logical, allocatable :: not_diag_dom(:)

        allocate(not_diag_dom(A%shape(1)))
        !$omp parallel do private(j,diag,sum)
        do i = 1,A%shape(1)
            do j = A%index_ptr(i),A%index_ptr(i+1)-1
                if (A%indices(j) == i) then
                    diag = abs(A%data(j))
                else
                    sum = sum + abs(A%data(j))
                end if
            end do
            if (sum > diag) then
                write(6,*) i,diag,sum
                not_diag_dom(i) = .true.
            else
                not_diag_dom(i) = .false.
            end if
        end do
        !$omp end parallel do

        if (any(not_diag_dom)) write(6,*) "A is not diagonally dominant!"
    end subroutine check_diag_dom
end module sparse_array_tools