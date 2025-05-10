program test_sparse_mv
    use kind_tools
    use constants_tools
    use sparse_array_tools
    use block_tools, only: block_diag_CS
    use omp_lib, only: omp_get_wtime

    integer, parameter :: n = 10000
    type(CSR_matrix) :: A,B
    type(block_diag_CS) :: A_block,B_block
    ! complex(dp), allocatable :: A_dense(:,:),B_dense(:,:)
    complex(dp),allocatable :: x(:), y_a(:), y_b(:)
    integer :: i,nnz_a,nnz_b,ptr_a,ptr_b
    real(dp) :: t_1,t_2

    nnz_a = n + 2*(n-1)
    nnz_b = n + n - 1

    call A%init([n,n],nnz_a)
    call B%init([n,n],nnz_b)

    A%data(1) = -2*i_
    A%indices(1) = 1
    A%data(2) = 1*i_
    A%indices(2) = 2
    A%index_ptr(1) = 1
    ptr_a = 3

    B%data(1) = -2*i_
    B%indices(1) = 1
    B%data(2) = 1*i_
    B%indices(2) = 2
    B%index_ptr(1) = 1
    ptr_b = 3

    do i=2,n-1
        A%index_ptr(i) = ptr_a
        A%data(ptr_a) = 1*i_
        A%indices(ptr_a) = i-1
        A%data(ptr_a+1) = -2*i_
        A%indices(ptr_a+1) = i
        A%data(ptr_a+2) = 1*i_
        A%indices(ptr_a+2) = i+1
        ptr_a = ptr_a + 3

        B%index_ptr(i) = ptr_b
        B%data(ptr_b) = -2*i_
        B%indices(ptr_b) = i
        B%data(ptr_b+1) = 1*i_
        B%indices(ptr_b+1) = i+1
        ptr_b = ptr_b + 2
    end do

    A%index_ptr(n) = ptr_a
    A%data(ptr_a) = 1*i_
    A%indices(ptr_a) = n-1
    A%data(ptr_a+1) = -2*i_
    A%indices(ptr_a+1) = n
    A%index_ptr(n+1) = ptr_a+2

    B%index_ptr(n) = ptr_b
    B%data(ptr_b) = -2*i_
    B%indices(ptr_b) = n
    B%index_ptr(n+1) = ptr_b+1

    write(6,*) A%index_ptr(n+1), A%nnz
    write(6,*) B%index_ptr(n+1), B%nnz


    allocate(x(n),y_a(n),y_b(n))
    x = 1.0_dp

    t_1 = omp_get_wtime()
    call CSR_mv(A,x,y_a)
    t_2 = omp_get_wtime()
    write(6,*) t_2-t_1

    t_1 = omp_get_wtime()
    call CSR_mv_sym(B,x,y_b)
    t_2 = omp_get_wtime()
    write(6,*) t_2-t_1

    x = y_a-y_b
    write(6,*) sqrt(dot_product(x,x))/sqrt(dot_product(y_a,y_a))

    call A%deall()
    call B%deall()


    call A_block%load("S_F.dat")
    call B_block%load("S_U.dat")

    call A_block%to_CS(A,.true.)
    call B_block%to_CS(B,.true.)

    deallocate(x,y_a,y_b)
    allocate(x(A%shape(1)),y_a(A%shape(1)),y_b(A%shape(1)))
    x = 1.0_dp


    t_1 = omp_get_wtime()
    do i=1,100
        call CSR_mv(A,x,y_a)
        !call mkl_zcsrgemv('N',A%shape(1),A%data,A%index_ptr,A%indices,x,y_a)
        ! call mkl_zcsrsymv('U',A%shape(1),B%data,B%index_ptr,B%indices,x,y_a)
    end do
    t_2 = omp_get_wtime()
    write(6,*) (t_2-t_1)/100_dp

    t_1 = omp_get_wtime()
    do i = 1,100
        call CSR_mv_sym(B,x,y_b)
        ! call mkl_zcsrsymv('U',A%shape(1),B%data,B%index_ptr,B%indices,x,y_b)
    end do
    t_2 = omp_get_wtime()
    write(6,*) (t_2-t_1)/100_dp

    x = y_a-y_b
    write(6,*) sqrt(dot_product(x,x))/sqrt(dot_product(y_a,y_a))

    ! call A%get_dense(A_dense)
    ! call B%get_dense(B_dense)

    ! do i = 1,n
    !     write(6,*) A_dense(i,:)
    ! end do

    ! write(6,*)
    ! do i = 1,n
    !     write(6,*) B_dense(i,:)
    ! end do

    ! write(6,*) y_a
    ! write(6,*) y_b

end program test_sparse_mv