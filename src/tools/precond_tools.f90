module precond_tools
    use kind_tools
    use sparse_array_tools
    use block_tools
    use PARDISO_tools, only: PARDISO_solver
    use ILU0_tools, only: ILU0
    use Jacobi_tools, only: Jacobi
    use input_tools, only: n_precond,n_relevant,couple_pq
    implicit none

    type, abstract, public :: block_PC
    contains
        procedure(setup), deferred :: setup
        procedure(update), deferred :: update
        procedure(solve), deferred :: solve
        procedure(cleanup), deferred :: cleanup
    end type block_PC

    abstract interface
        subroutine setup(this,H)
            import :: block_PC
            import :: block_CS
            class(block_PC), intent(inout) :: this
            type(block_CS), intent(in) :: H
        end subroutine setup

        subroutine update(this,H,first_in)
            import :: block_PC
            import :: block_CS
            class(block_PC), intent(inout) :: this
            type(block_CS), intent(in) :: H
            logical, optional, intent(in) :: first_in
        end subroutine update

        subroutine solve(this,x,b)
            import :: block_PC
            import :: dp
            class(block_PC), intent(inout) :: this
            complex(dp), intent(out) :: x(:)
            complex(dp), intent(inout) :: b(:)
        end subroutine solve

        subroutine cleanup(this)
            import :: block_PC
            class(block_PC), intent(inout) :: this
        end subroutine cleanup
    end interface

    type, extends(block_PC), public :: block_PQ_PC
        integer :: n_p
        integer :: n_q
        logical :: couple_pq
        type(CSR_matrix) :: H_pp
        type(CSR_matrix) :: H_pq
        type(block_diag_CS) :: H_qq_0
        type(CSR_matrix) :: H_qq
        type(PARDISO_solver) :: solv_pp
        ! type(ILU0) :: solv_pp
        ! type(PARDISO_solver), allocatable :: solv_qq(:)
        type(PARDISO_solver) :: solv_qq
        ! type(ILU0) :: solv_qq
        ! type(Jacobi) :: solv_qq
        integer, allocatable :: block_ptr(:)
        complex(dp), allocatable :: t_1(:),t_2(:),t_3(:)
    contains
        procedure :: setup => setup_block_PQ_PC
        procedure :: update => update_block_PQ_PC
        procedure :: solve => solve_block_PQ_PC
        procedure :: cleanup => cleanup_block_PQ_PC
    end type block_PQ_PC

    type, extends(block_PC), public :: block_Jacobi_PC
        integer :: n
        type(block_diag_CS) :: H_0_block
        type(CSR_matrix) :: H_0
        type(PARDISO_solver) :: solv
    contains
        procedure :: setup => setup_block_Jacobi_PC
        procedure :: update => update_block_Jacobi_PC
        procedure :: solve => solve_block_Jacobi_PC
        procedure :: cleanup => cleanup_block_Jacobi_PC
    end type block_Jacobi_PC
contains
    subroutine setup_block_PQ_PC(this,H)
        class(block_PQ_PC), intent(inout) :: this
        type(block_CS), intent(in) :: H

        integer :: i,ptr

        ! Set scalar components
        this%n_p = n_precond
        this%n_q = n_relevant-n_precond
        this%couple_pq = couple_pq

        call this%update(H,first_in=.true.)

        ! Setup pointers to p and q vectors
        allocate(this%block_ptr(0:this%n_q+1))
        this%block_ptr(0) = 1
        ptr = 1 + this%H_pp%shape(1)
        this%block_ptr(1) = ptr
        do i = 1,this%n_q
            ptr = ptr + this%H_qq_0%blocks(i)%shape(1)
            this%block_ptr(i+1) = ptr
        end do

        ! Allocate work arrays
        allocate(this%t_1(H%shape(1)),this%t_2(H%shape(1)),this%t_3(H%shape(1)))

    end subroutine setup_block_PQ_PC

    subroutine update_block_PQ_PC(this,H,first_in)
        class(block_PQ_PC), intent(inout) :: this
        type(block_CS), intent(in) :: H
        logical, optional, intent(in) :: first_in

        logical :: first
        integer :: i,j
        type(block_CS) :: H_pp_block
        type(block_CS) :: H_pq_block

        if (present(first_in)) then
            first = first_in
        else
            first = .false.
        end if

        ! Assemble the P-space Hamiltonian
        call H_pp_block%init([this%n_p,this%n_p])
        do j = 1,this%n_p
            do i = 1,this%n_p
                H_pp_block%blocks(i,j) = H%blocks(i,j)
            end do
        end do
        call H_pp_block%to_CS(this%H_pp,.true.)

        ! Assemble the 0:th order Hamiltonian in Q-space
        if (first) call this%H_qq_0%init(this%n_q)
        do i = 1,this%n_q
            this%H_qq_0%blocks(i) = H%blocks(this%n_p+i,this%n_p+i)
        end do
        call this%H_qq_0%to_CS(this%H_qq,.true.)

        ! Assemble the coupling from Q to P
        if (this%couple_pq) then
            call H_pq_block%init([this%n_p,this%n_q])
            do j = 1,this%n_q
                do i = 1,this%n_p
                    H_pq_block%blocks(i,j) = H%blocks(i,this%n_p+j)
                end do
            end do
            call H_pq_block%to_CS(this%H_pq,.true.)
        end if

        ! Factorize H_pp
        if (first) then
            call this%solv_pp%setup(this%H_pp%shape(1),this%H_pp%nnz,this%H_pp%data,&
                                    this%H_pp%index_ptr,this%H_pp%indices)
            ! call this%solv_pp%setup(this%H_pp)
        else
            ! call this%solv_pp%update(this%H_pp)
        end if
        call this%solv_pp%factor(this%H_pp%data,this%H_pp%index_ptr,this%H_pp%indices)

        ! Factorize H_qq_0
        if (first) then
            call this%solv_qq%setup(this%H_qq%shape(1),this%H_qq%nnz,this%H_qq%data,&
                                    this%H_qq%index_ptr,this%H_qq%indices)
            ! call this%solv_qq%setup(this%H_qq)
        else
            ! call this%solv_qq%update(this%H_qq)
        end if
        call this%solv_qq%factor(this%H_qq%data,this%H_qq%index_ptr,this%H_qq%indices)
        ! if (first) allocate(this%solv_qq(this%n_q))
        ! associate(H_qq => this%H_qq_0)
        ! do i = 1,this%n_q
        !     if (first) then
        !         call this%solv_qq(i)%setup(H_qq%blocks(i)%shape(1),H_qq%blocks(i)%nnz,H_qq%blocks(i)%data,&
        !                                 H_qq%blocks(i)%index_ptr,H_qq%blocks(i)%indices)
        !         ! call this%solv_qq(i)%setup(H_qq%blocks(i))
        !     else
        !         ! call this%solv_qq(i)%update(H_qq%blocks(i))
        !     end if
        !     call this%solv_qq(i)%factor(H_qq%blocks(i)%data,H_qq%blocks(i)%index_ptr,H_qq%blocks(i)%indices)
        ! end do
        ! end associate
    end subroutine update_block_PQ_PC

    subroutine solve_block_PQ_PC(this,x,b)
        class(block_PQ_PC), intent(inout) :: this
        complex(dp), intent(out) :: x(:)
        complex(dp), intent(inout) :: b(:)

        integer :: i,i_1,i_2,j_1,j_2

        associate(H_pp => this%H_pp,&
                  H_qq_0 => this%H_qq_0,&
                  H_pq => this%H_pq,&
                  t_1 => this%t_1(:),&
                  t_2 => this%t_2(:),&
                  t_3 => this%t_3(:))

        call zcopy(size(b),b,1,t_1,1)

        if (this%couple_pq) then
            i_1 = this%block_ptr(0)
            i_2 = this%block_ptr(1)-1

            ! Compute t_2_p = H_pp^-1*t_1_p
            call this%solv_pp%solve(H_pp%data,H_pp%index_ptr,H_pp%indices,t_2(i_1:i_2),t_1(i_1:i_2))
            ! call this%solv_pp%solve(t_2(i_1:i_2),t_1(i_1:i_2))

            ! Compute t_1_q = t_1_q - H_qp*t_2_p
            j_1 = i_1
            j_2 = i_2
            i_1 = this%block_ptr(1)
            i_2 = this%block_ptr(this%n_q+1)-1
            call CSR_T_mv(H_pq,t_2(j_1:j_2),t_2(i_1:i_2))

            !$omp parallel do
            do i = i_1,i_2
                t_1(i) = t_1(i) - t_2(i)
            end do
            !$omp end parallel do

            ! get very strange results (infs and nans) when using zaxpy, not sure why
            ! call zaxpy(i_2-i_1+1,(-1_dp,0_dp),t_2(i_1:i_2),1,t_1(i_1:i_2),1)
        end if

        ! Compute t_2_p = H_pp^-1*t_1_p
        i_1 = this%block_ptr(0)
        i_2 = this%block_ptr(1)-1
        call this%solv_pp%solve(H_pp%data,H_pp%index_ptr,H_pp%indices,t_2(i_1:i_2),t_1(i_1:i_2))
        ! call this%solv_pp%solve(t_2(i_1:i_2),t_1(i_1:i_2))

        ! Compute t_2_q = H_qq_0^-1*t_1_q
        i_1 = this%block_ptr(1)
        call this%solv_qq%solve(this%H_qq%data,this%H_qq%index_ptr,this%H_qq%indices,t_2(i_1:),t_1(i_1:))
        ! call this%solv_qq%solve(t_2(i_1:),t_1(i_1:))

        if (this%couple_pq) then
            i_1 = this%block_ptr(0)
            i_2 = this%block_ptr(1)-1
            j_1 = this%block_ptr(1)
            j_2 = this%block_ptr(this%n_q+1)-1

            ! Compute t_1_p = H_pq*t_2_q
            call CSR_mv(H_pq,t_2(j_1:j_2),t_1(i_1:i_2))
            call this%solv_pp%solve(H_pp%data,H_pp%index_ptr,H_pp%indices,t_3(i_1:i_2),t_1(i_1:i_2))
            ! call this%solv_pp%solve(t_3(i_1:i_2),t_1(i_1:i_2))

            !$omp parallel do
            do i = i_1,i_2
                t_2(i) = t_2(i)-t_3(i)
            end do
            !$omp end parallel do
            ! call zaxpy(i_2-i_1+1,(-1_dp,0_dp),t_3,1,t_2(i_1:i_2),1)
        end if

        call zcopy(size(x),t_2,1,x,1)
        end associate
    end subroutine solve_block_PQ_PC

    subroutine cleanup_block_PQ_PC(this)
        class(block_PQ_PC), intent(inout) :: this

        call this%H_pp%deall()
        call this%solv_pp%cleanup()

        call this%H_qq_0%deall()
        call this%H_qq%deall()
        call this%solv_qq%cleanup()
        deallocate(this%block_ptr,this%t_1,this%t_2,this%t_3)

        if (this%couple_pq) call this%H_pq%deall()
    end subroutine cleanup_block_PQ_PC

    subroutine setup_block_Jacobi_PC(this,H)
        class(block_Jacobi_PC), intent(inout) :: this
        type(block_CS), intent(in) :: H

        ! Set scalar components
        this%n = H%block_shape(1)

        ! Do first factorization
        call this%update(H,first_in=.true.)
    end subroutine setup_block_Jacobi_PC

    subroutine update_block_Jacobi_PC(this,H,first_in)
        class(block_Jacobi_PC), intent(inout) :: this
        type(block_CS), intent(in) :: H
        logical, optional, intent(in) :: first_in

        logical :: first
        integer :: i

        if (present(first_in)) then
            first = first_in
        else
            first = .false.
        end if

        ! Assemble the 0:th order Hamiltonian
        if (first) call this%H_0_block%init(this%n)
        do i = 1,this%n
            this%H_0_block%blocks(i) = H%blocks(i,i)
        end do
        call this%H_0_block%to_CS(this%H_0,.true.)

        ! Factorize H_0
        if (first) then
            call this%solv%setup(this%H_0%shape(1),this%H_0%nnz,this%H_0%data,&
                                    this%H_0%index_ptr,this%H_0%indices)
        end if
        call this%solv%factor(this%H_0%data,this%H_0%index_ptr,this%H_0%indices)
    end subroutine update_block_Jacobi_PC

    subroutine solve_block_Jacobi_PC(this,x,b)
        class(block_Jacobi_PC), intent(inout) :: this
        complex(dp), intent(out) :: x(:)
        complex(dp), intent(inout) :: b(:)

        call this%solv%solve(this%H_0%data,this%H_0%index_ptr,this%H_0%indices,x,b)
    end subroutine solve_block_Jacobi_PC

    subroutine cleanup_block_Jacobi_PC(this)
        class(block_Jacobi_PC), intent(inout) :: this

        call this%H_0%deall()
        call this%solv%cleanup()

    end subroutine cleanup_block_Jacobi_PC

end module precond_tools