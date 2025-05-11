module GMRES_tools
    use kind_tools
    use constants_tools, only: i_
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    use sparse_array_tools, only: CSR_dsymv
    use precond_tools
    implicit none

    type, public :: zFGMRES
        integer :: n
        procedure(mv_GMRES), pointer :: mv
        logical :: full
        real(dp), allocatable :: a(:) ! real part of matrix
        integer, allocatable :: ia(:), ja(:)
        real(dp), allocatable :: b(:) ! imaginary part of matrix. Could be chosen to be smaller than a, since only a small subset of matrix elements have imaginary part (from CAP, if I am clever)
        integer, allocatable :: ib(:), jb(:)
        real(dp) :: tol
        integer :: max_iter
        real(dp), allocatable :: tmp(:)
        real(dp), allocatable :: x_work(:),b_work(:)
        complex(dp), allocatable :: zwork(:)
        integer :: ipar(128)
        real(dp) :: dpar(128)
    contains
        procedure :: setup => setup_zFGMRES
        procedure :: update => update_zFGMRES
        procedure :: solve => solve_zFGMRES
        procedure :: cleanup => cleanup_zFGMRES
    end type zFGMRES

    abstract interface
        subroutine mv_GMRES(this,x,y)
            import :: zFGMRES
            import :: dp
            class(zFGMRES), intent(in) :: this
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine mv_GMRES
    end interface

    contains

        subroutine setup_zFGMRES(this,C,full,tol,max_iter)
            class(zFGMRES), intent(inout) :: this
            type(CSR_matrix), intent(in) :: C
            logical, intent(in) :: full
            real(dp), optional, intent(in) :: tol
            integer, optional, intent(in) :: max_iter

            if (present(tol)) then
                this%tol = tol
            else
                this%tol = 1e-14_dp
            end if

            if (present(max_iter)) then
                this%max_iter = max_iter
            else
                this%max_iter = 150
            end if

            this%full = full
            if (full) then
                this%mv => gemv_zFGMRES
            else
                this%mv => symv_zFGMRES
            end if

            this%a = real(C%data,kind=dp)
            this%b = aimag(C%data)

            this%ia = C%index_ptr
            ! this%ib = C%index_ptr

            this%ja = C%indices
            ! this%jb = C%indices

            this%n = C%shape(1)

            if (.not.allocated(this%tmp)) then
                allocate(this%tmp((2*this%max_iter+1)*2*this%n + this%max_iter*(this%max_iter+9)/2+1))
            end if

            if (.not.allocated(this%x_work)) then
                allocate(this%x_work(2*this%n),this%b_work(2*this%n),this%zwork(2*this%n))
            end if
        end subroutine setup_zFGMRES

        subroutine update_zFGMRES(this,C)
            class(zFGMRES), intent(inout) :: this
            type(CSR_matrix), intent(in) :: C

            if ( this%n /= C%shape(1)) then
                write(stderr,*) 'Error! Incompatible shape for C in zFGMRES update!'
                write(stderr,*) this%n, C%shape
                stop
            end if

            this%a = real(C%data,kind=dp)
            this%b = aimag(C%data)

            ! this%ia = C%index_ptr
            ! this%ib = C%index_ptr

            ! this%ja = C%indices
            ! this%jb = C%indices
        end subroutine update_zFGMRES

        subroutine symv_zFGMRES(this,x,y)
            class(zFGMRES), intent(in) :: this
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)

            integer :: i,j
            ! real(dp), allocatable :: x_r(:),x_i(:),y_r(:),y_i(:),temp(:)

            associate(n => this%n)

            ! allocate(x_r(n),x_i(n),y_r(n),y_i(n),temp(n))

            ! x_r = x(1:n)
            ! x_i = x(n+1:)

            ! ! Real part y_r = a*x_r - b*x_i
            ! call CSR_dsymv(n,this%a,this%ia,this%ja,1.0_dp,x(1:n),0.0_dp,y(1:n))
            ! call CSR_dsymv(n,this%b,this%ia,this%ja,-1.0_dp,x(n+1:),1.0_dp,y(1:n))

            ! ! Imaginary part y_i = b*x_r + a*x_i
            ! call CSR_dsymv(n,this%b,this%ia,this%ja,1.0_dp,x(1:n),0.0_dp,y(n+1:))
            ! call CSR_dsymv(n,this%a,this%ia,this%ja,1.0_dp,x(n+1:),1.0_dp,y(n+1:))

            ! call mkl_dcsrsymv('U',n,this%a,this%ia,this%ja,x_r,y_r)
            ! call mkl_dcsrsymv('U',n,this%b,this%ia,this%ja,x_i,temp)
            ! y_r = y_r - temp

            ! call mkl_dcsrsymv('U',n,this%a,this%ia,this%ja,x_i,y_i)
            ! call mkl_dcsrsymv('U',n,this%b,this%ia,this%ja,x_r,temp)
            ! y_i = y_i + temp

            ! y(1:n) = y_r
            ! y(n+1:) = y_i

            !$omp parallel private(j)
            !$omp do
            do i = 1,n
                y(i) = 0
                y(n+i) = 0
                do j = this%ia(i),this%ia(i+1)-1
                    y(i) = y(i) + this%a(j)*x(this%ja(j))
                    y(i) = y(i) - this%b(j)*x(n+this%ja(j))
                    y(n+i) = y(n+i) + this%a(j)*x(n+this%ja(j))
                    y(n+i) = y(n+i) + this%b(j)*x(this%ja(j))
                end do
            end do
            !$omp end do

            !$omp barrier

            !$omp do reduction(+:y)
            do i = 1,n
                do j = this%ia(i)+1,this%ia(i+1)-1
                    y(this%ja(j)) = y(this%ja(j)) + this%a(j)*x(i)
                    y(this%ja(j)) = y(this%ja(j)) - this%b(j)*x(n+i)
                    y(n+this%ja(j)) = y(n+this%ja(j)) + this%a(j)*x(n+i)
                    y(n+this%ja(j)) = y(n+this%ja(j)) + this%b(j)*x(i)
                end do
            end do
            !$omp end do
            !$omp end parallel

            end associate
        end subroutine symv_zFGMRES

        subroutine gemv_zFGMRES(this,x,y)
            class(zFGMRES), intent(in) :: this
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: y(:)

            integer :: i,j
            ! real(dp), allocatable :: x_r(:),x_i(:),y_r(:),y_i(:),temp(:)

            associate(n => this%n)

            ! allocate(x_r(n),x_i(n),y_r(n),y_i(n),temp(n))

            ! x_r = x(1:n)
            ! x_i = x(n+1:)

            ! ! Real part y_r = a*x_r - b*x_i
            ! call CSR_dsymv(n,this%a,this%ia,this%ja,1.0_dp,x(1:n),0.0_dp,y(1:n))
            ! call CSR_dsymv(n,this%b,this%ia,this%ja,-1.0_dp,x(n+1:),1.0_dp,y(1:n))

            ! ! Imaginary part y_i = b*x_r + a*x_i
            ! call CSR_dsymv(n,this%b,this%ia,this%ja,1.0_dp,x(1:n),0.0_dp,y(n+1:))
            ! call CSR_dsymv(n,this%a,this%ia,this%ja,1.0_dp,x(n+1:),1.0_dp,y(n+1:))

            ! call mkl_dcsrsymv('U',n,this%a,this%ia,this%ja,x_r,y_r)
            ! call mkl_dcsrsymv('U',n,this%b,this%ia,this%ja,x_i,temp)
            ! y_r = y_r - temp

            ! call mkl_dcsrsymv('U',n,this%a,this%ia,this%ja,x_i,y_i)
            ! call mkl_dcsrsymv('U',n,this%b,this%ia,this%ja,x_r,temp)
            ! y_i = y_i + temp

            ! y(1:n) = y_r
            ! y(n+1:) = y_i

            !$omp parallel do private(j)
            do i = 1,n
                y(i) = 0
                y(n+i) = 0
                do j = this%ia(i),this%ia(i+1)-1
                    y(i) = y(i) + this%a(j)*x(this%ja(j))
                    y(i) = y(i) - this%b(j)*x(n+this%ja(j))
                    y(n+i) = y(n+i) + this%a(j)*x(n+this%ja(j))
                    y(n+i) = y(n+i) + this%b(j)*x(this%ja(j))
                end do
            end do
            !$omp end parallel do

            end associate
        end subroutine gemv_zFGMRES

        subroutine solve_zFGMRES(this,x,b,precond)
            class(zFGMRES), intent(inout) :: this
            complex(dp), intent(out) :: x(:)
            complex(dp), intent(in) :: b(:)
            type(block_PC), intent(inout) :: precond

            integer :: ido,x_1,x_2,y_1,y_2,i,itercount
            real(dp), allocatable :: res(:)

            call precond%solve(x,b)
            !$omp parallel do
            do i = 1,this%n
                this%b_work(i) = real(b(i),kind=dp)
                this%b_work(this%n+i) = aimag(b(i))
                this%x_work(i) = real(x(i),kind=dp)
                this%x_work(this%n+i) = aimag(x(i))
            end do
            !$omp end parallel do

            ! this%x_work = this%b_work
            call dfgmres_init(2*this%n,this%x_work,this%b_work,ido,this%ipar,this%dpar,this%tmp)
            if (ido /= 0) then
                write(stderr,*) "Error! dfgmres_init exited with RCI_request: ", ido
                write(stderr,*) "Please consult dfgmres documentation"
                stop
            end if

            this%ipar(5) = this%max_iter
            this%ipar(15) = this%max_iter
            this%ipar(11) = 1 !Use preconditioner
            ! this%ipar(9) = 1 ! dfgmres performs residual test
            ! this%ipar(10) = 0 ! dfgmres does not request user to perform residual tests

            this%dpar(1) = this%tol ! set relative tolerance

            call dfgmres_check(2*this%n,this%x_work,this%b_work,ido,this%ipar,this%dpar,this%tmp)
            if (ido /= 0) then
                write(stderr,*) "Warning! dfgmres_check exited with RCI_request: ", ido
                write(stderr,*) "Please consult dfgmres documentation"
            end if

            ido = -1
            do while(ido /= 0)
                call dfgmres(2*this%n,this%x_work,this%b_work,ido,this%ipar,this%dpar,this%tmp)
                ! write(stdout,*) ido

                if (ido == 1) then
                    ! Apply matvec
                    x_1 = this%ipar(22)
                    x_2 = this%ipar(22) + 2*this%n - 1
                    y_1 = this%ipar(23)
                    y_2 = this%ipar(23) + 2*this%n - 1
                    call this%mv(this%tmp(x_1:x_2),this%tmp(y_1:y_2))

                else if (ido == 2) then
                    ! Check if residual less than tolerance
                    ! write(stdout,*) this%dpar(5),this%dpar(4)
                    if (this%dpar(5) <= this%dpar(4)) then
                        ! ido = 0
                        exit
                    end if

                else if (ido == 3) then
                    ! Apply precond
                    x_1 = this%ipar(22)
                    x_2 = this%ipar(22) + this%n

                    !$omp parallel do
                    do i = 1,this%n
                        this%zwork(i) = cmplx(this%tmp(x_1 + i - 1),this%tmp(x_2 + i - 1),kind=dp)
                    end do
                    !$omp end parallel do

                    call precond%solve(this%zwork(this%n+1:2*this%n),this%zwork(1:this%n))

                    y_1 = this%ipar(23)
                    y_2 = this%ipar(23) + this%n

                    !$omp parallel do
                    do i = 1,this%n
                        this%tmp(y_1 + i - 1) = real(this%zwork(this%n + i), kind=dp)
                        this%tmp(y_2 + i - 1) = aimag(this%zwork(this%n + i))
                    end do
                    !$omp end parallel do

                else if (ido == 4) then
                    ! Check norm of new vector
                    if (this%dpar(7) <= this%dpar(8)) then
                        write(stdout,*) this%dpar(7), this%dpar(8)
                        ! ido = 0
                        exit
                    end if

                else
                    exit
                end if
            end do

            if (ido < 0) then
                write(stderr,*) "Error! dfgmres exited with RCI_request: ", ido
                write(stderr,*) "Please consult dfgmres documentation"
                stop
            end if

            call dfgmres_get(2*this%n,this%x_work,this%b_work,ido,this%ipar,this%dpar,this%tmp,itercount)

            allocate(res(2*this%n))
            call this%mv(this%x_work,res)
            res = res-this%b_work

            if (ido /= 0) then
                write(stderr,*) "Error! dfgmres_get exited with RCI_request: ", ido
                write(stderr,*) "Please consult dfgmres documentation"
                stop
            end if

            !$omp parallel do
            do i = 1,this%n
                x(i) = cmplx(this%x_work(i),this%x_work(i + this%n),kind=dp)
            end do
            !$omp end parallel do
        end subroutine solve_zFGMRES

        subroutine cleanup_zFGMRES(this)
            class(zFGMRES), intent(inout) :: this

            deallocate(this%a,this%ia,this%ja,this%b)!,this%ib,this%jb)
            deallocate(this%b_work,this%x_work,this%zwork,this%tmp)
        end subroutine cleanup_zFGMRES

end module GMRES_tools