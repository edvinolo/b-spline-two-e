include 'mkl_pardiso.f90'
module PARDISO_tools
    use kind_tools
    use mkl_pardiso
    implicit none

    type,public :: PARDISO_solver
        type(MKL_PARDISO_HANDLE), dimension(:), allocatable  :: pt !Pointers used by PARDISO
        integer, dimension(:), allocatable :: iparm !PARDISO parameter array
        integer :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz !PARDISO integer parameters
        logical :: full, first
        integer, dimension(1) :: idum !dummy integer array
        complex(dp), dimension(1) :: ddum !dummy complex array

        contains
        procedure :: setup => PARDISO_setup
        procedure :: factor => PARDISO_factor
        procedure :: solve => PARDISO_solve
        procedure :: cleanup => PARDISO_cleanup
    end type PARDISO_solver

    type, extends(PARDISO_solver), public :: PARDISO_solver_sp
        complex(sp), allocatable :: a(:)
        complex(sp), allocatable :: x(:)
        complex(sp), allocatable :: b(:)
        complex(sp):: sdum(1) ! dummy complex single array
        contains
        procedure :: setup => PARDISO_setup_sp
        procedure :: factor => PARDISO_factor_sp
        procedure :: solve => PARDISO_solve_sp
        procedure :: cleanup => PARDISO_cleanup_sp
    end type PARDISO_solver_sp

contains
    subroutine PARDISO_set_params(this,n,nnz,full,use_sp)
        class(PARDISO_solver), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        logical, intent(in) :: full
        logical, intent(in) :: use_sp

        integer :: i !loop index


        this%n = n
        this%nnz = nnz
        this%nrhs = 1 !Hard code this for now, could change later
        this%full = full

        this%error  = 0 !initialize error flag
#if defined(DEBUG_PARDISO)
        this%msglvl = 1 !print statistical information
#else
        this%msglvl = 0 !print no statistical information
#endif
        if (this%full) then
            this%mtype = 13 ! complex general
        else
            this%mtype  = 6 ! complex symmetric
        end if
        this%maxfct = 1 !Maximum number of factors
        this%mnum = 1 !Number of matrix to solve


        !.. Initialize the internal solver memory pointer. This is only
        ! necessary for the FIRST call of the PARDISO solver.
        allocate (this%pt(64))
        do i = 1, 64
            this%pt(i)%dummy =  0
        end do

        !..
        !.. Set up PARDISO control parameters
        !..
        allocate(this%iparm(64))
        do i = 1, 64
            this%iparm(i) = 0
        end do

        !Please see the intel MKL online documentation for PARDISO for the complete meaning of all params
        this%iparm(1) = 1 !Do not use defaults for all iparm
        this%iparm(2) = 3 !Parallel fill in reordering
        if (this%full) then
            this%iparm(10) = 13 !Perturbation of small pivots by 1e-13
        else
            this%iparm(10) = 8 !Perturbation of small pivots by 1e-8
        endif
        if (this%full) then
            this%iparm(11) = 1 !1: Enable scaling vectors, must supply numerical values during phase 11
            this%iparm(13) = 1 !1: Enable weighted matching, must supply numerical values during phase 11
        else
            this%iparm(11) = 0 !1: Enable scaling vectors
            this%iparm(13) = 0 !1: Enable weighted matching
        end if
        this%iparm(21) = 1 !Enable Bunch-Kaufman pivoting for complex symmetric matrix
        if (this%full) then
            this%iparm(24) = 0 ! Use classical factorization algorithm, must be used if 11 and 13 are used
        else
            this%iparm(24) = 10 !Change to 10 for two-level factorization (must disable param 11 and 13 if used)
        end if
#if defined(DEBUG_PARDISO)
        this%iparm(27) = 1 !Enable matrix checker, could be useful for debugging, but probably disable later
#else
        this%iparm(27) = 0 !Enable matrix checker, could be useful for debugging, but probably disable later
#endif
        if (use_sp) then
            this%iparm(28) = 1 ! Use single precision
        end if

        this%first = .true.

    end subroutine PARDISO_set_params

    subroutine PARDISO_setup(this,n,nnz,a,ia,ja,full)
        class(PARDISO_solver), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        complex(dp), dimension(nnz), intent(in) :: a(nnz)
        integer, dimension(n), intent(in) :: ia(n+1)
        integer, dimension(nnz), intent(in) :: ja(nnz)
        logical, intent(in) :: full

        logical, parameter :: use_sp = .false.

        call PARDISO_set_params(this,n,nnz,full,use_sp)

        write(6,*) ""
        write(6,*) "Performing PARDISO setup..."
        this%phase = 11 !Analysis
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO setup error: ', this%error
        end if
        write(6,*) "Done!"
    end subroutine PARDISO_setup

    subroutine PARDISO_setup_sp(this,n,nnz,a,ia,ja,full)
        class(PARDISO_solver_sp), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        complex(dp), dimension(nnz), intent(in) :: a(nnz)
        integer, dimension(n), intent(in) :: ia(n+1)
        integer, dimension(nnz), intent(in) :: ja(nnz)
        logical, intent(in) :: full

        logical, parameter :: use_sp = .true.

        call PARDISO_set_params(this,n,nnz,full,use_sp)

        allocate(this%a(nnz), this%x(n),this%b(n))

        this%a = real(a,kind=sp)

        write(6,*) ""
        write(6,*) "Performing PARDISO setup..."
        this%phase = 11 !Analysis
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%sdum, this%sdum, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO setup error: ', this%error
        end if
        write(6,*) "Done!"
    end subroutine PARDISO_setup_sp

    subroutine PARDISO_factor(this,a,ia,ja)
        class(PARDISO_solver), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia(this%n+1)
        integer, intent(in) :: ja(this%nnz)

        write(6,*) ""
        write(6,*) "Performing PARDISO factorization..."
        if (this%full.and.(.not.this%first)) then
            this%phase = 12 !Analysis and factorization
        else
            this%phase = 22 !Numerical factorization
            this%first = .false.
        end if
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)
        write(6,*) "Done!"
    end subroutine PARDISO_factor

    subroutine PARDISO_factor_sp(this,a,ia,ja)
        class(PARDISO_solver_sp), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia(this%n+1)
        integer, intent(in) :: ja(this%nnz)

        this%a = real(a,kind=sp)

        write(6,*) ""
        write(6,*) "Performing PARDISO factorization..."
        if (this%full.and.(.not.this%first)) then
            this%phase = 12 !Analysis and factorization
        else
            this%phase = 22 !Numerical factorization
            this%first = .false.
        end if
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%sdum, this%sdum, this%error)
        write(6,*) "Done!"
    end subroutine PARDISO_factor_sp

    subroutine PARDISO_solve(this,a,ia,ja,x,b)
        class(PARDISO_solver), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia (this%n+1)
        integer, intent(in) :: ja(this%nnz)
        complex(dp), intent(inout) :: x(this%n)
        complex(dp), intent(inout) :: b(this%n)

        this%phase = 33 !Compute solution
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, b, x, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO solve error: ', this%error
        end if

    end subroutine PARDISO_solve

    subroutine PARDISO_solve_sp(this,a,ia,ja,x,b)
        class(PARDISO_solver_sp), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia (this%n+1)
        integer, intent(in) :: ja(this%nnz)
        complex(dp), intent(inout) :: x(this%n)
        complex(dp), intent(inout) :: b(this%n)

        this%a = real(a,kind=sp)
        ! this%x = real(x,kind=sp)
        this%b = real(b,kind=sp)

        this%phase = 33 !Compute solution
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%b, this%x, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO solve error: ', this%error
        end if

        x = real(this%x,kind=dp)

    end subroutine PARDISO_solve_sp

    subroutine PARDISO_cleanup(this)
        class(PARDISO_solver), intent(inout) :: this

        this%phase = -1 !Release internal memory for solver
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%ddum, this%idum, this%idum, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)

        deallocate(this%pt,this%iparm)

    end subroutine PARDISO_cleanup

    subroutine PARDISO_cleanup_sp(this)
        class(PARDISO_solver_sp), intent(inout) :: this

        this%phase = -1 !Release internal memory for solver
        call pardiso_64 (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%sdum, this%idum, this%idum, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%sdum, this%sdum, this%error)

        deallocate(this%pt,this%iparm,this%a,this%x,this%b)

    end subroutine PARDISO_cleanup_sp
end module PARDISO_tools