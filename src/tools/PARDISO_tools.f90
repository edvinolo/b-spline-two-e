include 'mkl_pardiso.f90'
module PARDISO_tools
    use kind_tools
    use mkl_pardiso
    implicit none

    type,public :: PARDISO_solver
        type(MKL_PARDISO_HANDLE), dimension(:), allocatable  :: pt !Pointers used by PARDISO
        integer, dimension(:), allocatable :: iparm !PARDISO parameter array
        integer :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz !PARDISO integer parameters
        integer, dimension(1) :: idum !dummy integer array
        complex(dp), dimension(1) :: ddum !dummy complex array

        contains
        procedure :: setup => PARDISO_setup
        procedure :: factor => PARDISO_factor
        procedure :: solve => PARDISO_solve
        procedure :: cleanup => PARDISO_cleanup
    end type PARDISO_solver

contains
    subroutine PARDISO_setup(this,n,nnz,a,ia,ja)
        class(PARDISO_solver), intent(inout) :: this
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        complex(dp), dimension(nnz), intent(in) :: A(nnz)
        integer, dimension(n), intent(in) :: ia(n+1)
        integer, dimension(nnz), intent(in) :: ja(nnz)

        integer :: i !loop index


        this%n = n
        this%nnz = nnz
        this%nrhs = 1 !Hard code this for now, could change later

        this%error  = 0 !initialize error flag
        this%msglvl = 0 !print no statistical information
        this%mtype  = 6 !complex symmetric
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
        this%iparm(10) = 8 !Perturbation of small pivots by 1e-8
        this%iparm(11) = 0 !1: Enable scaling vectors
        this%iparm(13) = 0 !1: Enable weighted matching
        this%iparm(21) = 1 !Enable Bunch-Kaufman pivoting for complex symmetric matrix
        this%iparm(24) = 10 !Change to 10 for two-level factorization (must disable param 11 and 13 if used)
        this%iparm(27) = 0 !Enable matrix checker, could be useful for debugging, but probably disable later

        ! this%phase = 12 !Analysis + numerical factorization
        ! call pardiso (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
            !   this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)

        write(6,*) ""
        write(6,*) "Performing PARDISO setup..."
        this%phase = 11 !Analysis
        call pardiso (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO setup error: ', this%error
        end if
        write(6,*) "Done!"
    end subroutine PARDISO_setup

    subroutine PARDISO_factor(this,a,ia,ja)
        class(PARDISO_solver), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia(this%n+1)
        integer, intent(in) :: ja(this%nnz)

        write(6,*) ""
        write(6,*) "Performing PARDISO factorization..."
        this%phase = 22 !Numerical factorization
        call pardiso (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)
        write(6,*) "Done!"
    end subroutine PARDISO_factor

    subroutine PARDISO_solve(this,a,ia,ja,x,b)
        class(PARDISO_solver), intent(inout) :: this
        complex(dp), intent(in) :: a(this%nnz)
        integer, intent(in) :: ia (this%n+1)
        integer, intent(in) :: ja(this%nnz)
        complex(dp), intent(inout) :: x(this%n)
        complex(dp), intent(inout) :: b(this%n)

        this%phase = 33 !Compute solution
        call pardiso (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, a, ia, ja, &
              this%idum, this%nrhs, this%iparm, this%msglvl, b, x, this%error)

        if (this%error.ne.0) then
            write(6,*) 'PARDISO solve error: ', this%error
        end if

    end subroutine PARDISO_solve

    subroutine PARDISO_cleanup(this)
        class(PARDISO_solver), intent(inout) :: this

        this%phase = -1 !Release internal memory for solver
        call pardiso (this%pt, this%maxfct, this%mnum, this%mtype, this%phase, this%n, this%ddum, this%idum, this%idum, &
              this%idum, this%nrhs, this%iparm, this%msglvl, this%ddum, this%ddum, this%error)

        deallocate(this%pt,this%iparm)

    end subroutine PARDISO_cleanup
end module PARDISO_tools