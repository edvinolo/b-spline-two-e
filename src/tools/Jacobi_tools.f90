module Jacobi_tools
    use kind_tools
    use sparse_array_tools, only: CSR_matrix
    implicit none

    type, public :: Jacobi
        integer :: n
        complex(dp), allocatable :: inv_diag(:)
    contains
        procedure :: setup => setup_Jacobi
        procedure :: update => update_Jacobi
        procedure :: factor => factor_Jacobi
        procedure :: solve_vec => solve_Jacobi
        procedure :: solve_mat => solve_Jacobi_mat
        generic :: solve => solve_vec,solve_mat
        procedure :: cleanup => cleanup_Jacobi
    end type Jacobi
contains

    subroutine setup_Jacobi(this,A)
        class(Jacobi), intent(inout) :: this
        type(CSR_matrix), intent(in) :: A

        this%n = A%shape(1)

        allocate(this%inv_diag(this%n))

        call this%factor(A)
    end subroutine setup_Jacobi

    subroutine update_Jacobi(this,A)
        class(Jacobi), intent(inout) :: this
        type(CSR_matrix), intent(in) :: A

        call this%factor(A)
    end subroutine update_Jacobi

    subroutine cleanup_Jacobi(this)
        class(Jacobi), intent(inout) :: this

        deallocate(this%inv_diag)
    end subroutine cleanup_Jacobi

    !Computes Jacobi preconditioner of matrix in csr format
    subroutine factor_Jacobi(this,A)
       class(Jacobi), intent(inout) :: this
       type(CSR_matrix), intent(in) :: A

        integer :: k,j

        do k = 1,this%n
            do j = A%index_ptr(k),a%index_ptr(k+1)-1
                if (A%indices(j) == k) then
                    this%inv_diag = 1.0_dp/A%data(j)
                    exit
                end if
            end do
        end do
    end subroutine factor_Jacobi

    !Subroutine to apply Jacobi preconditioner
    subroutine solve_Jacobi(this,sol,rhs)
        class(Jacobi), intent(in) :: this
        complex(dp), intent(out) :: sol(this%n)
        complex(dp), intent(in) :: rhs(this%n)

        integer i

        !$omp parallel do
        do i = 1,this%n
            sol(i) = rhs(i)*this%inv_diag(i)
        end do
        !$omp end parallel do
    end subroutine solve_Jacobi

    !Apply Jacobi preconditioner to multiple RHS
    subroutine solve_Jacobi_mat(this,nrhs,sol,rhs)
        class(Jacobi), intent(in) :: this
        integer, intent(in) :: nrhs
        complex(dp), intent(out) :: sol(this%n,nrhs)
        complex(dp), intent(in) :: rhs(this%n,nrhs)

        integer i,k

        !$omp parallel do
        do i=1,nrhs
            sol(:,i) = this%inv_diag*rhs(:,i)
        end do
        !$omp end parallel do
    end subroutine solve_Jacobi_mat
end module Jacobi_tools