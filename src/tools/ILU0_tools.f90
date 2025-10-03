module ILU0_tools
    use kind_tools
    use sparse_array_tools, only: CSR_matrix
    implicit none

    type, public :: ILU0
        type(CSR_matrix), pointer :: A
        integer :: n
        integer :: nnz
        integer, allocatable :: diag_ptr(:)
        complex(dp), allocatable :: lu_val(:)
        complex(dp), allocatable :: temp(:)
    contains
        procedure :: setup => setup_ILU0
        procedure :: update => update_ILU0
        procedure :: factor => factor_ILU0
        procedure :: solve_vec => solve_ILU0
        procedure :: solve_mat => solve_ILU0_mat
        generic :: solve => solve_vec,solve_mat
        procedure :: cleanup => cleanup_ILU0
    end type ILU0
contains

    subroutine setup_ILU0(this,A)
        class(ILU0), intent(inout) :: this
        type(CSR_matrix), target, intent(inout) :: A

        this%A => A
        this%n = A%shape(1)
        this%nnz = A%nnz

        allocate(this%diag_ptr(this%n),this%lu_val(this%nnz),this%temp(this%n))

        call this%factor()
    end subroutine setup_ILU0

    subroutine update_ILU0(this,A)
        class(ILU0), intent(inout) :: this
        type(CSR_matrix), target, intent(inout) :: A

        this%A => A
        call this%factor()
    end subroutine update_ILU0

    subroutine cleanup_ILU0(this)
        class(ILU0), intent(inout) :: this

        call this%A%deall()
        this%A => null()
        deallocate(this%diag_ptr,this%lu_val,this%temp)
    end subroutine cleanup_ILU0

    !Computes ILU_0 of matrix in csr format
    !Adapted from fortran code of Y. Saad in section 10.3 of "Iterative methods for Sparse Linear systems"
    subroutine factor_ILU0(this)
       class(ILU0), intent(inout) :: this

        integer :: k,j,jj,j1,j2,jrow,jwork
        integer, allocatable :: work(:)
        complex(dp) :: tl

        allocate(work(this%n), source = -1)

        associate(a => this%A)
        !Main loop over rows
        !Should start from 1 if setting diag to inverse, see Saad iterative linear systems book, sec 10.4
        ! call zcopy(this%n,a%data,1,this%lu_val,1)
        this%lu_val = a%data
        ! this%diag_ptr(1) = 1
        do k = 1,this%n
            !pointers to first and last non-zero column of row k
            j1 = a%index_ptr(k)
            j2 = a%index_ptr(k+1) - 1
            !write(6,*) 'j1, j2: ', j1, j2

            do j = j1,j2
                !write(6,*) j
                work(a%indices(j)) = j
            end do

            do j = j1,j2
                !work(j) = j
                jrow = a%indices(j)!get the current column index, should be less than k
                !write(6,*) 'jrow,k : ', jrow,k

                if (jrow.ge.k) then
                    !We have reached the diagonal
                    !write(6,*) 'Hej!'
                    this%diag_ptr(k) = j

                    if ((jrow.ne.k).or.(a%data(j).eq.0)) then
                        write(6,*) 'Error, encountered zero pivot in ILU(0) at row k = ', k
                        stop
                    end if

                    !Divide the diagonals so that this is does not have to be done each time in the solve routine
                    a%data(j) = 1.d0/a%data(j)

                    !reset entires of work
                    work(a%indices(j1:j2)) = -1

                    !break j loop
                    exit
                end if

                !divide by diagonal, update a_k,jrow
                tl = a%data(j)*a%data(this%diag_ptr(jrow)) !use this version if the inverses of diagonal elements are used
                !tl = values(j)/values(diag_ptr(jrow))
                a%data(j) = tl

                !update values to the right of a_k,jrow
                !a_k,m = a_k,m - a_k,jrow*a_jrow,m
                !write(6,*) diag_ptr
                do jj = this%diag_ptr(jrow) + 1, a%index_ptr(jrow+1)-1
                    jwork = work(a%indices(jj))
                    !write(6,*) "Col index: " ,col_index(jj), jj
                    if (jwork.ge.0) then
                        a%data(jwork) = a%data(jwork) - tl*a%data(jj)
                    end if
                end do

            end do

        end do
        end associate
    end subroutine factor_ILU0

    !Subroutine to perform the triangular solves on single RHS after ILU
    !Adapted from fortran code of Y. Saad in section 10.2 of "Iterative methods for Sparse Linear systems"
    !Assumes that inverse of diag of U has already been computed
    subroutine solve_ILU0(this,sol,rhs)
        class(ILU0), intent(in) :: this
        complex(dp), intent(out) :: sol(this%n)
        complex(dp), intent(in) :: rhs(this%n)

        integer i,k

        associate(a => this%A)
        !Forward solve L*sol = rhs
        !Think if it would be faster to store rhs/sol by row? When testing it looks like it, so change how it is recieved
        do i = 1,this%n
            sol(i) = rhs(i)
            do k = a%index_ptr(i), this%diag_ptr(i)-1
                !sol(i,:) = sol(i,:)-lu_values(k)*sol(col_index(k),:)
                sol(i) = sol(i)-a%data(k)*sol(a%indices(k))
            end do
        end do

        !Backward solve sol = inv(U)*sol
        do i = this%n,1,-1
            do k = this%diag_ptr(i) + 1, a%index_ptr(i+1) - 1
                !sol(i,:) = sol(i,:) - lu_values(k)*sol(col_index(k),:)
                sol(i) = sol(i) - a%data(k)*sol(a%indices(k))
            end do

            !Here is where it is assumed that the inverse of the diagonal has been computed
            !sol(i,:) = lu_values(diag_ptr(i))*sol(i,:)
            sol(i) = a%data(this%diag_ptr(i))*sol(i)
        end do

        ! call mkl_zcsrtrsv('L','N','U',this%n,this%lu_val,a%index_ptr,a%indices,rhs,this%temp)
        ! call mkl_zcsrtrsv('U','N','N',this%n,this%lu_val,a%index_ptr,a%indices,this%temp,sol)
        end associate

    end subroutine solve_ILU0

    !Subroutine to perform the triangular solves on n_RHS after ILU
    !Adapted from fortran code of Y. Saad in section 10.2 of "Iterative methods for Sparse Linear systems"
    !Assumes that inverse of diag of U has already been computed
    subroutine solve_ILU0_mat(this,nrhs,sol,rhs)
        class(ILU0), intent(in) :: this
        integer, intent(in) :: nrhs
        complex(dp), intent(out) :: sol(this%n,nrhs)
        complex(dp), intent(in) :: rhs(this%n,nrhs)

        integer i !,k

        associate(a => this%A)
        ! !Forward solve L*sol = rhs
        ! !Think if it would be faster to store rhs/sol by row? When testing it looks like it, so change how it is recieved
        ! do i = 1,this%n
        !     sol(i,:) = rhs(i,:)
        !     do k = a%index_ptr(i), this%diag_ptr(i)-1
        !         !sol(i,:) = sol(i,:)-lu_values(k)*sol(col_index(k),:)
        !         sol(:,i) = sol(:,i)-a%data(k)*sol(:,a%indices(k))
        !     end do
        ! end do

        ! !Backward solve sol = inv(U)*sol
        ! do i = this%n,1,-1
        !     do k = this%diag_ptr(i) + 1, a%index_ptr(i+1) - 1
        !         !sol(i,:) = sol(i,:) - lu_values(k)*sol(col_index(k),:)
        !         sol(:,i) = sol(:,i) - a%data(k)*sol(:,a%indices(k))
        !     end do

        !     !Here is where it is assumed that the inverse of the diagonal has been computed
        !     !sol(i,:) = lu_values(diag_ptr(i))*sol(i,:)
        !     sol(:,i) = a%data(this%diag_ptr(i))*sol(:,i)
        ! end do
        do i  = 1,nrhs
            call mkl_zcsrtrsv('L','N','U',this%n,this%lu_val,a%index_ptr,a%indices,rhs(:,i),this%temp)
            call mkl_zcsrtrsv('U','N','N',this%n,this%lu_val,a%index_ptr,a%indices,this%temp,sol(:,i))
        end do
        end associate
    end subroutine solve_ILU0_mat
end module ILU0_tools