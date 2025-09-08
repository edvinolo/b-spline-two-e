module essential_states
    use kind_tools
    use constants_tools
    use sparse_array_tools, only: CSR_matrix, CSR_mv, CSR_mv_sym
    use block_tools, only: block_diag_CS
    use orbital_tools, only: basis
    use diagonalization, only: do_diag_essential
    use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
    implicit none

    ! complex(dp), external :: zdotu
    interface
        complex(dp) function zdotu(N,zx,incx,zy,incy)
            import :: dp
            integer :: N
            complex(dp) :: zx(*)
            integer :: incx
            complex(dp) :: zy(*)
            integer :: incy
        end function zdotu
    end interface

    type, public :: essential_state
        integer :: block
        complex(dp), allocatable :: vector(:)
        complex(dp) :: energy
        integer :: ptr(2)
    contains
        procedure :: set_pointers => set_regular_pointers
        procedure :: projection
    end type essential_state

    type, extends(essential_state), public :: Floquet_essential_state
        integer :: Floquet_block
    contains
        procedure :: set_pointers => set_Floquet_pointers
    end type Floquet_essential_state

contains
    subroutine set_regular_pointers(this,bas)
        class(essential_state), intent(inout) :: this
        type(basis), intent(in) :: bas

        integer :: i,p_1,p_2

        p_1 = 1

        do i = 1, this%block-1
            p_1 = p_1 + bas%syms(i)%n_config
        end do

        p_2 = p_1 + bas%syms(this%block)%n_config - 1

        this%ptr = [p_1,p_2]
    end subroutine set_regular_pointers

    subroutine set_Floquet_pointers(this,bas)
        class(Floquet_essential_state), intent(inout) :: this
        type(basis), intent(in) :: bas

        integer :: i,p_1,p_2,n

        n = bas%n_states
        p_1 = (this%Floquet_block-1)*n + 1

        do i = 1, this%block-1
            p_1 = p_1 + bas%syms(i)%n_config
        end do

        p_2 = p_1 + bas%syms(this%block)%n_config - 1

        this%ptr = [p_1,p_2]
    end subroutine set_Floquet_pointers

    function projection(this,vec,S,full) result(res)
        class(essential_state), intent(in) :: this
        complex(dp), intent(in) :: vec(:)
        type(CSR_matrix), intent(in) :: S
        logical, intent(in) :: full
        complex(dp) :: res

        integer :: N
        complex(dp), allocatable :: temp(:)

        N = size(this%vector)
        allocate(temp(N))

        if (full) then
            call CSR_mv(S,this%vector,temp)
        else
            call CSR_mv_sym(S,this%vector,temp)
        end if

        res = zdotu(N,vec(this%ptr(1):this%ptr(2)),1,temp,1)
    end function projection

    subroutine setup_essential_states(states,n_ess,target_blocks,bas)
        class(essential_state), allocatable, intent(inout) :: states(:)
        integer, intent(in) :: n_ess
        integer, intent(in) :: target_blocks(n_ess)
        type(basis) :: bas

        integer :: i

        allocate(essential_state :: states(n_ess))

        do i = 1,n_ess
            states(i)%block = target_blocks(i)
            call states(i)%set_pointers(bas)
        end do
    end subroutine setup_essential_states

    subroutine setup_Floquet_essential_states(states,n_ess,target_blocks,n_blocks,target_Floquet_blocks,bas)
        class(essential_state), allocatable, intent(inout) :: states(:)
        integer, intent(in) :: n_ess
        integer, intent(in) :: target_blocks(n_ess)
        integer, intent(in) :: n_blocks(2)
        integer, intent(in) :: target_Floquet_blocks(n_ess)
        type(basis) :: bas

        integer :: i

        allocate(Floquet_essential_state :: states(n_ess))

        select type (states)
        type is (Floquet_essential_state)
            do i = 1,n_ess
                states(i)%block = target_blocks(i)
                states(i)%Floquet_block = target_Floquet_blocks(i) + n_blocks(1) + 1
                call states(i)%set_pointers(bas)
            end do
        class default
            write(stderr,*) "Essential state type must be Floquet_essential_state in setup_Floquet_essential_states"
            error stop
        end select
    end subroutine setup_Floquet_essential_states

    subroutine find_essential_states(states,n_ess,targets,H_0_block,S_Block)
        class(essential_state), intent(inout) :: states(n_ess)
        integer, intent(in) :: n_ess
        complex(dp), intent(in) :: targets(n_ess)
        type(block_diag_CS), intent(in) :: H_0_block
        type(block_diag_CS), intent(in) :: S_block

        integer :: i

        do i = 1,n_ess
            call do_diag_essential(H_0_block%blocks(states(i)%block),S_block%blocks(states(i)%block),&
                                    targets(i),states(i)%vector,states(i)%energy)
        end do
    end subroutine find_essential_states

    subroutine compute_H_eff(n_ess,projections,energies,H_eff)
        integer, intent(in) :: n_ess
        complex(dp), intent(in) :: projections(n_ess,n_ess)
        complex(dp), intent(in) :: energies(n_ess)
        complex(dp), intent(out) :: H_eff(n_ess,n_ess)

        integer :: i,j
        complex(dp) :: Lambda(n_ess,n_ess)
        complex(dp) :: W(n_ess,n_ess)
        complex(dp) :: W_inv(n_ess,n_ess)
        complex(dp) :: eigs(n_ess)
        complex(dp) :: H_temp(n_ess,n_ess)

        integer :: ipiv(n_ess)
        complex(dp), allocatable :: work(:)
        integer :: lwork
        real(dp), allocatable ::  rwork(:)
        complex(dp) :: vl(1),vr(1)
        integer :: info

        W = transpose(projections)
        W_inv = W

        ! Factorize
        call zgetrf(n_ess,n_ess,W_inv,n_ess,ipiv,info)
        if (info /= 0) then
            write(stderr,*) "Error, zgetrf exited with info: ", info
            stop
        end if

        ! Workspace query for inverse computation (lwork = -1)
        allocate(work(1))
        call zgetri(n_ess,W_inv,n_ess,ipiv,work,-1,info)
        if (info /= 0) then
            write(stderr,*) "Error, zgetri exited with info: ", info
            stop
        end if
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Compute the inverse
        call zgetri(n_ess,W_inv,n_ess,ipiv,work,lwork,info)
        if (info /= 0) then
            write(stderr,*) "Error, zgetri exited with info: ", info
            stop
        end if

        Lambda = (0.0_dp,0.0_dp)
        do i = 1,n_ess
            Lambda(i,i) = energies(i)
        end do

        H_eff = matmul(W,matmul(Lambda,W_inv))

        write(stdout,*) 'H_eff: '
        do i = 1,n_ess
            do j = 1,n_ess
                write(stdout,'(a)', advance='no') '('

                if (real(H_eff(i,j),kind = dp)<0) then
                    write(stdout,'(es25.17e3)', advance='no') real(H_eff(i,j),kind = dp)
                else
                    write(stdout,'(es24.17e3)', advance='no') real(H_eff(i,j),kind = dp)
                end if

                if (aimag(H_eff(i,j))<0) then
                    write(stdout,'(a)', advance='no') ','
                    write(stdout,'(es25.17e3)', advance='no') aimag(H_eff(i,j))
                else
                    write(stdout,'(a)', advance='no') ', '
                    write(stdout,'(es24.17e3)', advance='no') aimag(H_eff(i,j))
                end if

                write(stdout,'(a)', advance='no') ') '
            end do
            write(stdout,*)''
        end do
        write(stdout,*)''

        ! Workspace query for eigenvalues
        allocate(rwork(2*n_ess))
        H_temp = H_eff
        call zgeev('N','N',n_ess,H_temp,n_ess,eigs,vl,1,vr,1,work,-1,rwork,info)
        if (info /= 0) then
            write(stderr,*) "Error, zgeev exited with info: ", info
            stop
        end if
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Compute the eigenvalues
        call zgeev('N','N',n_ess,H_temp,n_ess,eigs,vl,1,vr,1,work,lwork,rwork,info)

        write(stdout,*) 'H_eff eigenvalues:'
        do i=1,n_ess
            write(stdout,*) eigs(i)
        end do

    end subroutine compute_H_eff
end module essential_states