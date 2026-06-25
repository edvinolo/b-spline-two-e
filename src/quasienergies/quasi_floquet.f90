module quasi_floquet
    use kind_tools
    use constants_tools
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_CS,block_diag_CS,AXPBY,APX,XPAY
    use orbital_tools
    use input_tools, only: gauge, n_blocks
    implicit none

contains

    subroutine setup_dip_floquet(dip,D)
        type(block_CS), intent(inout) :: dip
        type(block_CS), intent(out) :: D

        integer :: n,i,j,ptr,ptr_i,n_block

        n = dip%block_shape(1)*(n_blocks(1)+n_blocks(2) + 1)
        n_block = dip%block_shape(1)

        call D%init([n,n])

        ! Setup the shapes of the blocks
        do i = 1,n
            ptr_i = mod(i,n_block)
            if (ptr_i == 0) ptr_i = n_block
            do j = i,n
                ptr = mod(j,n_block)
                if (ptr == 0) ptr = n_block
                D%blocks(i,j)%shape = dip%blocks(ptr_i,ptr)%shape
                if (i /= j) then
                    D%blocks(j,i)%shape(1) =  D%blocks(i,j)%shape(2)
                    D%blocks(j,i)%shape(2) =  D%blocks(i,j)%shape(1)
                end if
            end do
        end do

        ptr = 1
        ptr_i = 1
        do i = -n_blocks(1),n_blocks(2)
            ! Add the parts in the A blocks (see MCR PRA 27 6 1983 for ref.)
            do j = 1,n_block-1
                D%blocks(ptr,ptr+1) = dip%blocks(j,j+1)
                ptr = ptr + 1
            end do
            ptr = ptr + 1

            ! Add the parts in the B blocks (I use B^T instead of B in MCR, beacause I order the blocks differently)
            if (i /= n_blocks(2)) then
                do j = 1,n_block-1
                    if (mod(j,2) == 0) then
                        call dip%blocks(j,j+1)%transp(D%blocks((ptr_i-1)*n_block + j+1,ptr_i*n_block+j))
                    else
                        D%blocks((ptr_i-1)*n_block + j,ptr_i*n_block+j+1) = dip%blocks(j,j+1)
                    end if
                end do
                ptr_i = ptr_i + 1
            end if
        end do

        ! do i = 1,n
        !     do j = 1,n
        !         write(6,'(l1)',advance='no') D%blocks(i,j)%arrays_allocated()
        !     end do
        !     write(6,*)''
        ! end do

        if (gauge == 'v') call D%scale(-i_)

    end subroutine setup_dip_floquet

    subroutine setup_H_floquet(H_block,bas, H_0, S_block, D, omega, shift, V_0, B_subset)
        type(block_CS), intent(inout) :: H_block
        type(basis), intent(in) :: bas
        type(block_diag_CS), intent(in) ::  H_0
        type(block_diag_CS), intent(in) :: S_block
        type(block_CS), intent(in) :: D
        real(dp), intent(in) :: omega
        complex(dp), intent(in) :: shift
        real(dp), intent(in) :: V_0
        logical, intent(in) :: B_subset

        complex(dp), allocatable :: shifts(:)
        integer :: n,i,j,ptr
        type(block_diag_CS) :: A_0
        type(block_diag_CS) :: H_0_block

        n = H_0%block_shape(1)*(n_blocks(1) + n_blocks(2) + 1)

        A_0 = H_0
        allocate(shifts(H_0%block_shape(1)),source = -shift)
        do i = 1, H_0%block_shape(1)
            if (bas%syms(i)%pi .eqv. .true.) then
                shifts(i) = shifts(i) - omega
            end if
        end do

        call APX(A_0,shifts,S_block,B_subset)

        deallocate(shifts)
        allocate(shifts(n),source = (0.0_dp,0.0_dp))
        call H_0_block%init(n)
        ptr = 1
        do i = -n_blocks(1), n_blocks(2)
            do j = 1,H_0%block_shape(1)
                H_0_block%blocks(ptr) = A_0%blocks(j)
                shifts(ptr) = 2*i*omega
                call H_0_block%blocks(ptr)%shift_B(shifts(ptr),S_block%blocks(j),B_subset)
                ptr = ptr + 1
            end do
        end do

        H_block = XPAY(H_0_block, D, cmplx(V_0*0.5_dp,kind=dp))
        call H_block%compute_shape()
    end subroutine setup_H_floquet

    subroutine setup_S_floquet(S_block,S)
        type(block_diag_CS), intent(inout) :: S_block
        type(CSR_matrix), intent(out) :: S

        type(block_diag_CS) :: temp_S
        integer :: i,n,ptr
        logical :: deall

        n = S_block%block_shape(1)*(n_blocks(1)+n_blocks(2) + 1)

        deall = .false.
        call temp_S%init(n)

        do i = 1, n
            ptr = mod(i,S_block%block_shape(1))
            if (ptr == 0) ptr = S_block%block_shape(1)
            temp_S%blocks(i) = S_block%blocks(ptr)
        end do

        call temp_S%to_CS(S,deall)
    end subroutine setup_S_floquet

end module quasi_floquet