module diagonalization
    use kind_tools
    use constants_tools
    use input_tools
    use dir_tools
    use orbital_tools, only: basis
    use sparse_array_tools, only: CSR_matrix
    use block_tools, only: block_diag_CS
    use PARDISO_tools, only: PARDISO_solver
    use eig_tools, only: drive_ARPACK_SI,sort_eig
    implicit none

    type, public :: eig_results
        integer :: dim
        integer :: n_eig
        complex(dp), allocatable :: eigs(:)
        complex(dp), allocatable :: vecs(:,:)
    contains
        procedure :: init => init_eig_results
        procedure :: cleanup => cleanup_eig_results
    end type eig_results

    type(CSR_matrix) :: H_i
    type(PARDISO_solver) :: solver

    type(eig_results), allocatable :: results(:)
contains

    subroutine init_eig_results(this,dim,n_eig)
        class(eig_results), intent(out) :: this
        integer, intent(in) :: dim
        integer, intent(in) :: n_eig

        this%dim = dim
        this%n_eig = n_eig

        allocate(this%eigs(n_eig),this%vecs(dim,n_eig))
    end subroutine init_eig_results

    subroutine cleanup_eig_results(this)
        class(eig_results), intent(inout) :: this

        if (allocated(this%eigs)) deallocate(this%eigs)
        if (allocated(this%vecs)) deallocate(this%vecs)
    end subroutine cleanup_eig_results

    subroutine do_diag(H_block,S_block,bas,res_dir)
        type(block_diag_CS), intent(in) :: H_block
        type(block_diag_CS), intent(in) :: S_block
        type(basis), intent(in) :: bas
        character(len=*), intent(in) :: res_dir

        integer :: i,j

        allocate(results(n_syms))

        do i = 1,n_syms
            j = blocks(i)
            call results(i)%init(bas%syms(j)%n_config,n_eigs(i))
            call diag_block(H_block%blocks(j),S_block%blocks(j),target_energies(i),n_eigs(i),results(i)%eigs,results(i)%vecs)
        end do

        call write_res(res_dir,bas)
    end subroutine do_diag

    subroutine do_diag_essential(H,S,target,vec,eig)
        type(CSR_matrix), intent(in) :: H
        type(CSR_matrix), intent(in) :: S
        complex(dp), intent(in) :: target
        complex(dp), allocatable, intent(inout) :: vec(:)
        complex(dp), intent(out) :: eig

        complex(dp) :: temp_eig(2)
        complex(dp), allocatable :: temp_vec(:,:)
        integer, parameter :: n_energies = 2

        allocate(vec(H%shape(1)), temp_vec(H%shape(1),n_energies))

        call diag_block(H,S,target,n_energies,temp_eig,temp_vec)

        eig = temp_eig(1)
        vec = temp_vec(:,1)
    end subroutine do_diag_essential

    subroutine diag_block(H,S,target,n_energies,energies,vectors)
        type(CSR_matrix), intent(in) :: H
        type(CSR_matrix), intent(in) :: S
        complex(dp), intent(in) :: target
        integer, intent(in) :: n_energies
        complex(dp), intent(out) :: energies(:)
        complex(dp), intent(out) :: vectors(:,:)

        H_i = H
        call H_i%shift_B(-target,S)

        call solver%setup(H_i%shape(1),H_i%nnz,H_i%data,H_i%index_ptr,H_i%indices)
        call solver%factor(H_i%data,H_i%index_ptr,H_i%indices)

        call drive_ARPACK_SI(H_i, solver, S, full, target, n_energies, energies, vectors)
        call sort_eig(n_energies,energies,vectors)

        call solver%cleanup()
    end subroutine diag_block

    subroutine write_res(res_dir,bas)
        character(len=*), intent(in) :: res_dir
        type(basis), intent(in) :: bas

        integer :: i,j,unit,ptr_1,ptr_2
        complex(dp), allocatable :: vec(:,:)
        open(file = res_dir//"energies.out", newunit = unit, action = 'write')

        do i = 1,n_syms
            j = blocks(i)
            write(unit,'(i4,i4,i4,l4)') j,bas%syms(j)%l,bas%syms(j)%m,bas%syms(j)%pi
            do j = 1,results(i)%n_eig
                write(unit,'(a)', advance='no') '('

                if (real(results(i)%eigs(j),kind = dp)<0) then
                    write(unit,'(es24.17e2)', advance='no') real(results(i)%eigs(j),kind = dp)
                else
                    write(unit,'(es23.17e2)', advance='no') real(results(i)%eigs(j),kind = dp)
                end if

                if (aimag(results(i)%eigs(j))<0) then
                    write(unit,'(es24.17e2)', advance='no') aimag(results(i)%eigs(j))
                else
                    write(unit,'(a)', advance='no') '+'
                    write(unit,'(es23.17e2)', advance='no') aimag(results(i)%eigs(j))
                end if

                write(unit,'(a)') 'j)'
            end do
        end do

        close(unit)

        if (store_diag_vecs) then
            allocate(vec(bas%n_states,n_syms),source = (0.0_dp,0.0_dp))
            open(file = res_dir//"vecs.dat", newunit = unit, action = 'write', form = 'unformatted')
            write(unit) bas%n_states
            write(unit) n_syms
            write(unit) 1 !n_calc = 1

            do i = 1,n_syms
                j = blocks(i)
                ptr_1 = bas%sym_ptr(j)
                ptr_2 = bas%sym_ptr(j+1)-1
                vec(ptr_1:ptr_2,i) = results(i)%vecs(:,1)
            end do

            write(unit) vec
            close(unit)
        end if
    end subroutine write_res

end module diagonalization