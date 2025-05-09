module dipole
    use bspline_tools
    use orbital_tools
    use mat_els
    implicit none

contains
    subroutine construct_dip_block_tensor(syms,q,b_splines,S,radial_dip,dip_block,compute)
        type(sym), dimension(2), intent(in) :: syms
        integer, intent(in) :: q
        type(b_spline), intent(in) :: b_splines
        double complex, dimension(:,:), intent(in) ::  S
        double complex, dimension(:,:), intent(in) :: radial_dip
        type(CSR_matrix), intent(out) :: dip_block
        logical, intent(in) :: compute

        double precision :: ang
        logical :: parity_allowed
        integer :: i,j,ptr,nnz
        type(config), dimension(2) :: confs,confs_ex

        ! Check if transitions are allowed between symmetry blocks
        parity_allowed = syms(1)%pi.neqv.syms(2)%pi
        ! Compute the Wigner-Eckart angular factor
        ang = (-1)**(syms(1)%l-syms(1)%m)*three_j(syms(1)%l,1,syms(2)%l,-syms(1)%m,q,syms(2)%m)
        if ((abs(ang)<5.d-16).or.(.not.parity_allowed).or.(.not.compute)) then
            nnz = 0
            call dip_block%init([syms(1)%n_config,syms(2)%n_config],nnz)
            return
        end if

        call init_dip_block(syms,b_splines,dip_block)

        !!$omp parallel do private(confs,confs_ex,ptr,j)
        do i=1,syms(1)%n_config
            confs(1) = syms(1)%configs(i)
            confs_ex(1) = confs(1)
            do ptr = dip_block%index_ptr(i),dip_block%index_ptr(i+1)-1
                j = dip_block%indices(ptr)
                confs(2) = syms(2)%configs(j)
                confs_ex(2)%l = [confs(2)%l(2),confs(2)%l(1)]
                confs_ex(2)%n = [confs(2)%n(2),confs(2)%n(1)]
                dip_block%data(ptr) = ang*dip_mat_neq(syms,confs,S,radial_dip)
            end do
        end do
        !!$omp end parallel do
    end subroutine construct_dip_block_tensor

    subroutine init_dip_block(syms,b_splines,dip_block)
        type(sym), dimension(2), intent(in) :: syms
        type(b_spline), intent(in) :: b_splines
        class(CS_matrix), intent(inout) :: dip_block

        integer :: i,j,nnz, ptr
        type(config), dimension(2) :: confs,confs_ex
        logical :: support,support_ex, nz

        nnz = 0
        do i=1,syms(1)%n_config
            confs(1) = syms(1)%configs(i)
            confs_ex(1) = confs(1)
            do j = 1,syms(2)%n_config
                confs(2) = syms(2)%configs(j)
                confs_ex(2)%l = [confs(2)%l(2),confs(2)%l(1)]
                confs_ex(2)%n = [confs(2)%n(2),confs(2)%n(1)]

                support = (abs(confs(1)%n(1)-confs(2)%n(1))<b_splines%k).and.(abs(confs(1)%n(2)-confs(2)%n(2))<b_splines%k)
                support_ex = (abs(confs_ex(1)%n(1)-confs_ex(2)%n(1))<b_splines%k).and.&
                            (abs(confs_ex(1)%n(2)-confs_ex(2)%n(2))<b_splines%k)

                if (support) then
                    if (ang_dip_red(syms(1)%l,syms(2)%l,confs) > 5.d-16) nnz = nnz + 1
                else if (support_ex) then
                    if (ang_dip_red(syms(1)%l,syms(2)%l,confs) > 5.d-16) nnz = nnz + 1
                end if
            end do
        end do

        call dip_block%init([syms(1)%n_config,syms(2)%n_config],nnz)

        ptr = 1
        do i=1,syms(1)%n_config
            confs(1) = syms(1)%configs(i)
            confs_ex(1) = confs(1)
            dip_block%index_ptr(i) = ptr
            do j = 1,syms(2)%n_config
                confs(2) = syms(2)%configs(j)
                confs_ex(2)%l = [confs(2)%l(2),confs(2)%l(1)]
                confs_ex(2)%n = [confs(2)%n(2),confs(2)%n(1)]

                support = (abs(confs(1)%n(1)-confs(2)%n(1))<b_splines%k).and.(abs(confs(1)%n(2)-confs(2)%n(2))<b_splines%k)
                support_ex = (abs(confs_ex(1)%n(1)-confs_ex(2)%n(1))<b_splines%k).and.&
                            (abs(confs_ex(1)%n(2)-confs_ex(2)%n(2))<b_splines%k)
                nz = .false.
                if (support) then
                    if (ang_dip_red(syms(1)%l,syms(2)%l,confs) > 5.d-16) nz = .true.
                else if (support_ex) then
                    if (ang_dip_red(syms(1)%l,syms(2)%l,confs) > 5.d-16) nz = .true.
                end if

                if (nz) then
                    dip_block%indices(ptr) = j
                    ptr = ptr + 1
                end if
            end do
        end do
        dip_block%index_ptr(syms(1)%n_config+1) = ptr
    end subroutine init_dip_block

    subroutine construct_dip_block_dense(syms,q,S,radial_dip,dip_block_dense)
        type(sym), dimension(2), intent(in) :: syms
        integer, intent(in) :: q
        double complex, dimension(:,:), intent(in) ::  S
        double complex, dimension(:,:), intent(in) :: radial_dip
        double complex, dimension(:,:), intent(out) :: dip_block_dense

        double precision :: ang
        logical :: parity_allowed
        integer :: i,j
        type(config), dimension(2) :: confs


        ! Check if transitions are allowed between symmetry blocks
        parity_allowed = syms(1)%pi.neqv.syms(2)%pi
        ! Compute the Wigner-Eckart angular factor
        ang = (-1)**(syms(1)%l-syms(1)%m)*three_j(syms(1)%l,1,syms(2)%l,-syms(1)%m,q,syms(2)%m)
        if ((abs(ang)<5.d-16).or.(.not.parity_allowed)) then
            dip_block_dense = 0.d0
            return
        end if

        dip_block_dense = 0.d0
        do i=1,syms(1)%n_config
            confs(1) = syms(1)%configs(i)
            do j = 1,syms(2)%n_config
                confs(2) = syms(2)%configs(j)
                dip_block_dense(i,j) = ang*dip_mat_neq(syms,confs,S,radial_dip)
            end do
        end do
    end subroutine construct_dip_block_dense
end module dipole