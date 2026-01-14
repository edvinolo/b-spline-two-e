module function_eval
    use kind_tools
    use bspline_tools
    use block_tools
    use orbital_tools
    use wigner_tools, only: CG,sph_harm
    use mat_els, only: setup_S
    use input_tools, only: k_GL, N_theta, N_phi, N_r, theta_limits, phi_limits, r_limits
    use stdlib_math, only: linspace
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    implicit none

    real(dp), allocatable :: r_vec(:)
    real(dp), allocatable :: theta_vec(:)
    real(dp), allocatable :: phi_vec(:)
    complex(dp), allocatable :: S(:,:) ! If ECS is not used then this could be real
    real(dp), allocatable :: density(:,:,:)
    real(dp), allocatable :: density_mat(:,:,:,:,:) !Should probably not use a regular array here, since it will be ginourmus.
    real(dp), allocatable :: R_n(:,:)
    logical, allocatable :: R_n_support(:,:)
    real(dp), allocatable :: rad_dens(:,:,:)
    complex(dp), allocatable :: Y_lm(:,:,:,:)
    complex(dp), allocatable :: Y_lm_conjg(:,:,:,:)
    real(dp), allocatable :: CG_table(:,:,:,:,:,:) ! Use hashmap instead? (probably not necessary)
    real(dp), allocatable :: ang_dens(:,:,:,:,:,:,:,:,:)
contains

    subroutine eval_from_states_1p(res_dir,root_dir,bas,splines)
        character(len=*), intent(in) :: res_dir
        character(len=*), intent(in) :: root_dir
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines

        complex(dp), allocatable :: vecs(:,:,:)
        integer :: unit,n_states,n_quasi,i,j,index
        integer (int32) n_calc
        character(len=:), allocatable :: res_file
        character(len=4) :: file

        open(file = root_dir//"vecs.dat", newunit = unit, action = 'read', form = 'unformatted')
        read(unit) n_states
        read(unit) n_quasi
        read(unit) n_calc

        allocate(vecs(n_states,n_quasi,n_calc),source=(0.0_dp,0.0_dp))
        read(unit) vecs

        index = 1
        do j = 1,n_calc
            do i = 1,n_quasi
                write(file,'(I4.4)') index
                res_file = res_dir//'rad_'//file//'.out'
                call write_func_1p(vecs(:,i,j),bas,splines,res_file)
                index = index + 1
            end do
        end do

    end subroutine eval_from_states_1p

    subroutine eval_from_states_2p(res_dir,root_dir,bas,splines)
        character(len=*), intent(in) :: res_dir
        character(len=*), intent(in) :: root_dir
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines

        complex(dp), allocatable :: vecs(:,:,:)
        integer :: unit,n_states,n_quasi,i,j,index
        integer (int32) n_calc
        character(len=:), allocatable :: res_file
        character(len=4) :: file

        open(file = root_dir//"vecs.dat", newunit = unit, action = 'read', form = 'unformatted')
        read(unit) n_states
        read(unit) n_quasi
        read(unit) n_calc

        allocate(vecs(n_states,n_quasi,n_calc),source=(0.0_dp,0.0_dp))
        read(unit) vecs


        call setup_two_el_density(splines,bas)

        index = 1
        do j = 1,n_calc
            do i = 1,n_quasi
                write(file,'(I4.4)') index
                res_file = res_dir//'density_'//file//'.out'
                call write_density_2p(vecs(:,i,j),bas,splines,res_file)
                index = index + 1
            end do
        end do

    end subroutine eval_from_states_2p

    subroutine write_func_1p(state,bas,splines,res_file)
        complex(dp), intent(in) :: state(:)
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines
        character(len=*), intent(in) :: res_file

        complex(dp), allocatable :: coeffs(:)
        complex(dp) :: val

        integer :: i,j,ptr,i_r,unit,index

        if (.not.allocated(r_vec)) allocate(r_vec(N_r))
        r_vec = linspace(r_limits(1),r_limits(2),N_r)
        allocate(coeffs(splines%n))


        open(file = res_file, newunit = unit, action = 'write')
        write(unit,'(I4,I4)') N_r,bas%n_sym
        do i_r = 1,N_r
            write(unit,'(es24.17e3)') r_vec(i_r)
        end do

        ptr = 1
        do i = 1,bas%n_sym
            write(unit,'()')
            write(unit,'(I4,I4)') bas%syms(i)%l,bas%syms(i)%m

            coeffs = (0.0_dp,0.0_dp)
            do j = 1,bas%syms(i)%n_config
                index = bas%syms(i)%configs(j)%n(1) + 1
                coeffs(index) = state(ptr)
                ptr = ptr + 1
            end do

            do i_r = 1,N_r
                val = eval_rad_1p(coeffs,splines,r_vec(i_r))
                write(unit,'(a)', advance='no') '('

                if (real(val,kind = dp)<0) then
                    write(unit,'(es25.17e3)', advance='no') real(val,kind = dp)
                else
                    write(unit,'(es24.17e3)', advance='no') real(val,kind = dp)
                end if

                if (aimag(val)<0) then
                    write(unit,'(es25.17e3)', advance='no') aimag(val)
                else
                    write(unit,'(a)', advance='no') '+'
                    write(unit,'(es24.17e3)', advance='no') aimag(val)
                end if

                write(unit,'(a)') 'j) '
            end do

        end do

        close(unit)
    end subroutine write_func_1p

    subroutine write_func_1p_from_vals(vals,bas,res_file,lm_basis)
        complex(dp), intent(in) :: vals(:,:)
        type(basis), intent(in) :: bas
        character(len=*), intent(in) :: res_file
        logical, intent(in) :: lm_basis

        integer :: i,i_r,unit
        complex(dp) :: val

        open(file = res_file, newunit = unit, action = 'write')
        write(unit,'(a)') "# First row gives number of radial points (N_r) and number of channels (n_sym)."
        write(unit,'(a)') "# The next N_r rows gives the radial points. Then follows a blank row."
        write(unit,'(a)') "# Then comes the radial function of each channel. First row is the channel labels"
        write(unit,'(a)') "# and the next N_r rows are the values. The channels are separated by blank rows."

        write(unit,'(I4,I4)') N_r,bas%n_sym
        do i_r = 1,N_r
            write(unit,'(es24.17e3)') r_vec(i_r)
        end do

        do i = 1,bas%n_sym
            write(unit,'()')
            if (lm_basis) then
                write(unit,'(I4,I4)') bas%syms(i)%l,bas%syms(i)%m
            else
                write(unit,'(I4,I4)') i,i
            end if

            do i_r = 1,N_r
                val = vals(i_r,i)
                if (isnan(real(val,kind=dp)).or.isnan(aimag(val))) then
                    write(stdout,*) val,r_vec(i_r),i
                end if
                write(unit,'(a)', advance='no') '('

                if (real(val,kind = dp)<0) then
                    write(unit,'(es25.17e3)', advance='no') real(val,kind = dp)
                else
                    write(unit,'(es24.17e3)', advance='no') real(val,kind = dp)
                end if

                if (aimag(val)<0) then
                    write(unit,'(es25.17e3)', advance='no') aimag(val)
                else
                    write(unit,'(a)', advance='no') '+'
                    write(unit,'(es24.17e3)', advance='no') aimag(val)
                end if

                write(unit,'(a)') 'j) '
            end do

        end do

        close(unit)
    end subroutine write_func_1p_from_vals

    subroutine array_func_1p(state,bas,splines,arr_out)
        complex(dp), intent(in) :: state(:)
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines
        complex(dp), intent(out) :: arr_out(:,:)

        complex(dp), allocatable :: coeffs(:)
        complex(dp) :: val

        integer :: i,j,ptr,i_r,index

        if (.not.allocated(r_vec)) allocate(r_vec(N_r))
        r_vec = linspace(r_limits(1),r_limits(2),N_r)
        allocate(coeffs(splines%n))

        ptr = 1
        do i = 1,bas%n_sym

            coeffs = (0.0_dp,0.0_dp)
            do j = 1,bas%syms(i)%n_config
                index = bas%syms(i)%configs(j)%n(1) + 1
                coeffs(index) = state(ptr)
                ptr = ptr + 1
            end do

            do i_r = 1,N_r
                val = eval_rad_1p(coeffs,splines,r_vec(i_r))
                arr_out(i_r,i) = val
            end do

        end do
    end subroutine array_func_1p

    subroutine write_density_2p(state,bas,splines,res_file)
        complex(dp), intent(in) :: state(:)
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines
        character(len=*), intent(in) :: res_file

        integer :: i_r,i_theta,i_phi
        integer :: unit

        do i_phi = 1,N_phi
            do i_theta = 1,N_theta
                !$omp parallel do schedule(dynamic)
                do i_r = 1,N_r
                    write(stdout,'(a,i0,a,i0,a,i0)') 'i_phi: ',i_phi,' i_theta: ',i_theta,' i_r: ',i_r
                    density(i_r,i_theta,i_phi) = get_density_point_2p(state,bas,splines,i_r,i_theta,i_phi)
                end do
                !$omp end parallel do
            end do
        end do


        open(file = res_file, newunit = unit, action = 'write')
        write(unit,'(a)') "# First line is N_r,N_theta,N_phi. Following lines are r,theta,phi,density"
        write(unit,'(I4,I4,I4)') N_r,N_theta,N_phi

        do i_phi = 1,N_phi
            do i_theta = 1,N_theta
                do i_r = 1,N_r
                    write(unit,'(es24.17e3,a)',advance='no') r_vec(i_r), ' '
                    write(unit,'(es24.17e3,a)',advance='no') theta_vec(i_theta), ' '
                    write(unit,'(es24.17e3,a)',advance='no') phi_vec(i_phi), ' '
                    write(unit,'(es24.17e3)') density(i_r,i_theta,i_phi)
                end do
            end do
        end do

        close(unit)
    end subroutine write_density_2p

    function eval_rad_1p(coeffs,splines,r) result(res)
        complex(dp), intent(in) :: coeffs(:)
        type(b_spline), intent(in) :: splines
        real(dp), intent(in) :: r
        complex(dp) :: res

        real(dp), allocatable :: re_coeffs(:),im_coeffs(:)
        real(dp) :: re_res,im_res

        allocate(re_coeffs(size(coeffs)),im_coeffs(size(coeffs)))
        re_coeffs = real(coeffs,kind = dp)
        im_coeffs = aimag(coeffs)

        re_res = splines%eval_d(r,re_coeffs)
        im_res = splines%eval_d(r,im_coeffs)

        res = cmplx(re_res,im_res,kind = dp)
    end function eval_rad_1p

    subroutine setup_two_el_density(splines,bas)
        type(b_spline), intent(in) :: splines
        type(basis), intent(in) :: bas

        allocate(r_vec(N_r),theta_vec(N_theta),phi_vec(N_phi))
        r_vec = linspace(r_limits(1),r_limits(2),N_r)
        theta_vec = linspace(theta_limits(1),theta_limits(2),N_theta)
        phi_vec = linspace(phi_limits(1),phi_limits(2),N_phi)

        allocate(density(N_r,N_theta,N_phi),source=0.0_dp)
        allocate(S(splines%n_b,splines%n_b),source=(0.0_dp,0.0_dp))

        write(stdout,*) "Setting up two-electron density variables..."

        call setup_S(splines,k_GL,S)
        call setup_R_n(splines)
        call setup_rad_dens(splines)
        call setup_Y_lm(bas)
        call setup_CG_table(bas)
        call setup_ang_dens(bas)

        write(stdout,*) "Done with two-electron density setup!"

    end subroutine setup_two_el_density

    subroutine setup_R_n(splines)
        type(b_spline), intent(in) :: splines

        integer :: i,j,iv
        real(dp), allocatable :: coeffs(:)

        allocate(R_n(splines%n_b,N_r),source=0.0_dp)
        allocate(R_n_support(splines%n_b,N_r),source=.false.)
        allocate(coeffs(splines%n),source = 0.0_dp)

        do i = 1,N_r
            iv = splines%find_interval(r_vec(i))
            do j = 1,splines%n_b
                coeffs(j+1) = 1.0_dp
                R_n(j,i) = splines%eval_d(r_vec(i),coeffs)/r_vec(i)
                coeffs(j+1) = 0.0_dp
                R_n_support(j,i) = splines%support(iv,j+1)
            end do
        end do
    end subroutine setup_R_n

    subroutine setup_Y_lm(bas)
        type(basis), intent(in) :: bas

        integer :: i,j,l_i,m_i
        integer :: l_1p

        l_1p = bas%max_l_1p

        allocate(Y_lm(-l_1p:l_1p,0:l_1p,N_theta,N_phi),source=(0.0_dp,0.0_dp))
        allocate(Y_lm_conjg(-l_1p:l_1p,0:l_1p,N_theta,N_phi),source=(0.0_dp,0.0_dp))

        do i = 1,N_phi
            do j = 1,N_theta
                do l_i = 0,l_1p
                    do m_i = -l_i,l_i
                        Y_lm(m_i,l_i,j,i) = sph_harm(l_i,m_i,theta_vec(j),phi_vec(i))
                        Y_lm_conjg(m_i,l_i,j,i) = conjg(Y_lm(m_i,l_i,j,i))
                    end do
                end do
            end do
        end do
    end subroutine setup_Y_lm

    subroutine setup_CG_table(bas)
        type(basis), intent(in) :: bas

        integer :: i,l_i,m_i,l_1,l_2,m_1,m_2
        integer :: l_1p,l_2p

        l_1p = bas%max_l_1p
        l_2p = bas%max_L

        allocate(CG_table(-l_1p:l_1p,-l_1p:l_1p,0:l_1p,0:l_1p,-l_2p:l_2p,0:l_2p),source = 0.0_dp)

        do i = 1,bas%n_sym
            l_i = bas%syms(i)%l
            m_i = bas%syms(i)%m
            do l_2 = 0,l_1p
                do l_1 = 0,l_1p
                    do m_2 = -l_2,l_2
                        do m_1 = -l_1,l_1
                            CG_table(m_1,m_2,l_1,l_2,m_i,l_i) = CG(l_1,l_2,l_i,m_1,m_2,m_i)
                        end do
                    end do
                end do
            end do
        end do

    end subroutine setup_CG_table

    subroutine setup_ang_dens(bas)
        type(basis), intent(in) :: bas

        integer :: l_1p,l_2p
        integer :: l_a_1,l_a_2,l_b_1,M_a,L_a,M_b,L_b
        integer :: i_theta,i_phi

        l_1p = bas%max_l_1p
        l_2p = bas%max_L

        allocate(ang_dens(0:l_1p,0:l_1p,0:l_1p,-l_2p:l_2p,0:l_2p,-l_2p:l_2p,0:l_2p,N_theta,N_phi),source = 0.0_dp)

        do i_phi = 1,N_phi
            do i_theta = 1,N_theta
                do L_b = 0,l_2p
                    do M_b = -L_b,L_b
                        do L_a = 0,l_2p
                            do M_a = -L_a,L_a
                                do l_b_1 = 0,l_1p
                                    do l_a_2 = 0,l_1p
                                        do l_a_1 = 0,l_1p
                                            ang_dens(l_a_1,l_a_2,l_b_1,M_a,L_a,M_b,L_b,i_theta,i_phi) = &
                                            ang_m_sum([l_a_1,l_a_2,l_b_1],[L_a,M_a],[L_b,M_b],i_theta,i_phi)
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine setup_ang_dens

    function ang_m_sum(l,L_a,L_b,i_theta,i_phi) result(res)
        integer, intent(in) :: l(3)
        integer, intent(in) :: L_a(2)
        integer, intent(in) :: L_b(2)
        integer, intent(in) :: i_theta
        integer, intent(in) :: i_phi
        real(dp) :: res

        integer :: m_a_1,m_a_2,m_b_1
        complex(dp) :: temp

        temp = 0.0_dp
        do m_a_1 = -l(1),l(1)
            do m_a_2 = -l(2),l(2)
                do m_b_1 = -l(3),l(3)
                    temp = temp + CG_table(m_a_1,m_a_2,l(1),l(2),L_a(2),L_a(1))&
                            *CG_table(m_b_1,m_a_2,l(3),l(2),L_b(2),L_b(1))&
                            *Y_lm_conjg(m_a_1,l(1),i_theta,i_phi)*Y_lm(m_b_1,l(3),i_theta,i_phi)
                end do
            end do
        end do

        res = real(temp,kind=dp)

    end function ang_m_sum

    subroutine setup_rad_dens(splines)
        type(b_spline), intent(in) :: splines

        integer :: i,j,i_r

        allocate(rad_dens(splines%n_b,splines%n_b,N_r),source = 0.0_dp)

        do i_r = 1,N_r
            do i = 1,splines%n_b
                do j = 1,splines%n_b
                    if (R_n_support(i,i_r).and.R_n_support(j,i_r)) rad_dens(j,i,i_r) = R_n(j,i_r)*R_n(i,i_r)
                end do
            end do
        end do
    end subroutine setup_rad_dens

    function get_density_point_2p(state,bas,splines,i_r,i_theta,i_phi) result(res)
        complex(dp), intent(in) :: state(:)
        type(basis), intent(in) :: bas
        type(b_spline), intent(in) :: splines
        integer, intent(in) :: i_r
        integer, intent(in) :: i_theta
        integer, intent(in) :: i_phi
        real(dp) :: res

        type(sym) :: syms(2)
        type(config) :: confs(2)
        integer :: ptr_j, ptr_i, i, j, k, l
        complex(dp) :: dens
        complex(dp) :: temp

        ! Should perhaps reorganize this to utilize symmetry
        ! Can turn the two inner loops into a function, that return when ptr_i == ptr_j
        ! Taking appropriate care for symmetry in the expressions
        ptr_i = 0
        temp = 0

        do i = 1,bas%n_sym
            syms(2) = bas%syms(i)
            do k = 1,bas%syms(i)%n_config
                confs(2) = bas%syms(i)%configs(k)
                ptr_i = ptr_i + 1
                ptr_j = 0
                do j = 1,bas%n_sym
                    syms(1) = bas%syms(j)
                    do l = 1,bas%syms(j)%n_config
                        confs(1) = bas%syms(j)%configs(l)
                        ptr_j = ptr_j + 1
                        ! if (abs(state(ptr_i)*state(ptr_j)) < 1.0e-16_dp) cycle
                        if (all(confs(1)%l == confs(2)%l).or.all(confs(1)%l == confs(2)%l([2,1]))) then
                            if (all(abs(confs(1)%n - confs(2)%n) < splines%k).or.&
                                    all(abs(confs(1)%n - confs(2)%n([2,1])) < splines%k)) then
                                dens = density_mat_el(syms,confs,i_r,i_theta,i_phi,splines)
                                temp = temp + conjg(state(ptr_i))*state(ptr_j)*dens
                            end if
                        end if
                        ! write(stdout,*) i,k,j,l,dens
                    end do
                end do
            end do
        end do

        res = real(temp,kind=dp)
    end function get_density_point_2p

    function density_mat_el(syms,confs,i_r,i_theta,i_phi,splines) result(res)
        type(sym), intent(in) :: syms(2)
        type(config), intent(in) :: confs(2)
        integer, intent(in) :: i_r
        integer, intent(in) :: i_theta
        integer, intent(in) :: i_phi
        type(b_spline), intent(in) :: splines
        complex(dp) :: res

        type(config) :: confs_temp(2)

        res = 0.0_dp
        confs_temp%eqv = .false.

        res = res + partial_density_mat_el(syms,confs,i_r,i_theta,i_phi,splines)

        confs_temp(1) = confs(1)
        confs_temp(2)%l = confs(2)%l([2,1])
        confs_temp(2)%n = confs(2)%n([2,1])
        res = res + (-1)**(syms(2)%l+sum(confs(2)%l))*partial_density_mat_el(syms,confs_temp,i_r,i_theta,i_phi,splines)

        confs_temp(2) = confs(2)
        confs_temp(1)%l = confs(1)%l([2,1])
        confs_temp(1)%n = confs(1)%n([2,1])
        res = res + (-1)**(syms(1)%l+sum(confs(1)%l))*partial_density_mat_el(syms,confs_temp,i_r,i_theta,i_phi,splines)

        confs_temp(1)%l = confs(1)%l([2,1])
        confs_temp(1)%n = confs(1)%n([2,1])
        confs_temp(2)%l = confs(2)%l([2,1])
        confs_temp(2)%n = confs(2)%n([2,1])
        res = res + (-1)**(syms(1)%l + syms(2)%l + sum(confs(1)%l) + sum(confs(2)%l))&
        *partial_density_mat_el(syms,confs_temp,i_r,i_theta,i_phi,splines)

    end function density_mat_el

    function partial_density_mat_el(syms,confs,i_r,i_theta,i_phi,splines) result(res)
        type(sym), intent(in) :: syms(2)
        type(config), intent(in) :: confs(2)
        integer, intent(in) :: i_r
        integer, intent(in) :: i_theta
        integer, intent(in) :: i_phi
        type(b_spline), intent(in) :: splines
        complex(dp) :: res

        integer :: l_a_1,l_a_2,l_b_1,l_b_2
        integer :: n_a_1,n_a_2,n_b_1,n_b_2

        l_a_1 = confs(1)%l(1)
        l_a_2 = confs(1)%l(2)
        l_b_1 = confs(2)%l(1)
        l_b_2 = confs(2)%l(2)

        n_a_1 = confs(1)%n(1)
        n_a_2 = confs(1)%n(2)
        n_b_1 = confs(2)%n(1)
        n_b_2 = confs(2)%n(2)


        if ((l_a_2 /= l_b_2).or.(abs(n_a_2-n_b_2)>= splines%k)) then
            res = 0.0_dp
            return
        end if

        if (.not.(R_n_support(n_a_1,i_r).and.R_n_support(n_b_1,i_r))) then
            res = 0.0_dp
            return
        end if

        res = ang_dens(l_a_1,l_a_2,l_b_1,syms(1)%m,syms(1)%l,syms(2)%m,syms(2)%l,i_theta,i_phi)&
                *rad_dens(n_a_1,n_b_1,i_r)*S(n_a_2,n_b_2)

        ! write(stdout,*) res, maxval(abs(ang_dens)), maxval(abs(rad_dens))
    end function partial_density_mat_el
end module function_eval