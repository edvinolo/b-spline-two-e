module function_eval
    use kind_tools
    use bspline_tools
    use block_tools
    use orbital_tools
    use stdlib_math, only: linspace
    implicit none

contains

    subroutine eval_from_states_1p(res_dir,root_dir,bas,r_limits,N_r,splines)
        character(len=*), intent(in) :: res_dir
        character(len=*), intent(in) :: root_dir
        type(basis), intent(in) :: bas
        real(dp), intent(in) :: r_limits(2)
        integer(int32), intent(in) :: N_r
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
                call write_func_1p(vecs(:,i,j),bas,r_limits,N_r,splines,res_file)
                index = index + 1
            end do
        end do

    end subroutine eval_from_states_1p

    subroutine write_func_1p(state,bas,r_limits,N_r,splines,res_file)
        complex(dp), intent(in) :: state(:)
        type(basis), intent(in) :: bas
        real(dp), intent(in) :: r_limits(2)
        integer(int32), intent(in) :: N_r
        type(b_spline), intent(in) :: splines
        character(len=*), intent(in) :: res_file

        real(dp), allocatable :: r_vec(:)
        complex(dp), allocatable :: coeffs(:)
        complex(dp) :: val

        integer :: i,j,ptr,i_r,unit,index

        r_vec = linspace(r_limits(1),r_limits(2),N_r)
        allocate(coeffs(splines%n))


        open(file = res_file, newunit = unit, action = 'write')
        write(unit,'(I4,I4)') N_r,bas%n_sym
        do i_r = 1,N_r
            write(unit,'(es24.17e3)') r_vec(i_r)
        end do

        ptr = 1
        do i = 1,bas%n_sym
            write(unit,*)''
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

                write(unit,'(a)', advance='no') 'j) '
                write(unit,*)''
            end do

        end do

        close(unit)
    end subroutine write_func_1p

    function eval_rad_1p(coeffs,splines,r) result(res)
        complex(dp), intent(in) :: coeffs(:)
        type(b_spline), intent(in) :: splines
        real(dp), intent(in) :: r
        complex(dp) :: res

        real(dp), allocatable :: re_coeffs(:),im_coeffs(:)
        real(dp) :: re_res,im_res

        re_coeffs = real(coeffs,kind = dp)
        im_coeffs = aimag(coeffs)

        re_res = splines%eval_d(r,re_coeffs)
        im_res = splines%eval_d(r,im_coeffs)

        res = cmplx(re_res,im_res,kind = dp)
    end function eval_rad_1p
end module function_eval