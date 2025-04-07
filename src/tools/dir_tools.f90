module dir_tools
    implicit none

contains
    function is_dir(dir) result(res)
        character(len=:), allocatable, intent(in) :: dir
        logical :: res

        inquire(file=dir,exist=res)
    end function is_dir

    subroutine mkdir(dir)
        character(len=:), allocatable, intent(in) :: dir

        call execute_command_line("mkdir " // dir)
    end subroutine mkdir

    subroutine make_res_dir(root_dir,res_path)
        character(len=:), allocatable, intent(in) :: root_dir
        character(len=:), allocatable, intent(out) ::  res_path

        character(len=4) :: res_dir
        integer :: i

        write(6,*) root_dir

        do i = 1,9999
            write(res_dir,'(I4.4)') i
            res_path = root_dir // res_dir
            write(6,*) res_path
            if (.not.is_dir(res_path)) then
                call mkdir(res_path)
                exit
            end if
        end do
    end subroutine make_res_dir
end module dir_tools