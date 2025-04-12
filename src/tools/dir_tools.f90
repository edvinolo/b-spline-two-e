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

    subroutine make_res_dir(root_dir_in,res_path)
        character(len=:), allocatable, intent(in) :: root_dir_in
        character(len=:), allocatable, intent(out) ::  res_path

        character(len=:), allocatable :: root_dir
        character(len=4) :: res_dir
        integer :: i
        integer :: length

        length = len(root_dir_in)
        if(root_dir_in(length:length)=='/') then
            root_dir = root_dir_in
        else
            root_dir = root_dir_in // '/'
        end if

        do i = 1,9999
            write(res_dir,'(I4.4)') i
            res_path = root_dir // res_dir // '/'
            if (.not.is_dir(res_path)) then
                call mkdir(res_path)
                write(6,*)
                write(6,*) "Output will be written to ", res_path
                write(6,*)
                exit
            end if
        end do
    end subroutine make_res_dir
end module dir_tools