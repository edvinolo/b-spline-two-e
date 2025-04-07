program test_dir
    use dir_tools, only: make_res_dir
    implicit none

    character(len=:), allocatable :: res_dir,root_dir

    root_dir = 'TESTOUTPUT/'

    call make_res_dir(root_dir,res_dir)
end program test_dir