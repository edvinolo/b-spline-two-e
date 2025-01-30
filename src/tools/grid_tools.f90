module grid_tools
    implicit none

contains

subroutine generate_grid(k,m,Z,h_max,r_max,grid)
    integer, intent(in) :: k,m,Z
    double precision, intent(in) :: h_max,r_max
    double precision,dimension(:),allocatable, intent(inout) :: grid

    double precision, dimension(:), allocatable :: temp_grid

    integer :: i,N_grid
    double precision :: h,next

    ! Minumum grid spacing
    h = 2.d0**(-m)

    ! Compute size of grid array (too big reduce after filling)
    N_grid = 2*k + m + ceiling(r_max/h)
    allocate(grid(N_grid))

    ! 
    grid(1:k) = 0.d0

    ! Small linear spacing
    do i = k+1,k+m
        grid(i) = grid(i-1) + h
    end do

    ! increase grid spacing until it reaches h_max
    i=k+m
    do
        next = grid(i)*(1.d0 + h)
        if ((next-grid(i)).ge.h_max) exit 
        grid(i+1) = next
        i = i + 1
    end do

    ! Linear spacing with h_max
    do while(grid(i) < z*r_max)
        grid(i+1) = grid(i) + h_max
        i = i + 1
    end do

    ! 
    grid(i+1:i+k-1) = grid(i)

    N_grid = i+k-1

    call move_alloc(grid,temp_grid)

    grid = temp_grid(1:N_grid)

end subroutine generate_grid

end module grid_tools