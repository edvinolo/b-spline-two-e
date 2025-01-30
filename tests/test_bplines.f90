program test_bsplines
    use grid_tools
    use bspline_tools
    implicit none

    integer :: k,m,Z,i,j
    double precision :: h_max,r_max
    double precision :: dr
    double precision, dimension(:), allocatable :: grid,r,test_grid
    double complex,dimension(:),allocatable :: c
    double complex :: val
    double precision,dimension(:),allocatable :: vals

    type(b_spline) :: splines,spline_test

    k = 5
    m = 3
    Z = 1
    h_max = 0.5d0
    r_max = 10.d0
    call generate_grid(k,m,Z,h_max,r_max,grid)



    call splines%init(k,grid)

    write(6,*) splines%knots
    write(6,*) splines%breakpoints

    call INTERV(grid,size(grid),2.d0,z,m)
    write(6,*) z,m
    write(6,*) grid(z),grid(z+1)

    call INTERV(grid,size(grid),1.d0,z,m)

    write(6,*) z,m
    write(6,*) grid(z),grid(z+1)

    allocate(c(splines%n),r(1101))

    r = 0.d0
    dr = 0.01d0
    do i = 2,1101
        r(i) = r(i-1) + dr
    end do

    c = dcmplx(0.d0,0.d0)
    c(10) = dcmplx(1.d0,1.d0)

    open(unit= 1, file = 'spline.dat')
    do i = 1,1101
        val = splines%eval_z(r(i),c)
        write(1,*) r(i),real(val,kind=8),aimag(val)
    end do
    close(1)

    allocate(test_grid(1:19))
    test_grid = 0.d0
    do i = 6,15
        test_grid(i) = test_grid(i-1) + 0.1d0
    end do
    test_grid(16:19) = test_grid(15)

    deallocate(r)
    allocate(r(1001))

    r = 0.d0
    dr = 0.001d0
    do i = 2,1001
        r(i) = r(i-1) + dr
    end do

    call spline_test%init(5,test_grid)
    deallocate(c)
    allocate(c(spline_test%n))
    c = dcmplx(0.d0,0.d0)

    open(unit= 1, file = 'spline_test.dat')
    allocate(vals(spline_test%n))
    do i = 1,1001
        do j= 1,spline_test%n
            c = dcmplx(0.d0,0.d0)
            c(j) = 1.d0
            vals(j) = real(spline_test%eval_z(r(i),c),kind=8)
        end do
        write(1,*) r(i),vals
    end do
    close(1)

end program test_bsplines