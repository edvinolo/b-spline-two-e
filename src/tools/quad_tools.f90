module quad_tools
    use stdlib_quadrature, only: gauss_legendre
    implicit none

    type, public :: gau_leg
        integer :: N
        double precision, dimension(:,:), allocatable :: x
        double precision, dimension(:,:), allocatable :: w
    contains
        procedure :: init => init_gau_leg
    end type gau_leg
contains

subroutine setup_GL(N,a,b,x,w)
    integer, intent(in) :: N
    double precision, intent(in) :: a
    double precision, intent(in) :: b
    double precision, dimension(N), intent(inout) :: x
    double precision, dimension(N), intent(inout) :: w

    double precision, dimension(2) :: limits

    limits(1) = a
    limits(2) = b

    call gauss_legendre(x,w,limits)
end subroutine setup_GL

subroutine init_gau_leg(this,breakpoints,N)
    class(gau_leg), intent(inout) :: this
    double precision, dimension(:), intent(in) :: breakpoints
    integer, intent(in) :: N

    integer :: i,j,N_break
    double precision :: a,b
    double precision, dimension(N) :: x_temp,w_temp

    this%N = N

    N_break = size(breakpoints)

    allocate(this%x(N,N_break-1),this%w(N,N_break-1))

    do i = 1, N_break-1
        a = breakpoints(i)
        b = breakpoints(i+1)
        call setup_GL(N,a,b,x_temp,w_temp)
        this%x(:,i) = x_temp
        this%w(:,i) = w_temp
    end do

end subroutine init_gau_leg

end module quad_tools