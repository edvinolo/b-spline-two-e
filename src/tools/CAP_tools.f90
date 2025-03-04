module CAP_tools
    implicit none

    type, public :: CAP
        integer :: order
        double precision :: r_0
        double complex :: eta
    contains
        procedure :: init => init_CAP
        procedure :: V => V_CAP
    end type CAP
contains
    subroutine init_CAP(this,order,r_0,eta)
        class(CAP), intent(inout) :: this
        integer, intent(in) :: order
        double precision, intent(in) :: r_0
        double complex, intent(in) :: eta

        this%order = order
        this%r_0 = r_0
        this%eta = eta
    end subroutine init_CAP

    pure function V_CAP(this,r) result(retval)
        class(CAP), intent(in) :: this
        double precision, intent(in) :: r
        double complex :: retval

        if (r>=this%r_0) then
            retval = dcmplx(0.d0,-1.d0)*this%eta*(r-this%r_0)**this%order
        else
            retval = 0.d0
        end if
    end function V_cap
end module CAP_tools