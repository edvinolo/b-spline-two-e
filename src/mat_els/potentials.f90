module potentials
    implicit none

    type, abstract :: sph_pot
    contains
        procedure (spherical_V), deferred, pass :: v
    end type sph_pot

    abstract interface
        function spherical_V(this,r,l) result(retval)
            import :: sph_pot
            class(sph_pot), intent(in) :: this
            double precision, intent(in) :: r
            integer, intent(in) :: l
            double precision :: retval
        end function spherical_V
    end interface

    type, extends(sph_pot), public :: hydrogenic
        integer :: Z
    contains
        procedure :: init => init_hydrogenic
        procedure :: V => V_hydrogenic
    end type hydrogenic

contains

    subroutine init_hydrogenic(this,Z)
        class(hydrogenic), intent(inout) :: this
        integer, intent(in) ::  Z

        this%Z = Z
    end subroutine init_hydrogenic

    function V_hydrogenic(this,r,l) result(retval)
        class(hydrogenic), intent(in) :: this
        double precision, intent(in) :: r
        integer, intent(in) :: l

        double precision :: retval

        retval = 0.5d0*l*(l+1)/r**2 - this%Z/r
    end function V_hydrogenic

end module potentials