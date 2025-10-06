module fields
    use kind_tools
    use constants_tools
    use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
    implicit none

    type, public :: pulse
        real(dp), allocatable :: env_params(:)
        real(dp), allocatable :: carrier_params(:)
        real(dp) :: A_0
        logical :: pol(3)
        procedure(envelope), pointer :: env => null()
        procedure(carrier), pointer :: carr => null()
    contains
        procedure :: init => init_pulse
        procedure :: eval_A
    end type pulse

    abstract interface
        function envelope(this,t) result(res)
            import :: dp
            import :: pulse
            class(pulse), intent(in) :: this
            real(dp), intent(in) :: t
            real(dp) :: res
        end function envelope

        function carrier(this,t) result(res)
            import :: dp
            import :: pulse
            class(pulse), intent(in) :: this
            real(dp), intent(in) :: t
            real(dp) :: res
        end function carrier
    end interface

contains
    subroutine init_pulse(this,A_0,t_0,env_params,carrier_params,pol,env_name,carrier_name)
        class(pulse), intent(out) :: this
        real(dp), intent(in) :: A_0
        real(dp), intent(in) :: t_0
        real(dp), intent(in) :: env_params(:)
        real(dp), intent(in) :: carrier_params(:)
        logical, intent(in) :: pol(3)
        character(len=*), intent(in) :: env_name
        character(len=*), intent(in) :: carrier_name

        this%env_params = [t_0,env_params]
        this%carrier_params = [t_0,carrier_params]
        this%A_0 = A_0
        this%pol = pol

        if (env_name == 'cos2') then
            this%env => cos2_env
        else if (env_name == 'gaus') then
            this%env => gaus_env
        else if (env_name == 'gaus_chirp') then
            this%env => gaus_chirp_env
            if (carrier_name /= "sin_lin_chirp") then
                write(stderr,*) "Error, gaus_chirp envelope must be used with sin_lin_chirp carrier."
                write(stderr,*) "Envelope: ", env_name
                write(stderr,*) "Carrier: ", carrier_name
                stop
            end if
        else if (env_name == 'super_gaus') then
            this%env => super_gaus_env
        else if (env_name == 'flat_top') then
            this%env => flat_top_env
        else
            write(stderr,*) "Error, unkown envelope type: ", env_name
            stop
        end if

        if (carrier_name == 'sin') then
            this%carr => sin_carrier
        else if(carrier_name == 'sin_lin_chirp') then
            this%carr => sin_lin_chirp_carrier
            if (env_name /= "gaus_chirp") then
                write(stderr,*) "Error, gaus_chirp envelope must be used with sin_lin_chirp carrier."
                write(stderr,*) "Envelope: ", env_name
                write(stderr,*) "Carrier: ", carrier_name
                stop
            end if
        else
            write(stderr,*) "Error, unkown carrier type: ", carrier_name
            stop
        end if

    end subroutine init_pulse

    function eval_A(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        res = this%A_0*this%carr(t)*this%env(t)
    end function

    ! ---------------------------------------------------------------
    ! Carriers
    ! ---------------------------------------------------------------
    function sin_carrier(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%carrier_params(1),&
                omega => this%carrier_params(2),&
                phi => this%carrier_params(3)&
                )

        res = sin(omega*(t-t_0)-phi)
        end associate
    end function sin_carrier

    function sin_lin_chirp_carrier(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        real(dp) :: omega_t

        associate(t_0 => this%carrier_params(1),&
            omega => this%carrier_params(2),&
            phi => this%carrier_params(3),&
            FWHM => this%env_params(2),&
            beta => this%env_params(3)&
        )
        omega_t = omega + 4.0_dp*log_2/(beta + 1.0_dp/beta)*(t-t_0)/FWHM**2
        res = sin(omega_t*(t-t_0) - phi)
        end associate
    end function

    ! ---------------------------------------------------------------
    ! Envelopes
    ! ---------------------------------------------------------------
    function cos2_env(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                tau => this%env_params(2)&
                )

        if (abs(t_0-t) < 0.5_dp*tau) then
            res = cos(pi*(t-t_0)/tau)**2
        else
            res = 0
        end if
        end associate
    end function cos2_env

    function gaus_env(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                FWHM => this%env_params(2)&
                )
            res = exp(-half_log_2*(2*(t-t_0)/FWHM)**2)
        end associate

    end function gaus_env

    function gaus_chirp_env(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        real(dp) :: FWHM_beta,amp_beta

        associate(t_0 => this%env_params(1),&
                FWHM => this%env_params(2),&
                beta => this%env_params(3)&
        )

        FWHM_beta = sqrt(1.0_dp + beta**2)*FWHM
        amp_beta = (1.0_dp + beta**2)**(-1.0_dp/4.0_dp)
        res = amp_beta*exp(-half_log_2*(2*(t-t_0)/FWHM_beta)**2)
        end associate
    end function gaus_chirp_env

    function super_gaus_env(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                FWHM => this%env_params(2),&
                n => this%env_params(3)&
                )
            res = exp(-half_log_2*(2*(t-t_0)/FWHM)**(2*n))
        end associate
    end function super_gaus_env

    function flat_top_env(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                tau => this%env_params(2)&
        )
            if (abs(t-t_0)<0.5_dp*tau) then
                res = 1.0_dp
            else
                res = 0.0_dp
            end if
        end associate
    end function flat_top_env
end module fields