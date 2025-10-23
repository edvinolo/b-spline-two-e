module fields
    use kind_tools
    use constants_tools
    use iso_fortran_env, only: stderr => error_unit, stdout => output_unit
    implicit none

    type, public :: pulse
        character(len=:), allocatable :: env_name
        real(dp), allocatable :: env_params(:)
        character(len=:), allocatable :: carrier_name
        real(dp), allocatable :: carrier_params(:)
        real(dp) :: intensity
        real(dp) :: A_0
        real(dp) :: E_0
        real(dp) :: U_p
        logical :: pol(3)
        integer :: id
        procedure(envelope), pointer :: env => null()
        procedure(carrier), pointer :: carr => null()
        procedure(envelope), pointer :: env_d => null()
        procedure(carrier), pointer :: carr_d => null()
    contains
        procedure :: init => init_pulse
        procedure :: print => print_pulse
        procedure :: eval_A
        procedure :: eval_E
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
    subroutine init_pulse(this,intensity,t_0,env_params,carrier_params,pol,env_name,carrier_name,id)
        class(pulse), intent(out) :: this
        real(dp), intent(in) :: intensity
        real(dp), intent(in) :: t_0
        real(dp), intent(in) :: env_params(:)
        real(dp), intent(in) :: carrier_params(:)
        character(len=1), intent(in) :: pol
        character(len=*), intent(in) :: env_name
        character(len=*), intent(in) :: carrier_name
        integer, intent(in) :: id

        this%env_params = [t_0,env_params]
        this%carrier_params = [t_0,carrier_params]
        this%intensity = intensity
        this%A_0 = A_0_au(this%intensity,this%carrier_params(2))
        this%E_0 = E_0_au(this%intensity)
        this%U_p = this%A_0**2/4.0_dp

        this%pol=.false.
        if (pol == 'x') this%pol(1) = .true.
        if (pol == 'y') this%pol(2) = .true.
        if (pol == 'z') this%pol(3) = .true.

        this%env_name = trim(env_name)
        if (env_name == 'cos2') then
            this%env => cos2_env
            this%env_d => cos2_env_d
        else if (env_name == 'gaus') then
            this%env => gaus_env
            this%env_d => gaus_env_d
        else if (env_name == 'gaus_chirp') then
            this%env => gaus_chirp_env
            this%env_d => gaus_chirp_env_d
            if (carrier_name /= "sin_lin_chirp") then
                write(stderr,*) "Error, gaus_chirp envelope must be used with sin_lin_chirp carrier."
                write(stderr,*) "Envelope: ", env_name
                write(stderr,*) "Carrier: ", carrier_name
                stop
            end if
        else if (env_name == 'super_gaus') then
            this%env => super_gaus_env
            this%env_d => super_gaus_env_d
        else if (env_name == 'flat_top') then
            this%env => flat_top_env
            this%env => flat_top_env_d
        else
            write(stderr,*) "Error, unkown envelope type: ", env_name
            stop
        end if

        this%carrier_name = trim(carrier_name)
        if (carrier_name == 'sin') then
            this%carr => sin_carrier
            this%carr_d => sin_carrier_d
        else if(carrier_name == 'sin_lin_chirp') then
            this%carr => sin_lin_chirp_carrier
            this%carr_d => sin_lin_chirp_carrier_d
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

        this%id = id

        call this%print()
    end subroutine init_pulse

    subroutine print_pulse(this)
        class(pulse), intent(in) :: this

        integer :: j

        write(stdout,*)
        write(stdout,'(a,i0)') 'Pulse id: ', this%id
        write(stdout,'(a,a)') 'env_name: ', this%env_name
        write(stdout,'(a)', advance='no') 'env_params: '
        do j = 2,size(this%env_params)
            write(stdout,'(es25.17e3)', advance='no') this%env_params(j)
        end do
        write(stdout,*)
        write(stdout,'(a,a)') 'carrier_name: ', this%carrier_name
        write(stdout,'(a)', advance='no') 'carrier_params: '
        do j = 2,size(this%carrier_params)
            write(stdout,'(es25.17e3)', advance='no') this%carrier_params(j)
        end do
        write(stdout,*)
        write(stdout,'(a,es25.17e3)') 'intensity: ', this%intensity
        write(stdout,'(a,es25.17e3)') 'A_0: ', this%A_0
        write(stdout,'(a,es25.17e3)') 'E_0: ', this%E_0
        write(stdout,'(a,es25.17e3)') 'U_p: ', this%U_p
        write(stdout,'(a,es25.17e3)') 't_0: ', this%env_params(1)
        write(stdout,'(a,l2,l2,l2)') 'pol: ', this%pol
        write(stdout,*)
    end subroutine print_pulse

    function eval_A(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        res = this%A_0*this%carr(t)*this%env(t)
    end function

    function eval_E(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        res = -this%A_0*(this%carr_d(t)*this%env(t) + this%carr(t)*this%env_d(t))
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

        real(dp) :: phase_t

        associate(t_0 => this%carrier_params(1),&
            omega => this%carrier_params(2),&
            phi => this%carrier_params(3),&
            FWHM => this%env_params(2),&
            beta => this%env_params(3)&
        )
        ! omega_t = omega + 4.0_dp*log_2/(beta + 1.0_dp/beta)*(t-t_0)/FWHM**2 ! Instantaneous freq.
        phase_t = omega*(t-t_0) + 2.0_dp*log_2/(beta + 1.0_dp/beta)*(t-t_0)**2/FWHM**2 - phi
        res = sin(phase_t)
        end associate
    end function sin_lin_chirp_carrier

    ! ---------------------------------------------------------------
    ! Carrier derivatives
    ! ---------------------------------------------------------------

    function sin_carrier_d(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%carrier_params(1),&
                omega => this%carrier_params(2),&
                phi => this%carrier_params(3)&
                )

        res = omega*cos(omega*(t-t_0)-phi)
        end associate
    end function sin_carrier_d

    function sin_lin_chirp_carrier_d(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        real(dp) :: omega_t,phase_t

        associate(t_0 => this%carrier_params(1),&
            omega => this%carrier_params(2),&
            phi => this%carrier_params(3),&
            FWHM => this%env_params(2),&
            beta => this%env_params(3)&
        )
        omega_t = omega + 4.0_dp*log_2/(beta + 1.0_dp/beta)*(t-t_0)/FWHM**2 ! Instantaneous freq.
        phase_t = omega*(t-t_0) + 2.0_dp*log_2/(beta + 1.0_dp/beta)*(t-t_0)**2/FWHM**2 - phi
        res = cos(omega_t*(t-t_0) - phi)*(omega_t)
        end associate
    end function sin_lin_chirp_carrier_d

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

    ! ---------------------------------------------------------------
    ! Envelope derivatives
    ! ---------------------------------------------------------------
    function cos2_env_d(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                tau => this%env_params(2)&
                )

        if (abs(t_0-t) < 0.5_dp*tau) then
            res = -2*pi/tau*cos(pi*(t-t_0)/tau)*sin(pi*(t-t_0)/tau)
        else
            res = 0
        end if
        end associate

        ! Danger!!! The derivative is discontinous
    end function cos2_env_d

    function gaus_env_d(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                FWHM => this%env_params(2)&
                )
            res = -log_2*(t-t_0)*((2/FWHM)**2)*exp(-half_log_2*(2*(t-t_0)/FWHM)**2)
        end associate

    end function gaus_env_d

    function gaus_chirp_env_d(this,t) result(res)
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
        res = -amp_beta*log_2*(t-t_0)*(2/FWHM_beta)**2*exp(-half_log_2*(2*(t-t_0)/FWHM_beta)**2)
        end associate
    end function gaus_chirp_env_d

    function super_gaus_env_d(this,t) result(res)
        class(pulse), intent(in) :: this
        real(dp), intent(in) :: t
        real(dp) :: res

        associate(t_0 => this%env_params(1),&
                FWHM => this%env_params(2),&
                n => this%env_params(3)&
                )
            res = -log_2*n*((t-t_0)**(2*n-1))*((2/FWHM)**(2*n))*exp(-half_log_2*(2*(t-t_0)/FWHM)**(2*n))
        end associate
    end function super_gaus_env_d

    function flat_top_env_d(this,t) result(res)
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
        ! This is perhaps not strictly correct (the derivative is not defined at the discontinuities)
    end function flat_top_env_d
end module fields