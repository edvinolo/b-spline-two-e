module constants_tools
use kind_tools
implicit none

public
    real(dp), parameter :: length_au = 5.291772109e-11_dp
    real(dp), parameter :: au_to_eV = 27.211386245988_dp
    real(dp), parameter :: time_au = 2.4188843265857e-17_dp
    real(dp), parameter :: fs_to_au = 1e-15_dp/time_au
    real(dp), parameter :: alpha = 0.007297352562787135_dp
    real(dp), parameter :: c_SI = 299792458.0_dp
    real(dp), parameter :: e = 1.60217662e-19_dp
    real(dp), parameter :: hbar = 1.0545718e-34_dp
    real(dp), parameter :: I_au = 3.51e16_dp
    real(dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
    real(dp), parameter :: sqrt_2 = sqrt(2.0_dp)
    real(dp), parameter :: inv_sqrt_2 = 1.0_dp/sqrt_2
    real(dp), parameter :: log_2 = log(2.0_dp)
    real(dp), parameter :: half_log_2 = 0.5_dp*log(2.0_dp)
    complex(dp), parameter :: i_ = (0.0_dp,1.0_dp)

contains

    pure function I_W_cm2(E_0) result(res)
        real(dp), intent(in) :: E_0
        real(dp) :: res

        res = E_0**2*I_au
    end function I_W_cm2

    pure function E_0_au(I) result(res)
        real(dp), intent(in) :: I
        real(dp) :: res

        res = sqrt(I/I_au)
    end function E_0_au

    pure function A_0_au(I,omega) result(res)
        real(dp), intent(in) :: I
        real(dp), intent(in) :: omega
        real(dp) :: res

        res = sqrt(I/I_au)/omega
    end function A_0_au

end module constants_tools