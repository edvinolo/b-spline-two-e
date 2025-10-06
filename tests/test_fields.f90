program test_fields
    use fields
    use stdlib_math, only: linspace
    implicit none

    real(dp), allocatable :: t(:)
    real(dp), allocatable :: A_t(:)
    real(dp), allocatable :: env_params(:), carrier_params(:)
    type(pulse) :: pulse_pair(2)

    integer :: i,j,unit

    character(len=30) :: env_name(2)
    character(len=30) :: carrier_name(2)

    real(dp), parameter :: t_0(2) = [0,1000_dp]
    real(dp), parameter :: omega(2) = [1.0_dp,0.057_dp]
    real(dp), parameter :: tau(2) = [50.0_dp,300.0_dp]
    real(dp), parameter :: A_0(2) = [0.01_dp,0.02_dp]
    real(dp), parameter :: phi(2) = [0.0_dp,0.5_dp*pi]
    logical, parameter :: pol(3) = [.false.,.false.,.true.]
    real(dp), parameter :: beta = 1.0_dp

    integer(int32), parameter :: N_t = 10000
    real(dp), parameter :: t_1 = -1.0e2_dp
    real(dp), parameter :: t_2 = 2.0e3_dp

    env_name(1) = 'cos2'
    env_name(2) = 'gaus_chirp'

    carrier_name(1) = 'sin'
    carrier_name(2) = 'sin_lin_chirp'

    do i = 1,2
        if (i == 1) then 
            env_params = [tau(i)]
        else 
            env_params = [tau(i),beta]
        end if
        carrier_params = [omega(i),phi(i)]
        call pulse_pair(i)%init(A_0(i),t_0(i),env_params,carrier_params,pol,env_name(i), carrier_name(i))
    end do

    t = linspace(t_1,t_2,N_t)
    allocate(A_t(N_t), source = 0.0_dp)
    do i = 1,N_t
        do j = 1,2
            A_t(i) = A_t(i) + pulse_pair(j)%eval_A(t(i))
        end do
    end do

    open(file = 'test_fields.dat',newunit = unit, action = 'write')
    do i = 1,N_t
        write(unit,'(es25.17e3,a,es25.17e3)') t(i),' ', A_t(i)
    end do
    close(unit)

end program test_fields