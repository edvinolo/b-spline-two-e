module propagator_CN
    use kind_tools
    use constants_tools
    use sparse_array_tools, only: CSR_matrix,CSR_mv,CSR_mv_sym,mv
    use block_tools, only: block_CS, block_diag_CS, XPAY, AXPBY
    use GMRES_tools, only: zFGMRES, mv_GMRES
    use precond_tools, only: block_Jacobi_PC
    use fields, only: pulse
    implicit none

    logical, private :: pol_active(3)
    complex(dp), private :: A_t(3), E_t(3)
    type(pulse), allocatable, private :: pulses(:)
    integer, private :: n_pulses

    type(CSR_matrix), private :: A
    type(CSR_matrix), private :: B
    type(block_diag_CS), private :: A_block
    type(block_diag_CS), private :: B_block
    type(CSR_matrix), private :: dipole(3)
    type(block_CS), private :: D
    logical, private :: full
    integer, private :: n

    complex(dp), allocatable, private :: temp_vec(:,:)
    complex(dp), allocatable, private :: temp_x(:),temp_y(:),temp(:)
    procedure(mv), pointer, private :: matvec => null()

    type(zFGMRES), private :: solver
    type(block_Jacobi_PC), private :: precond

    complex(dp), private :: forward_factor,backward_factor

contains
    ! Steps needed to implment Crank-Nicholson propagator (S+0.5*1j*H*dt)^(-1)(S-0.5*1j*H*dt)
    !
    ! The right parenthesis can be implemented as several independent matrix vector products and vector additions,
    ! So there is no need to explicitly form H(t). For the left one, it might only be necessary to form the matrix if a direct solver is used,
    ! provided that the "diagonal" preconditioner is sufficiently good (which is hopefully true). In that case, need to ensure that the iterative solver
    ! has access to the correct matvec.

    subroutine init_CN(pulses_in,pol_active_in,dt_in,H_0_block,S_block,dip,full_in,B_subset)
        type(pulse), intent(in) :: pulses_in(:)
        logical, intent(in) :: pol_active_in(3)
        real(dp), intent(in) :: dt_in
        type(block_diag_CS), intent(inout) :: H_0_block
        type(block_diag_CS), intent(inout) :: S_block
        type(block_CS), intent(inout) :: dip(-1:1)
        logical, intent(in) :: full_in
        logical, intent(in) :: B_subset

        n_pulses = size(pulses_in)
        pulses = pulses_in
        pol_active = pol_active_in

        call set_factors(dt_in)

        n = H_0_block%shape(1)
        allocate(temp_vec(n,3),temp_x(n),temp_y(n),temp(n))
        call setup_A_and_B(H_0_block,S_block,B_subset)

        call setup_dipole(dip)

        full = full_in
        call set_matvec()

        call setup_solver()
    end subroutine init_CN

    subroutine eval_A_t(t)
        real(dp), intent(in) :: t

        integer :: i,j

        A_t = 0
        do i = 1,3
            if (pol_active(i)) then
                do j = 1,n_pulses
                    if (pulses(j)%pol(i)) then
                        A_t(i) = A_t(i) + pulses(j)%eval_A(t)
                    end if
                end do
            end if
        end do
    end subroutine eval_A_t

    subroutine eval_E_t(t)
        real(dp), intent(in) :: t

        integer :: i,j

        E_t = 0
        do i = 1,3
            if (pol_active(i)) then
                do j = 1,n_pulses
                    if (pulses(j)%pol(i)) then
                        E_t(i) = E_t(i) + pulses(j)%eval_E(t)
                    end if
                end do
            end if
        end do
    end subroutine eval_E_t

    function get_A_t(t) result(res)
        real(dp), intent(in) :: t
        real(dp) :: res(3)

        call eval_A_t(t)
        res = real(A_t,kind=dp)
    end function get_A_t

    function get_E_t(t) result(res)
        real(dp), intent(in) :: t
        real(dp) :: res(3)

        call eval_E_t(t)
        res = real(E_t,kind=dp)
    end function get_E_t

    subroutine CN_step(t,x,y,A_ret,E_ret)
        real(dp), intent(in) :: t
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        real(dp), intent(out) :: A_ret(3)
        real(dp), intent(out) :: E_ret(3)

        call eval_A_t(t)
        call eval_E_t(t)
        call CN_forward(x,temp)
        call CN_backward(temp,y)
        A_ret = real(A_t,kind=dp)
        E_ret = real(E_t,kind=dp)
    end subroutine CN_step

    subroutine CN_forward(x,y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) ::  y(:)

        integer :: i
        complex(dp) :: factor

        ! A = S - 0.5*1j*dt*H_0 is constant, so can be computed before the propagation and used in every time step
        ! Matvec for y = A*x
        ! GEMV for y = y + V_i(t)x, for each active polarization

        call matvec(A,x,y)

        do i = 1,3
            if (pol_active(i)) then
                call matvec(dipole(i),x,temp_vec(:,i))
                factor = A_t(i)*forward_factor
                call zaxpy(size(x),factor,temp_vec(:,i),1,y,1)
            end if
        end do

    end subroutine CN_forward

    subroutine CN_backward(x,y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        ! Here, an iterative solver is assumed, with a fixed preconditioner
        ! B = S + 0.5*1j*dt*H_0 is constant, can be computed before the propagation and used in every time step
        ! The preconditioner could be (an approximation to) S^(-1) or B^(-1)
        ! Since S is real, one could perhaps reduce the memory footprint by utilizing that, but maybe unneccessary
        ! Matvec that should be used by the iterative solver would then be
        ! Matvec for y = B*x
        ! GEMV for y = y + V_i(t)x, for each active polarization

        call solver%solve(y,x,precond)
    end subroutine CN_backward

    subroutine CN_backward_matvec(this,x,y)
        class(zFGMRES), intent(in) :: this
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)

        integer :: i
        complex(dp) :: factor

        do i = 1,this%n
            temp_x(i) = cmplx(x(i),x(i+this%n),kind=dp)
        end do

        call matvec(B,temp_x,temp_y)

        do i = 1,3
            if (pol_active(i)) then
                call matvec(dipole(i),temp_x,temp_vec(:,i))
                factor = A_t(i)*backward_factor
                call zaxpy(this%n,factor,temp_vec(:,i),1,temp_y,1)
            end if
        end do

        do i = 1,this%n
            y(i) = real(temp_y(i),kind=dp)
            y(this%n + i) = aimag(temp_y(i))
        end do
    end subroutine CN_backward_matvec

    subroutine set_matvec()
        if (full) then
            matvec => CSR_mv
        else
            matvec => CSR_mv_sym
        end if
    end subroutine set_matvec

    subroutine set_factors(dt)
        real(dp), intent(in) :: dt
        forward_factor = -0.5_dp*i_*dt
        backward_factor = -forward_factor
    end subroutine set_factors

    subroutine setup_A_and_B(H_0_block,S_Block,B_subset)
        type(block_diag_CS), intent(inout) :: H_0_block
        type(block_diag_CS), intent(inout) :: S_block
        logical, intent(in) :: B_subset

        A_block = forward_factor*H_0_block
        B_block = backward_factor*H_0_block

        call H_0_block%deall()

        call A_block%shift_B((1.0_dp,0.0_dp),S_block,B_subset)
        call B_block%shift_B((1.0_dp,0.0_dp),S_block,B_subset)

        ! call S_block%deall()

        call A_block%to_CS(A,.true.)
        call B_block%to_CS(B,.false.)
    end subroutine setup_A_and_B

    subroutine setup_dipole(dip)
        type(block_CS), intent(inout) :: dip(-1:1)

        integer :: i

        if (pol_active(1)) then
            ! X = -/sqrt(2)*(T_1 - T_-1)
            D = AXPBY(dip(1), -(inv_sqrt_2,0.0_dp), dip(-1), (inv_sqrt_2,0.0_dp))
            call D%to_CS(dipole(1),.false.)
        end if

        if (pol_active(2)) then
            ! Y = i/sqrt(2)*(T_1 + T_-1)
            D = AXPBY(dip(1), i_*inv_sqrt_2, dip(-1), i_*inv_sqrt_2)
            call D%to_CS(dipole(2),.false.)
        end if

        if (pol_active(3)) then
            ! Z = T_0
            call dip(0)%to_CS(dipole(3),.false.)
        end if

        do i = -1,1
            call dip(i)%deall()
        end do
    end subroutine setup_dipole

    subroutine setup_solver()
        call precond%setup_diag(B_block)
        call B_block%deall()
        call solver%setup(B,full,store_matrix=.false.)
        solver%mv => CN_backward_matvec
    end subroutine setup_solver

    subroutine cleanup_CN()
        call A_block%deall()
        call B_block%deall()

        call A%deall()
        call B%deall()

        call precond%cleanup()
        call solver%cleanup()

        deallocate(temp_vec,temp_x,temp_y,temp)
    end subroutine cleanup_CN

end module propagator_CN