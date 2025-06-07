module orbital_tools
    implicit none

    type, public :: orbital
        integer :: l
        integer :: m
        logical :: pi !parity of orbital, even = .false., odd = .true.
    end type orbital

    type, extends(orbital), public :: sym
        integer :: n_config
        type(config), dimension(:), allocatable :: configs
    end type

    type, public :: config
        integer, dimension(2) :: n
        integer, dimension(2) :: l
        logical :: eqv
    end type config

    type, public :: basis
        integer :: n_sym
        type(sym), dimension(:), allocatable :: syms
    contains
        procedure :: init => init_basis
        procedure :: store => store_basis
        procedure :: load => load_basis
    end type basis
contains

    pure function pi(orbs) result(res)
        type(orbital), dimension(2),  intent(in) :: orbs
        logical :: res

        res = parity([orbs(1)%pi,orbs(2)%pi])
    end function pi

    pure function consistent(orbs,term,eqv) result(res)
        type(orbital), dimension(2), intent(in) :: orbs
        type(sym), intent(in) :: term
        logical, intent(in) :: eqv
        logical :: res

        logical :: tri, par

        if (eqv.and.(mod(term%l,2)/=0)) then
            res = .false.
            return
        end if

        tri = (abs(orbs(1)%l-orbs(2)%l) <= term%l).and.(term%l<= orbs(1)%l + orbs(2)%l)
        if (.not.tri) then
            res = .false.
            return
        end if

        par = (pi(orbs) .eqv. term%pi)
        if(.not.par) then
            res = .false.
            return
        end if

        res = .true.
    end function consistent

    pure function dip_allowed(term_1,term_2,pol) result(res)
        type(sym), intent(in) :: term_1
        type(sym), intent(in) :: term_2
        character, intent(in) :: pol
        logical :: res

        ! Parity must change
        if (term_1%pi.neqv.term_2%pi) then
            res = .false.
            return
        end if

        ! L = 0 <--> 0 is not allowed
        if((term_1%l==0).and.(term_2%l==0)) then
            res = .false.
            return
        end if

        ! Delta_L = 0,\pm 1
        if (((term_1%l-term_2%l)>1).or.(term_1%l-term_2%l<-1)) then
            res = .false.
            return
        end if

        if ((pol=='x').or.(pol=='y')) then
            ! Delta_M_L = \pm 1
            if (((term_1%m-term_2%m)>1).or.(term_1%m-term_2%m<-1)) then
                res = .false.
                return
            end if
        else if (pol=='z') then
            ! Delta_M_L = 0
            if (term_1%m/=term_2%m) then
                res = .false.
                return
            ! Delta_M_L != 0 if Delta_L = 0
            else if ((term_1%l==term_2%l)) then
                res = .false.
                return
            end if
        end if

        res = .true.
    end function dip_allowed

    pure function count_configs(term,max_l_1p,n_b,k_spline,max_n_b,n_all_l,l_2_max,eigs) result(res)
        type(sym), intent(in) :: term
        integer, intent(in) :: max_l_1p
        integer, intent(in) :: n_b
        integer, intent(in) :: k_spline
        integer, intent(in) :: max_n_b
        integer, intent(in) :: n_all_l
        integer, intent(in) :: l_2_max
        double complex, dimension(:,:), allocatable, intent(in) :: eigs
        type(config), dimension(:), allocatable :: res

        integer :: i,j,n_i,n_j,ptr
        logical :: eqv
        double complex :: energy_1p
        type(orbital), dimension(2) :: orbs
        type(config) :: conf
        type(config), dimension(:), allocatable :: temp_list

        orbs%m = 0
        ptr = 1
        energy_1p = 0.d0
        allocate(temp_list((max_l_1p+1)**2*n_b**2))
        do i = 0,max_l_1p
            orbs(1)%l = i
            orbs(1)%pi = (mod(i,2)/=0)
            do j = 0,i
                ! if (j+term%l > max_l_1p+1) cycle !Limit the angular momentum of second electron
                orbs(2)%l = j
                orbs(2)%pi = (mod(j,2)/=0)
                do n_i = min(i+1,k_spline-1),n_b
                    if ((n_i>n_all_l).and.(j>l_2_max)) cycle !limit angular momentum of second electron when the other is far away
                    if (j==i) then
                        do n_j = min(j+1,k_spline-1),min(n_i,max_n_b)
                            eqv = (i==j).and.(n_i==n_j)
                            energy_1p = eigs(n_i,i) + eigs(n_j,j)
                            if (consistent(orbs,term,eqv)) then!.and.(real(energy_1p)<1.5d1)) then
                                conf%n = [n_i,n_j]
                                conf%l = [i,j]
                                conf%eqv = eqv
                                temp_list(ptr) = conf
                                ptr = ptr + 1
                            end if
                        end do
                    else
                        do n_j = min(j+1,k_spline-1),min(n_b,max_n_b)
                            eqv = (i==j).and.(n_i==n_j)
                            energy_1p = eigs(n_i,i) + eigs(n_j,j)
                            if (consistent(orbs,term,eqv)) then!.and.(real(energy_1p)<1.5d1)) then
                                conf%n = [n_i,n_j]
                                conf%l = [i,j]
                                conf%eqv = eqv
                                temp_list(ptr) = conf
                                ptr = ptr + 1
                            end if
                        end do
                    end if
                end do
            end do
        end do

        allocate(res(ptr-1))
        res = temp_list(1:ptr-1)
    end function count_configs

    pure function count_terms(max_L,z_pol) result(res)
        integer, intent(in) :: max_L
        logical, intent(in) :: z_pol
        integer :: res

        !Include ony those that are reachable from S^e block.
        if (z_pol) then
            res = max_L + 1
        else
            res = (max_L +1)**2
        end if
    end function count_terms

    pure subroutine init_basis(this,max_L,max_l_1p,n_b,k_spline,max_n_b,n_all_l,l_2_max,z_pol,eigs)
        class(basis), intent(inout) :: this
        integer, intent(in) ::  max_L
        integer, intent(in) :: max_l_1p
        integer, intent(in) :: n_b
        integer, intent(in) :: k_spline
        integer, intent(in) :: max_n_b
        integer, intent(in) :: n_all_l
        integer, intent(in) :: l_2_max
        logical, intent(in) :: z_pol
        double complex, dimension(:,:), allocatable, intent(in) :: eigs

        integer :: l,m,p,ptr
        logical :: par
        this%n_sym = count_terms(max_L,z_pol)
        allocate(this%syms(this%n_sym))

        this%syms(1)%l = 0
        this%syms(1)%m = 0
        this%syms(1)%pi = .false.
        this%syms(1)%configs = count_configs(this%syms(1),max_l_1p,n_b,k_spline,max_n_b,n_all_l,l_2_max,eigs)
        this%syms(1)%n_config = size(this%syms(1)%configs)

        ptr = 2

        if (z_pol) then
            do l = 1,max_L
                this%syms(ptr)%l = l
                this%syms(ptr)%m = 0
                this%syms(ptr)%pi = (mod(l,2)/=0)
                !write(6,*) this%syms(ptr)%l,this%syms(ptr)%m,this%syms(ptr)%pi
                this%syms(ptr)%configs = count_configs(this%syms(ptr),max_l_1p,n_b,k_spline,&
                                                        max_n_b,n_all_l,l_2_max,eigs)
                this%syms(ptr)%n_config = size(this%syms(ptr)%configs)
                ptr = ptr + 1
            end do
        else
            do l = 1,max_L
                do p = 0,1
                    if (p==0) par = .false.
                    if (p==1) par = .true.
                    do m = -l,l
                        if (abs(mod(m,2))/=p) cycle
                        this%syms(ptr)%l = l
                        this%syms(ptr)%m = m
                        this%syms(ptr)%pi = par
                        !write(6,*) this%syms(ptr)%l,this%syms(ptr)%m,this%syms(ptr)%pi
                        this%syms(ptr)%configs = count_configs(this%syms(ptr),max_l_1p,n_b,k_spline,&
                                                                max_n_b,n_all_l,l_2_max,eigs)
                        this%syms(ptr)%n_config = size(this%syms(ptr)%configs)
                        ptr = ptr + 1
                    end do
                end do
            end do
        end if
    end subroutine init_basis

    subroutine store_basis(this,loc)
        class(basis), intent(in) :: this
        character(len=*), intent(in) ::  loc

        integer :: i, unit
        character(len=:), allocatable :: format_header,format_data

        open(file = loc//"basis.dat", newunit = unit, action = 'write', form = "unformatted")
        write(unit) this%n_sym
        do i = 1, this%n_sym
            write(unit) this%syms(i)%l
            write(unit) this%syms(i)%m
            write(unit) this%syms(i)%pi
            write(unit) this%syms(i)%n_config
            write(unit) this%syms(i)%configs
        end do
        close(unit)

        format_header = "(a5,a5,a5,a5,a9)"
        format_data = "(i5,i5,i5,l5,i9)"
        open(file = loc//"basis_info.txt", newunit = unit, action = 'write')
        write(unit,format_header) "# Sym", "L", "M", "Pi", "states"
        do i = 1,this%n_sym
            write(unit,format_data) i, this%syms(i)%l, this%syms(i)%m, this%syms(i)%pi, this%syms(i)%n_config
        end do
        close(unit)
    end subroutine store_basis

    subroutine load_basis(this,loc)
        class(basis), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        integer :: i, unit

        open(file = loc//"basis.dat", newunit = unit, action = 'read', form = "unformatted")
        read(unit) this%n_sym
        allocate(this%syms(this%n_sym))
        do i = 1, this%n_sym
            read(unit) this%syms(i)%l
            read(unit) this%syms(i)%m
            read(unit) this%syms(i)%pi
            read(unit) this%syms(i)%n_config
            allocate(this%syms(i)%configs(this%syms(i)%n_config))
            read(unit) this%syms(i)%configs
        end do
        close(unit)
    end subroutine load_basis
end module orbital_tools