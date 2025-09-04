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
        integer :: max_l_1p
        integer :: max_L
        logical :: two_el
        integer :: n_sym
        integer :: n_states
        type(sym), allocatable :: syms(:)
        integer, allocatable :: sym_ptr(:)
    contains
        procedure :: init => init_basis
        procedure :: compute_size => compute_basis_size
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

        integer :: i,j,n_i,n_j,ptr,n_orbitals
        logical :: eqv
        double complex :: energy_1p
        type(orbital), dimension(2) :: orbs
        type(config) :: conf
        type(config), dimension(:), allocatable :: temp_list
        integer, dimension(:,:), allocatable :: orbital_list,temp_orbital

        orbs%m = 0
        ptr = 1
        energy_1p = 0.d0
        n_orbitals = (max_l_1p+1)*n_b
        allocate(temp_list(n_orbitals**2),temp_orbital(2,n_orbitals))

        do i = 0,max_l_1p
            do n_i = min(i+1,k_spline-1),n_b
                temp_orbital(1,ptr) = n_i
                temp_orbital(2,ptr) = i
                ptr = ptr + 1
            end do
        end do

        n_orbitals = ptr-1
        allocate(orbital_list(2,n_orbitals))
        orbital_list = temp_orbital(:,:n_orbitals)

        ptr = 1
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

        ! ptr = 1
        ! do i = 1,n_orbitals
        !     orbs(1)%l = orbital_list(2,i)
        !     orbs(1)%pi = (mod(orbs(1)%l,2)/=0)
        !     do j = i,n_orbitals
        !         ! if (j+term%l > max_l_1p+1) cycle !Limit the angular momentum of second electron
        !         orbs(2)%l = orbital_list(2,j)
        !         orbs(2)%pi = (mod(orbs(2)%l,2)/=0)
        !         eqv = (i==j)
        !         if (consistent(orbs,term,eqv)) then
        !             conf%n = [orbital_list(1,i),orbital_list(1,j)]
        !             conf%l = [orbs(1)%l,orbs(2)%l]
        !             conf%eqv = eqv
        !             temp_list(ptr) = conf
        !             ptr = ptr + 1
        !         end if
        !     end do
        ! end do

        allocate(res(ptr-1))
        res = temp_list(1:ptr-1)
    end function count_configs

    pure function count_configs_1p(term,n_b,k_spline) result(res)
        type(sym), intent(in) :: term
        integer, intent(in) :: n_b
        integer, intent(in) :: k_spline
        type(config), dimension(:), allocatable :: res

        integer :: n_i,ptr,n_orbitals
        type(config) :: conf
        type(config), dimension(:), allocatable :: temp_list

        ptr = 1
        n_orbitals = n_b
        allocate(temp_list(n_orbitals))

        ptr = 1
        do n_i = min(term%l+1,k_spline-1),n_b
            conf%n = [n_i,-1]
            conf%l = [term%l,-1]
            conf%eqv = .false.
            temp_list(ptr) = conf
            ptr = ptr + 1
        end do

        allocate(res(ptr-1))
        res = temp_list(1:ptr-1)
    end function count_configs_1p

    pure function count_terms(max_L,z_pol,two_el) result(res)
        integer, intent(in) :: max_L
        logical, intent(in) :: z_pol
        logical, intent(in) :: two_el
        integer :: res

        !Include ony those that are reachable from S^e block.
        if (z_pol) then
            res = max_L + 1
        else
            if (two_el) then
                res = (max_L +1)**2
            else
                res = max_L + 1 + max_L*(max_L + 1)/2
            end if
        end if
    end function count_terms

    pure subroutine init_basis(this,max_L,max_l_1p,n_b,k_spline,max_n_b,n_all_l,l_2_max,z_pol,two_el,eigs)
        class(basis), intent(inout) :: this
        integer, intent(in) ::  max_L
        integer, intent(in) :: max_l_1p
        integer, intent(in) :: n_b
        integer, intent(in) :: k_spline
        integer, intent(in) :: max_n_b
        integer, intent(in) :: n_all_l
        integer, intent(in) :: l_2_max
        logical, intent(in) :: z_pol
        logical, intent(in) :: two_el
        double complex, dimension(:,:), allocatable, intent(in) :: eigs

        integer :: l,m,p,ptr
        logical :: par

        this%max_l_1p = max_l_1p
        this%max_L = max_L
        this%two_el = two_el
        this%n_sym = count_terms(max_L,z_pol,two_el)
        allocate(this%syms(this%n_sym))

        this%syms(1)%l = 0
        this%syms(1)%m = 0
        this%syms(1)%pi = .false.
        if (two_el) then
            this%syms(1)%configs = count_configs(this%syms(1),max_l_1p,n_b,k_spline,max_n_b,n_all_l,l_2_max,eigs)
        else
            this%syms(1)%configs = count_configs_1p(this%syms(1),n_b,k_spline)
        end if
        this%syms(1)%n_config = size(this%syms(1)%configs)

        ptr = 2

        if (z_pol) then
            do l = 1,max_L
                this%syms(ptr)%l = l
                this%syms(ptr)%m = 0
                this%syms(ptr)%pi = (mod(l,2)/=0)
                !write(6,*) this%syms(ptr)%l,this%syms(ptr)%m,this%syms(ptr)%pi
                if (two_el) then
                    this%syms(ptr)%configs = count_configs(this%syms(ptr),max_l_1p,n_b,k_spline,&
                                                        max_n_b,n_all_l,l_2_max,eigs)
                else
                    this%syms(ptr)%configs = count_configs_1p(this%syms(ptr),n_b,k_spline)
                end if
                this%syms(ptr)%n_config = size(this%syms(ptr)%configs)
                ptr = ptr + 1
            end do
        else
            do l = 1,max_L
                if (two_el) then
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
                else
                    par = (mod(l,2)/=0)
                    do m = -l,l,2
                        this%syms(ptr)%l = l
                        this%syms(ptr)%m = m
                        this%syms(ptr)%pi = par
                        this%syms(ptr)%configs = count_configs_1p(this%syms(ptr),n_b,k_spline)
                        this%syms(ptr)%n_config = size(this%syms(ptr)%configs)
                        ptr = ptr + 1
                    end do
                end if
            end do
        end if
    end subroutine init_basis

    subroutine compute_basis_size(this)
        class(basis), intent(inout) :: this

        integer :: i,ptr

        this%n_states = sum(this%syms%n_config)

        if (allocated(this%sym_ptr)) deallocate(this%sym_ptr)
        allocate(this%sym_ptr(this%n_sym+1))

        ptr = 1
        do i = 1,this%n_sym
            this%sym_ptr(i) = ptr
            ptr = ptr + this%syms(i)%n_config
        end do
        this%sym_ptr(this%n_sym+1) = ptr

    end subroutine compute_basis_size

    subroutine store_basis(this,loc)
        class(basis), intent(inout) :: this
        character(len=*), intent(in) ::  loc

        integer :: i, unit
        character(len=:), allocatable :: format_header,format_data

        call this%compute_size()

        open(file = loc//"basis.dat", newunit = unit, action = 'write', form = "unformatted")
        write(unit) this%max_l_1p
        write(unit) this%max_L
        write(unit) this%two_el
        write(unit) this%n_sym
        write(unit) this%n_states
        write(unit) this%sym_ptr
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
        read(unit) this%max_l_1p
        read(unit) this%max_L
        read(unit) this%two_el
        read(unit) this%n_sym
        read(unit) this%n_states

        allocate(this%sym_ptr(this%n_sym+1))
        read(unit) this%sym_ptr

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