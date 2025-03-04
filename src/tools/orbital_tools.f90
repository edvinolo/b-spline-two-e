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

    pure function count_configs(term,max_l,n_b,eigs) result(res)
        type(sym), intent(in) :: term
        integer, intent(in) :: max_l
        integer, intent(in) :: n_b
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
        allocate(temp_list((max_l+1)**2*n_b**2))
        do i = 0,max_l
            orbs(1)%l = i
            orbs(1)%pi = (mod(i,2)/=0)
            do j = 0,i
                orbs(2)%l = j
                orbs(2)%pi = (mod(j,2)/=0)
                do n_i = 1,n_b
                    if (j==i) then
                        do n_j = 1,n_i
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
                        do n_j = 1,n_b
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

    pure function count_terms(max_L) result(res)
        integer, intent(in) :: max_L
        integer :: res

        integer :: i

        res = 1 !L = 0,Pi = -1 is not possible.
        do i = 1,max_L
            res = res + 2*(2*i + 1)
        end do
    end function count_terms

    pure subroutine init_basis(this,max_L,n_b,eigs)
        class(basis), intent(inout) :: this
        integer, intent(in) ::  max_L
        integer, intent(in) :: n_b
        double complex, dimension(:,:), allocatable, intent(in) :: eigs

        integer :: l,m,p,ptr
        logical :: par
        this%n_sym = count_terms(max_L)
        allocate(this%syms(this%n_sym))

        this%syms(1)%l = 0
        this%syms(1)%m = 0
        this%syms(1)%pi = .false.
        this%syms(1)%configs = count_configs(this%syms(1),max_L,n_b,eigs)
        this%syms(1)%n_config = size(this%syms(1)%configs)

        ptr = 2
        do l = 1,max_L
            do p = 0,1
                if (p==0) par = .false.
                if (p==1) par = .true.
                do m = -l,l
                    this%syms(ptr)%l = l
                    this%syms(ptr)%m = m
                    this%syms(ptr)%pi = par
                    !write(6,*) this%syms(ptr)%l,this%syms(ptr)%m,this%syms(ptr)%pi
                    this%syms(ptr)%configs = count_configs(this%syms(ptr),max_L,n_b,eigs)
                    this%syms(ptr)%n_config = size(this%syms(ptr)%configs)
                    ptr = ptr + 1
                end do
            end do
        end do
    end subroutine init_basis
end module orbital_tools