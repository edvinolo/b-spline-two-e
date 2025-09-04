module wigner_tools
    use kind_tools
    use constants_tools
    use orbital_tools
    implicit none

    interface
        pure function gsl_sf_coupling_3j(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc) result(retval) bind(C)
            use iso_c_binding
            integer(c_int), intent(in), value :: two_ja,two_jb,two_jc,two_ma,two_mb,two_mc
            real(c_double) :: retval
        end function gsl_sf_coupling_3j

        pure function gsl_sf_coupling_6j(two_ja,two_jb,two_jc,two_jd,two_je,two_jf) result(retval) bind(C)
            use iso_c_binding
            integer(c_int), intent(in), value :: two_ja,two_jb,two_jc,two_jd,two_je,two_jf
            real(c_double) :: retval
        end function gsl_sf_coupling_6j

        pure function gsl_sf_legendre_sphPlm(l,m,x) result(retval) bind(C,name="gsl_sf_legendre_sphPlm")
            use iso_c_binding
            integer(c_int), intent(in), value :: l,m
            real(c_double), intent(in), value :: x
            real(c_double) :: retval
        end function gsl_sf_legendre_sphPlm
    end interface
contains

    ! Wrapper for the GSL Wigner-3j function
    pure function three_j(j_a,j_b,j_c,m_a,m_b,m_c) result(retval)
        use iso_c_binding
        integer, intent(in) :: j_a,j_b,j_c,m_a,m_b,m_c
        real(dp) :: retval

        integer(c_int) :: two_ja,two_jb,two_jc,two_ma,two_mb,two_mc
        two_ja = int(2*j_a,kind=c_int)
        two_jb = int(2*j_b,kind=c_int)
        two_jc = int(2*j_c,kind=c_int)
        two_ma = int(2*m_a,kind=c_int)
        two_mb = int(2*m_b,kind=c_int)
        two_mc = int(2*m_c,kind=c_int)

        retval = gsl_sf_coupling_3j(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc)
    end function three_j

    pure function six_j(j_a,j_b,j_c,j_d,j_e,j_f) result(retval)
        use iso_c_binding
        integer, intent(in) :: j_a,j_b,j_c,j_d,j_e,j_f
        real(dp) :: retval

        integer(c_int) :: two_ja,two_jb,two_jc,two_jd,two_je,two_jf
        two_ja = int(2*j_a,kind=c_int)
        two_jb = int(2*j_b,kind=c_int)
        two_jc = int(2*j_c,kind=c_int)
        two_jd = int(2*j_d,kind=c_int)
        two_je = int(2*j_e,kind=c_int)
        two_jf = int(2*j_f,kind=c_int)

        retval = gsl_sf_coupling_6j(two_ja,two_jb,two_jc,two_jd,two_je,two_jf)
    end function six_j

    pure function sph_harm_wrapper(l,m,theta,phi) result(res)
        use iso_c_binding
        integer, intent(in) :: l,m
        real(dp), intent(in) :: theta,phi
        complex(dp) :: res

        integer(c_int) :: l_c,m_c
        real(c_double) :: x

        l_c = int(l,kind=c_int)
        m_c = int(m,kind=c_int)
        x = real(cos(theta),kind=c_double)

        res = exp(i_*m*phi)*gsl_sf_legendre_sphPlm(l_c,m_c,x)
    end function sph_harm_wrapper

    pure function sph_harm(l,m,theta,phi) result(res)
        integer, intent(in) :: l,m
        real(dp), intent(in) :: theta,phi
        complex(dp) :: res

        if ((l < 0).or.(l<abs(m))) then
            res = (0.0_dp,0.0_dp)
        else if (m<0) then
            res = (-1)**m*conjg(sph_harm_wrapper(l,-m,theta,phi))
        else
            res = sph_harm_wrapper(l,m,theta,phi)
        end if
    end function sph_harm

    pure function CG(j_a,j_b,j_c,m_a,m_b,m_c) result(res)
        integer, intent(in) :: j_a,j_b,j_c,m_a,m_b,m_c
        real(dp) :: res

        res = (-1)**(j_b-j_a-m_c)*sqrt(2*j_c + 1.0_dp)*three_j(j_a,j_b,j_c,m_a,m_b,-m_c)
    end function CG

    pure function three_j0(j_a,j_b,j_c) result(res)
        integer, intent(in) :: j_a,j_b,j_c
        real(dp) :: res

        res = three_j(j_a,j_b,j_c,0,0,0)

    end function three_j0

    pure function C_red_mat(k,a,b) result(res)
        integer, intent(in) :: k,a,b
        real(dp) :: res

        res = (-1)**a*sqrt(real((2*a+1)*(2*b+1),kind=dp))*three_j0(a,k,b)
    end function C_red_mat

    pure function X_k(k,orbs) result(res)
        integer, intent(in) :: k
        type(orbital), dimension(4), intent(in) :: orbs
        real(dp) :: res

        if ((mod(orbs(1)%l+k+orbs(3)%l,2) /= 0).or.(mod(orbs(2)%l+k+orbs(4)%l,2) /= 0)) then
            res = 0.d0
        else
            res = (-1)**k*C_red_mat(k,orbs(1)%l,orbs(3)%l)*C_red_mat(k,orbs(2)%l,orbs(4)%l)
        end if
    end function X_k

    pure function ang_k_LS(k,confs,L) result(res)
        integer, intent(in) :: k
        type(config), dimension(2), intent(in) :: confs
        integer, intent(in) :: L
        real(dp) :: res

        if ((mod(confs(1)%l(1)+k+confs(2)%l(1),2) /= 0).or.(mod(confs(1)%l(2)+k+confs(2)%l(2),2) /= 0)) then
            res = 0.d0
        else
            res = (-1)**(confs(1)%l(2)+confs(2)%l(1)+L)*six_j(confs(1)%l(1),confs(1)%l(2),L,confs(2)%l(2),confs(2)%l(1),k) &
                *C_red_mat(k,confs(1)%l(1),confs(2)%l(1))*C_red_mat(k,confs(1)%l(2),confs(2)%l(2))
        end if
    end function ang_k_LS

end module wigner_tools