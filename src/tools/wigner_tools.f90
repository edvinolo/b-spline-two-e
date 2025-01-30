module wigner_tools
    implicit none

    interface
        function gsl_sf_coupling_3j(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc) result(retval) bind(C)
            use iso_c_binding
            integer(c_int), intent(in),value :: two_ja,two_jb,two_jc,two_ma,two_mb,two_mc
            real(c_double) :: retval
        end function gsl_sf_coupling_3j
    end interface
contains

    ! Wrapper for the GSL Wigner-3j function
    function three_j(j_a,j_b,j_c,m_a,m_b,m_c) result(retval)
        use iso_c_binding
        integer, intent(in) :: j_a,j_b,j_c,m_a,m_b,m_c
        double precision :: retval

        integer(c_int) :: two_ja,two_jb,two_jc,two_ma,two_mb,two_mc
        two_ja = 2*j_a
        two_jb = 2*j_b
        two_jc = 2*j_c
        two_ma = 2*m_a
        two_mb = 2*m_b
        two_mc = 2*m_c

        retval = gsl_sf_coupling_3j(two_ja,two_jb,two_jc,two_ma,two_mb,two_mc)
    end function three_j

end module wigner_tools