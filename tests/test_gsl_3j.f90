program test_gsl_3j
    !use fgsl
    use wigner_tools
    implicit none
    integer :: l_a,l_b,l_c,m_a,m_b,m_c
    !double precision :: three_j
    double precision :: result

    l_a = 1
    l_b = 1
    l_c = 0
    m_a = 0
    m_b = 0
    m_c = 0

    !three_j = fgsl_sf_coupling_3j(l_a,l_b,l_c,m_a,m_b,m_c)
    result = three_j(l_a,l_b,l_c,m_a,m_b,m_c)

    write(6,*) result

end program test_gsl_3j