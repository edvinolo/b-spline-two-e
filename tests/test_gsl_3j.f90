program test_gsl_3j
    use wigner_tools
    use orbital_tools, only: config
    implicit none
    integer :: l_a,l_b,l_c,m_a,m_b,m_c
    integer :: i,j
    !double precision :: three_j
    double precision :: result
    type(config), dimension(2) :: configs

    l_a = 1
    l_b = 1
    l_c = 0
    m_a = 0
    m_b = 0
    m_c = 0

    !three_j = fgsl_sf_coupling_3j(l_a,l_b,l_c,m_a,m_b,m_c)
    result = three_j(l_a,l_b,l_c,m_a,m_b,m_c)

    write(6,*) result

    result = six_j(0,0,0,0,0,0)
    write(6,*) result

    do i=1,8
        do j=1,8
            write(6,*) six_j(i,i,0,j,j,0), i,j
        end do
    end do

    write(6,*)

    do i=1,4
        do j=1,4
            write(6,*) i,j,three_j(i,1,j,-i,1,j-1)
        end do
    end do

    write(6,*)

    configs(1)%l = [0,0]
    configs(2)%l=[2,2]
    result = ang_k_LS(2,configs,0)
    write(6,*) six_j(0,0,0,2,2,2),C_red_mat(2,0,2),result

    write(6,*) sph_harm(2,2,0.1d0,0.2d0),sph_harm(2,-2,0.1d0,0.2d0)
    write(6,*) CG(5,4,1,2,-1,1)

end program test_gsl_3j