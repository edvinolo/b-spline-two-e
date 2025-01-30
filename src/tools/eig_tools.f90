module eig_tools
    use stdlib_linalg_lapack, only: ggev
    use stdlib_linalg_blas, only: dotc,dotu,gemv
    use iso_fortran_env, only: stderr => error_unit
    implicit none

contains
    ! Normalize vecs according to v.T*B*v = id
    ! Requires copying B, since it is destroyed by ggev
    subroutine eig_general(A,B,eigs,vecs)
        double complex, dimension(:,:), intent(inout) :: A
        double complex, dimension(:,:), intent(inout) :: B
        double complex, dimension(:), intent(out) :: eigs
        double complex, dimension(:,:), intent(out) :: vecs

        integer :: N
        character :: jobvl,jobvr
        double complex, dimension(:), allocatable :: alpha,beta
        double complex, dimension(1,1) :: vl
        double complex, dimension(:), allocatable :: work
        integer :: lwork
        double precision, dimension(:), allocatable :: rwork
        integer :: info
        double complex, dimension(:,:), allocatable :: B_copy
        integer :: i
        double complex, dimension(:), allocatable :: temp_v

        N = size(eigs)
        jobvl = 'N'
        jobvr = 'V'
        allocate(alpha(N),beta(N),work(1))
        lwork = -1 !First ggev call is workspace query
        allocate(rwork(8*N))
        info = 0
        allocate(B_copy(N,N))
        B_copy = B

        call ggev(jobvl,jobvr,N,A,N,B,N,alpha,beta,vl,N,vecs,N,work,lwork,rwork,info)

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call ggev(jobvl,jobvr,N,A,N,B,N,alpha,beta,vl,N,vecs,N,work,lwork,rwork,info)

        if (info /= 0) then
            write(stderr,*) "zggev exited with info: ", info
            write(stderr,*) "Please consult lapack documentation for explanation."
            stop
        end if

        eigs = alpha/beta !Should check for zeros in beta if properly doing this

        allocate(temp_v(N))
        ! Normalize the eigenvectors
        do i=1,N
            call gemv('N',N,N,dcmplx(1.d0,0.d0),B_copy,N,vecs(:,i),1,dcmplx(0.d0,0.d0),temp_v,1)
            vecs(:,i) = vecs(:,i)/sqrt(dotu(N,vecs(:,i),1,temp_v,1))
        end do

    end subroutine eig_general
end module eig_tools