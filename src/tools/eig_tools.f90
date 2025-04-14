module eig_tools
    use sparse_array_tools, only: CSR_matrix,CSC_matrix
    use stdlib_linalg_lapack, only: ggev
    use stdlib_linalg_blas, only: dotu,gemv
    use stdlib_sorting, only: sort_index
    use iso_fortran_env, only: stderr => error_unit
    use omp_lib, only: omp_get_wtime
    implicit none

contains
    ! Normalize vecs according to v.T*B*v = id
    ! Requires copying B, since it is destroyed by ggev
    subroutine eig_general(A,B,eigs,vecs)
        double complex, dimension(:,:), intent(inout) :: A
        double complex, dimension(:,:), intent(in) :: B
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

        call zggev3(jobvl,jobvr,N,A,N,B_copy,N,alpha,beta,vl,N,vecs,N,work,lwork,rwork,info)

        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        call zggev3(jobvl,jobvr,N,A,N,B_copy,N,alpha,beta,vl,N,vecs,N,work,lwork,rwork,info)

        if (info /= 0) then
            write(stderr,*) "zggev exited with info: ", info
            write(stderr,*) "Please consult lapack documentation for explanation."
            stop
        end if

        eigs = alpha/beta !Should check for zeros in beta if properly doing this

        allocate(temp_v(N))
        ! Normalize the eigenvectors
        do i=1,N
            call gemv('N',N,N,dcmplx(1.d0,0.d0),B,N,vecs(:,i),1,dcmplx(0.d0,0.d0),temp_v,1)
            vecs(:,i) = vecs(:,i)/sqrt(dotu(N,vecs(:,i),1,temp_v,1))
        end do

    end subroutine eig_general

    subroutine B_normalize(B,vecs)
        double complex, dimension(:,:), intent(in) :: B
        double complex, dimension(:,:), intent(out) ::  vecs

        integer :: i,N
        integer, dimension(2) :: N_2
        double complex, dimension(:), allocatable :: temp_V

        N_2 = shape(B)
        N = N_2(1)
        allocate(temp_V(N))
        ! Normalize the eigenvectors
        do i=1,N
            call gemv('N',N,N,dcmplx(1.d0,0.d0),B,N,vecs(:,i),1,dcmplx(0.d0,0.d0),temp_v,1)
            vecs(:,i) = vecs(:,i)/sqrt(dotu(N,vecs(:,i),1,temp_v,1))
        end do

    end subroutine B_normalize

    ! Sort eigenvalues and eigenvector based on the real part of the eigenvalues
    subroutine sort_eigvals(n_b,eigs)
        integer, intent(in) :: n_b
        double complex, dimension(n_b), intent(inout) :: eigs

        integer, dimension(n_b) :: index
        double precision, dimension(n_b) :: eigs_real

        eigs_real = real(eigs)

        call sort_index(eigs_real,index)

        eigs = eigs(index)
    end subroutine sort_eigvals

    ! Sort eigenvalues and eigenvector based on the real part of the eigenvalues
    subroutine sort_eig(n_b,eigs,vecs)
        integer, intent(in) :: n_b
        double complex, dimension(n_b), intent(inout) :: eigs
        double complex, dimension(n_b,n_b), intent(inout) ::  vecs

        integer, dimension(n_b) :: index
        double precision, dimension(n_b) :: eigs_real

        eigs_real = real(eigs)

        call sort_index(eigs_real,index)

        eigs = eigs(index)
        vecs = vecs(:,index)
    end subroutine sort_eig

    subroutine FEAST(H,S,A,B,mid,rad,eigs,vecs,M)
        double complex, dimension(:,:), intent(in) :: H
        double complex, dimension(:,:), intent(in) :: S
        type(CSR_matrix), intent(inout) :: A
        type(CSR_matrix), intent(inout) :: B
        double complex, intent(in) :: mid
        double precision, intent(in) :: rad
        double complex, dimension(:), allocatable, intent(inout) :: eigs
        double complex, dimension(:,:), allocatable, intent(inout) :: vecs
        integer, intent(out) :: M

        double precision :: t_1,t_2
        double precision :: epsout
        integer :: M_0,loop,info,i
        double precision, dimension(:), allocatable :: res
        integer, dimension(64) :: fpm

        call feastinit(fpm)
        fpm(1) = 1
        write(6,*) fpm(42)
        if (fpm(42)==1)then
            stop
        end if
        !fpm(42) = 0

        M_0 = 10
        if (allocated(eigs)) deallocate(eigs)
        if (allocated(vecs)) deallocate(vecs)
        allocate(res(M_0),eigs(M_0),vecs(A%shape(1),M_0))

        write(6,*) "Finding eigenvalues using FEAST"
        t_1 = omp_get_wtime()
        call zfeast_scsrgv('F',A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
                            ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
        !call zfeast_gcsrgv(A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
        !                    ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
        !call zfeast_sygv('F',A%shape(1),H,A%shape(1),S,A%shape(1) &
        !                    ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
        t_2 = omp_get_wtime()
        write(6,*) "Time for FEAST (s): ", t_2-t_1

        do i = 1,M
            write(6,*) eigs(i)
        end do

    end subroutine FEAST
end module eig_tools