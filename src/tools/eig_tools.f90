module eig_tools
    use kind_tools
    use sparse_array_tools, only: CSR_matrix,CSR_mv,CSR_mv_sym,mv
    use stdlib_linalg_lapack, only: ggev
    use stdlib_linalg_blas, only: gemv
    use stdlib_sorting, only: sort_index
    use iso_fortran_env, only: stdout => output_unit, stderr => error_unit
    use PARDISO_tools, only: PARDISO_solver
    use precond_tools, only: block_PC,block_PQ_PC,block_Jacobi_PC
    use ILU0_tools, only: ILU0
    use GMRES_tools, only: zFGMRES
    use omp_lib, only: omp_get_wtime
    implicit none

    ! complex(dp), external :: zdotu
    interface
        complex(dp) function zdotu(N,zx,incx,zy,incy)
            import :: dp
            integer :: N
            complex(dp) :: zx(*)
            integer :: incx
            complex(dp) :: zy(*)
            integer :: incy
        end function zdotu
    end interface

    integer, parameter :: default_ARPACK_max_iter = 300
    real(dp), parameter :: default_ARPACK_tol = 1e-12_dp

    interface drive_ARPACK_SI
        module procedure :: drive_ARPACK_GMRES, drive_ARPACK_PARDISO, drive_ARPACK_precond
    end interface

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
            vecs(:,i) = vecs(:,i)/sqrt(zdotu(N,vecs(:,i),1,temp_v,1))
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
            vecs(:,i) = vecs(:,i)/sqrt(zdotu(N,vecs(:,i),1,temp_v,1))
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
        complex(dp), dimension(n_b), intent(inout) :: eigs
        complex(dp), dimension(:,:), intent(inout) ::  vecs

        integer, dimension(n_b) :: index
        real(dp), dimension(n_b) :: eigs_real

        eigs_real = real(eigs,kind=dp)

        call sort_index(eigs_real,index)

        eigs = eigs(index)
        vecs = vecs(:,index)
    end subroutine sort_eig

    ! Sort eigenvalues and eigenvector based on the distance from target
    subroutine sort_eig_target(n_b,eigs,vecs,target)
        integer, intent(in) :: n_b
        complex(dp), dimension(n_b), intent(inout) :: eigs
        complex(dp), dimension(:,:), intent(inout) ::  vecs
        complex(dp), intent(in) :: target

        integer, dimension(n_b) :: index
        real(dp), dimension(n_b) :: eigs_real

        eigs_real = abs(eigs-target)

        call sort_index(eigs_real,index)

        eigs = eigs(index)
        vecs = vecs(:,index)
    end subroutine sort_eig_target

#if defined(WITH_FEAST)
    subroutine FEAST(A,B,sym,full,mid,rad,eigs,vecs,M,max_tries)
        type(CSR_matrix), intent(inout) :: A
        type(CSR_matrix), intent(inout) :: B
        logical, intent(in) :: sym
        logical, intent(in) :: full
        double complex, intent(in) :: mid
        double precision, intent(in) :: rad
        double complex, dimension(:), allocatable, intent(inout) :: eigs
        double complex, dimension(:,:), allocatable, intent(inout) :: vecs
        integer, intent(out) :: M
        integer, optional :: max_tries

        double precision :: t_1,t_2
        double precision :: epsout
        integer :: M_0,loop,info,i,max
        double precision, dimension(:), allocatable :: res
        integer, dimension(64) :: fpm
        character :: uplo

        double complex, dimension(:), allocatable :: temp_V

        if (present(max_tries)) then
            max = max_tries
        else
            max = 3
        end if

        if (full) then
            uplo = 'F'
        else
            uplo = 'U'
        end if

        call feastinit(fpm)
        fpm(1) = 1 ! Feast prints info to screen
        fpm(14) = 2 ! Feast does a stochastic estimate of the number of enclosed eigenvalues

        M_0 = 5
        if (allocated(eigs)) deallocate(eigs)
        if (allocated(vecs)) deallocate(vecs)
        allocate(res(M_0),eigs(M_0),vecs(A%shape(1),M_0))

        write(stdout,*) "Estimating number of eigenvalues using FEAST"
        t_1 = omp_get_wtime()
        if (sym) then
            call zfeast_scsrgv(uplo,A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
                            ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
        else
            call zfeast_gcsrgv(A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
                               ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
        end if
        t_2 = omp_get_wtime()
        write(stdout,*) "Time for FEAST eigenvalue estimation (s): ", t_2-t_1

        if (M>10) then
            M_0 = int(M*1.5_dp)
        else
            M_0 = 2*M
        end if

        write(stdout,'(a,i0)')  'FEAST estimate for # of enclosed eigenvalues:', M

        write(stdout,*) "Finding eigenvalues using FEAST"
        t_1 = omp_get_wtime()
        i = 0
        do while (i <= max)
            call feastinit(fpm)
            fpm(1) = 1 ! Feast prints info to screen
            deallocate(res,eigs,vecs)
            allocate(res(M_0),eigs(M_0),vecs(A%shape(1),M_0))
            if (sym) then
                call zfeast_scsrgv(uplo,A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
                                ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
            else
                call zfeast_gcsrgv(A%shape(1),A%data,A%index_ptr,A%indices,B%data,B%index_ptr,B%indices &
                                ,fpm,epsout,loop,mid,rad,M_0,eigs,vecs,M,res,info)
            end if

            if (info == 3) then
                i = i + 1
                if (M_0>10) then
                    M_0 = int(M_0*1.5_dp)
                else
                    M_0 = 2*M_0
                end if
                write(stdout,'(a)') 'FEAST exited with info code 3 (work subspace too small)'
                write(stdout,'(a,i0)') 'Increasing M_0 and retrying. New M_0: ', M_0
            else
                exit
            end if

        end do
        t_2 = omp_get_wtime()
        write(stdout,*) "Time for FEAST (s): ", t_2-t_1

        do i = 1,M
            write(stdout,*) eigs(i)
        end do

        allocate(temp_v(B%shape(1)))
        ! Normalize the eigenvectors.
        ! It seems like FEAST already normalizes them correctly, but I don't see anything in the docs that mentions this.
        do i=1,M
            ! Should use a sparse gemv routine here!
            ! call gemv('N',B%s,N,dcmplx(1.d0,0.d0),B,N,vecs(:,i),1,dcmplx(0.d0,0.d0),temp_v,1)
            if (full) then
                call CSR_mv(B,vecs(:,i),temp_v)
            else
                call CSR_mv_sym(B,vecs(:,i),temp_v)
            end if
            vecs(:,i) = vecs(:,i)/sqrt(zdotu(B%shape(1),vecs(:,i),1,temp_v,1))
        end do

    end subroutine FEAST
#endif

    subroutine drive_ARPACK_PARDISO(A,solver,B,full,shift,n_eigs,eigs,vecs,v0,max_iter,tol)
        type(CSR_matrix), intent(in) :: A
        type(PARDISO_solver), intent(inout) :: solver
        type(CSR_matrix), intent(in) :: B
        logical, intent(in) :: full
        complex(dp), intent(in) :: shift
        integer, intent(in) :: n_eigs
        complex(dp), intent(out) :: eigs(:)
        complex(dp), intent(out) :: vecs(:,:)
        complex(dp), optional, intent(in) :: v0(:)
        integer, optional, intent(in) :: max_iter
        real(dp), optional, intent(in) :: tol

        ! Actual values of optional agruments
        complex(dp) :: resid(B%shape(1))
        integer :: max_i
        real(dp) :: r_tol

        ! Variables used by znaupd
        integer :: info, ido, n, ncv, ldv, lworkl
        character :: bmat
        character(len=2) :: which
        complex(dp), allocatable :: V(:,:)
        integer :: iparam(11), ipntr(14)
        complex(dp), allocatable :: workd(:)
        complex(dp), allocatable :: workl(:)
        real(dp), allocatable :: rwork(:)

        integer :: x_1,x_2,y_1,y_2
        complex(dp), allocatable :: temp(:)

        ! Select correct routine for B matvecs
        procedure(mv), pointer :: mv_B => null()
        if (full) then
            mv_B => CSR_mv
        else
            mv_B => CSR_mv_sym
        end if

        ! Set optional arguments
        if (present(v0)) then
            resid = v0
            info = 1 ! User supplied starting vector
        else
            resid = 0
            info = 0 ! Random starting vector
        end if

        if (present(max_iter)) then
            max_i = max_iter
        else
            max_i = default_ARPACK_max_iter
        end if

        if (present(tol)) then
            r_tol = tol
        else
            r_tol = default_ARPACK_tol
        end if

        ! Set variables used by RCI calls
        ido = 0
        bmat = 'G'
        n = B%shape(1)
        which = 'LM'
        ncv = 2*n_eigs ! Recomended by ARPACK documentation
        allocate(V(n,ncv))
        ldv = n
        iparam(1) = 1
        iparam(3) = max_i
        iparam(4) = 1
        iparam(7) = 3
        allocate(workd(3*n))
        lworkl = 3*ncv**2 + 5*ncv
        allocate(workl(lworkl))
        allocate(rwork(ncv))

        ! Allocate temp vector used during RCI
        allocate(temp(n))

        write(stdout,*) ""
        write(stdout,*) "Finding eigenvalues with znaupd..."

        do while (ido /= 99)
            call znaupd(ido, bmat, n, which, n_eigs, r_tol, resid, ncv, V, ldv,&
                        iparam, ipntr, workd, workl, lworkl, rwork, info)

            if (ido == -1) then
                ! Perform (A-shift*B)^-1*B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),temp)
                call solver%solve(A%data,A%index_ptr,A%indices,workd(y_1:y_2),temp)

            else if (ido == 1) then
                ! Perform (A-shift*B)^-1*B*x = y, with B*x stored in workd(ipntr(3))
                x_1 = ipntr(3)
                x_2 = ipntr(3) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call solver%solve(A%data,A%index_ptr,A%indices,workd(y_1:y_2),workd(x_1:x_2))

            else if (ido == 2) then
                ! Perform B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),workd(y_1:y_2))

            else if (ido == 3) then
                write(stderr,*) "ido = 3 returned by znaupd, but user supplied shifts not implemented."
                stop
            end if
        end do

        call postprocess_ARPACK(eigs,vecs,B,mv_B,temp,V,ldv,shift,bmat,n,which,n_eigs,r_tol,resid,ncv,iparam,ipntr,&
                                workd,workl,lworkl,rwork,info)
    end subroutine drive_ARPACK_PARDISO

    subroutine drive_ARPACK_precond(precond,B,full,shift,n_eigs,eigs,vecs,v0,max_iter,tol)
        class(block_PC), intent(inout) :: precond
        type(CSR_matrix), intent(in) :: B
        logical, intent(in) :: full
        complex(dp), intent(in) :: shift
        integer, intent(in) :: n_eigs
        complex(dp), intent(out) :: eigs(:)
        complex(dp), intent(out) :: vecs(:,:)
        complex(dp), optional, intent(in) :: v0(:)
        integer, optional, intent(in) :: max_iter
        real(dp), optional, intent(in) :: tol

        ! Actual values of optional agruments
        complex(dp) :: resid(B%shape(1))
        integer :: max_i
        real(dp) :: r_tol

        ! Variables used by znaupd
        integer :: info, ido, n, ncv, ldv, lworkl
        character :: bmat
        character(len=2) :: which
        complex(dp), allocatable :: V(:,:)
        integer :: iparam(11), ipntr(14)
        complex(dp), allocatable :: workd(:)
        complex(dp), allocatable :: workl(:)
        real(dp), allocatable :: rwork(:)

        integer :: x_1,x_2,y_1,y_2
        complex(dp), allocatable :: temp(:)

        ! Select correct routine for B matvecs
        procedure(mv), pointer :: mv_B => null()
        if (full) then
            mv_B => CSR_mv
        else
            mv_B => CSR_mv_sym
        end if

        ! Set optional arguments
        if (present(v0)) then
            resid = v0
            info = 1 ! User supplied starting vector
        else
            resid = 0
            info = 0 ! Random starting vector
        end if

        if (present(max_iter)) then
            max_i = max_iter
        else
            max_i = default_ARPACK_max_iter
        end if

        if (present(tol)) then
            r_tol = tol
        else
            r_tol = default_ARPACK_tol
        end if

        ! Set variables used by RCI calls
        ido = 0
        bmat = 'G'
        n = B%shape(1)
        which = 'LM'
        ncv = 2*n_eigs ! Recomended by ARPACK documentation
        allocate(V(n,ncv))
        ldv = n
        iparam(1) = 1
        iparam(3) = max_i
        iparam(4) = 1
        iparam(7) = 3
        allocate(workd(3*n))
        lworkl = 3*ncv**2 + 5*ncv
        allocate(workl(lworkl))
        allocate(rwork(ncv))

        ! Allocate temp vector used during RCI
        allocate(temp(n))

        write(stdout,*) ""
        write(stdout,*) "Finding eigenvalues with znaupd..."

        do while (ido /= 99)
            call znaupd(ido, bmat, n, which, n_eigs, r_tol, resid, ncv, V, ldv,&
                        iparam, ipntr, workd, workl, lworkl, rwork, info)

            if (ido == -1) then
                ! Perform (A-shift*B)^-1*B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),temp)
                call precond%solve(workd(y_1:y_2),temp)

            else if (ido == 1) then
                ! Perform (A-shift*B)^-1*B*x = y, with B*x stored in workd(ipntr(3))
                x_1 = ipntr(3)
                x_2 = ipntr(3) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call precond%solve(workd(y_1:y_2),workd(x_1:x_2))

            else if (ido == 2) then
                ! Perform B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),workd(y_1:y_2))

            else if (ido == 3) then
                write(stderr,*) "ido = 3 returned by znaupd, but user supplied shifts not implemented."
                stop
            end if
        end do

        call postprocess_ARPACK(eigs,vecs,B,mv_B,temp,V,ldv,shift,bmat,n,which,n_eigs,r_tol,resid,ncv,iparam,ipntr,&
                                workd,workl,lworkl,rwork,info)
    end subroutine drive_ARPACK_precond

    subroutine drive_ARPACK_GMRES(GMRES,precond,B,full,shift,n_eigs,eigs,vecs,v0,max_iter,tol)
        type(zFGMRES), intent(inout) :: GMRES
        class(*), intent(inout) :: precond
        type(CSR_matrix), intent(in) :: B
        logical, intent(in) :: full
        complex(dp), intent(in) :: shift
        integer, intent(in) :: n_eigs
        complex(dp), intent(out) :: eigs(:)
        complex(dp), intent(out) :: vecs(:,:)
        complex(dp), optional, intent(in) :: v0(:)
        integer, optional, intent(in) :: max_iter
        real(dp), optional, intent(in) :: tol

        ! Actual values of optional agruments
        complex(dp) :: resid(B%shape(1))
        integer :: max_i
        real(dp) :: r_tol

        ! Variables used by znaupd
        integer :: info, ido, n, ncv, ldv, lworkl
        character :: bmat
        character(len=2) :: which
        complex(dp), allocatable :: V(:,:)
        integer :: iparam(11), ipntr(14)
        complex(dp), allocatable :: workd(:)
        complex(dp), allocatable :: workl(:)
        real(dp), allocatable :: rwork(:)

        integer :: x_1,x_2,y_1,y_2
        complex(dp), allocatable :: temp(:)

        ! Select correct routine for B matvecs
        procedure(mv), pointer :: mv_B => null()
        if (full) then
            mv_B => CSR_mv
        else
            mv_B => CSR_mv_sym
        end if

        ! Set optional arguments
        if (present(v0)) then
            resid = v0
            info = 1 ! User supplied starting vector
        else
            resid = 0
            info = 0 ! Random starting vector
        end if

        if (present(max_iter)) then
            max_i = max_iter
        else
            max_i = default_ARPACK_max_iter
        end if

        if (present(tol)) then
            r_tol = tol
        else
            r_tol = default_ARPACK_tol
        end if

        ! Set variables used by RCI calls
        ido = 0
        bmat = 'G'
        n = B%shape(1)
        which = 'LM'
        ncv = 2*n_eigs ! Recomended by ARPACK documentation
        allocate(V(n,ncv))
        ldv = n
        iparam(1) = 1
        iparam(3) = max_i
        iparam(4) = 1
        iparam(7) = 3
        allocate(workd(3*n))
        lworkl = 3*ncv**2 + 5*ncv
        allocate(workl(lworkl))
        allocate(rwork(ncv))

        ! Allocate temp vector used during RCI
        allocate(temp(n))

        write(stdout,*) ""
        write(stdout,*) "Finding eigenvalues with znaupd..."

        do while (ido /= 99)
            call znaupd(ido, bmat, n, which, n_eigs, r_tol, resid, ncv, V, ldv,&
                        iparam, ipntr, workd, workl, lworkl, rwork, info)

            if (ido == -1) then
                ! Perform (A-shift*B)^-1*B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),temp)
                call GMRES%solve(workd(y_1:y_2),temp,precond)

            else if (ido == 1) then
                ! Perform (A-shift*B)^-1*B*x = y, with B*x stored in workd(ipntr(3))
                x_1 = ipntr(3)
                x_2 = ipntr(3) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call GMRES%solve(workd(y_1:y_2),workd(x_1:x_2),precond)

            else if (ido == 2) then
                ! Perform B*x = y
                x_1 = ipntr(1)
                x_2 = ipntr(1) + n - 1
                y_1 = ipntr(2)
                y_2 = ipntr(2) + n - 1
                call mv_B(B,workd(x_1:x_2),workd(y_1:y_2))

            else if (ido == 3) then
                write(stderr,*) "ido = 3 returned by znaupd, but user supplied shifts not implemented."
                stop
            end if
        end do

        call postprocess_ARPACK(eigs,vecs,B,mv_B,temp,V,ldv,shift,bmat,n,which,n_eigs,r_tol,resid,ncv,iparam,ipntr,&
                                workd,workl,lworkl,rwork,info)
    end subroutine drive_ARPACK_GMRES

    subroutine postprocess_ARPACK(eigs,vecs,B,mv_B,temp,V,ldv,shift,bmat,n,which,n_eigs,r_tol,resid,ncv,iparam,ipntr,&
                                workd,workl,lworkl,rwork,info)
        complex(dp), intent(out) :: eigs(:)
        complex(dp), intent(out) :: vecs(:,:)
        type(CSR_matrix), intent(in) :: B
        complex(dp), intent(inout) :: temp(:)
        procedure(mv), pointer, intent(in) :: mv_B
        complex(dp), intent(inout) :: V(:,:)
        integer, intent(in) :: ldv
        complex(dp), intent(in) :: shift
        character, intent(in) :: bmat
        integer, intent(in) :: n
        character(len=2), intent(in) :: which
        integer, intent(in) :: n_eigs
        real(dp), intent(in) :: r_tol
        complex(dp), intent(inout) :: resid(:)
        integer, intent(in) :: ncv
        integer, intent(in) :: iparam(11)
        integer, intent(inout) :: ipntr(14)
        complex(dp), intent(in) :: workd(:)
        complex(dp), intent(inout) :: workl(:)
        integer, intent(in) :: lworkl
        real(dp), intent(in) :: rwork(:)
        integer, intent(inout) :: info

        ! Variables used by zneupd
        logical :: rvec
        character :: howmny
        logical, allocatable :: select(:)
        integer :: nconv
        complex(dp), allocatable :: d(:)
        complex(dp), allocatable :: workev(:)

        integer :: i

        if (info /= 0) then
            ! ARPACK did not converge, or some other error occurred
            write(stderr,'(a,i0)') "znaupd exited with info = ", info
            write(stderr,*) "Please consult znaupd documentation in ARPACK-ng for explanation"
            stop
        end if

        ! ARPACK converged, extract eigenvalues and eigenvectors
        nconv = iparam(5)
        rvec = .true.
        howmny = 'A'
        allocate(select(ncv),d(nconv),workev(2*ncv))
        call zneupd(rvec, howmny, select, d, V, ldv, shift, workev, bmat, n, which, n_eigs, r_tol, resid, &
                    ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, rwork, info)

        if (info /= 0) then
            ! Some error occurred in zneupd
            write(stderr,'(a,i0)') "zneupd exited with info = ", info
            write(stderr,*) "Please consult zneupd documentation in ARPACK-ng for explanation"
            stop
        end if

        eigs = d(1:n_eigs)
        vecs = V(:,1:n_eigs)

        ! Apply B normalization to eigenvectors (complex-syms)
        do i=1,n_eigs
            call mv_B(B,vecs(:,i),temp)
            vecs(:,i) = vecs(:,i)/sqrt(zdotu(n,vecs(:,i),1,temp,1))
        end do

        write(stdout,*) "Done!"
        write(stdout,*) "Eigenvalues: "
        do i = 1,n_eigs
            write(stdout,*) eigs(i)
        end do
        write(stdout,*) ""

    end subroutine postprocess_ARPACK
end module eig_tools