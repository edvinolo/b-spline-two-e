module SLEPc_tools
    #include <slepc/finclude/slepceps.h>
    use slepceps
    use kind_tools
    use sparse_array_tools, only: CSR_matrix
    use eig_tools, only: default_EPS_tol => default_ARPACK_tol

    implicit none
    PetscErrorCode :: ierr
    PetscMPIInt :: size,rank

    EPS :: eps
    ST :: st
    KSP :: ksp
    PC :: pc
    Mat :: F

    real(dp), parameter :: default_KSP_tol = 1.0e-14_dp


    #if defined(PETSC_USE_64BIT_INDICES)
        #define MPIU_INT MPI_Long
    #else
        #define MPIU_INT MPI_INT
    #endif

    contains

    subroutine CSR_to_PETSc(csr_in,full,A)
        type(CSR_matrix), intent(inout) :: CSR_in
        logical, intent(in) :: full
        Mat, intent(inout) :: A

        integer :: i,k
        integer, allocatable :: ja(:)

        write(6,*) 'Hello from rank: ', rank

        call CSR_bcast_metadata(CSR_in)

        PetscCallA(MatCreate(PETSC_COMM_WORLD,A,ierr))
        PetscCallA(Matsetsizes(A,PETSC_DECIDE,PETSC_DECIDE,CSR_in%shape(1),CSR_in%shape(2),ierr))

        if (full) then
            PetscCallA(MatSetType(A,MATMPIAIJ,ierr))
            !PetscCallA(MatMPIAIJSetPreallocationCSR(A,CSR%index_ptr,CSR_in%indices,CSR_in%data,ierr))
        else
            PetscCallA(MatSetType(A,MATMPISBAIJ,ierr))
            !PetscCallA(MatMPISBAIJSetPreallocationCSR(A,CSR%index_ptr,CSR_in%indices,CSR_in%data,ierr))
        end if

        PetscCallA(MatSetFromOptions(A,ierr))
        PetscCallA(MatSetUp(A,ierr))

        if (rank == 0) then
            allocate(ja(CSR_in%nnz))

            !ia = CSR_in%index_ptr - 1
            ja = CSR_in%indices - 1

            write(6,*) 'Assembly of A ...'
            do i = 1,CSR_in%ptr_size()-1
                k = CSR_in%index_ptr(i+1)-CSR_in%index_ptr(i)
                PetscCallA(MatSetValues(A,1,[i-1],k,ja(CSR_in%index_ptr(i):CSR_in%index_ptr(i+1)-1),&
                        CSR_in%data(CSR_in%index_ptr(i):CSR_in%index_ptr(i+1)-1),INSERT_VALUES,ierr))
                ! do k = CSR_in%index_ptr(i),CSR_in%index_ptr(i+1)-1
                !     if ((ja(k))>= CSR_in%shape(2)) write(6,*) i-1,ja(k), CSR_in%shape
                !     PetscCallA(MatSetValues(A,1,[i-1],1,[ja(k)],[CSR_in%data(k)],INSERT_VALUES,ierr))
                ! end do
            end do
            write(6,*) 'Done!'
        end if
        PetscCallA(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr))
        write(6,*) 'Hello again from rank: ', rank
        PetscCallA(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr))
    end subroutine CSR_to_PETSc

    subroutine CSR_bcast_metadata(A)
        type(CSR_matrix), intent(inout) :: A

        PetscMPIInt :: bcast_count,root

        bcast_count = 2
        root = 0
        PetscCallMPIA(MPI_Bcast(A%shape,bcast_count,MPIU_INT,root,PETSC_COMM_WORLD,ierr))

        bcast_count = 1
        PetscCallMPIA(MPI_Bcast(A%nnz,bcast_count,MPIU_INT,root,PETSC_COMM_WORLD,ierr))
    end subroutine CSR_bcast_metadata

    subroutine setup_SLEPc_EPS(A,B,tol)
        Mat, intent(in) :: A
        Mat, intent(in) :: B
        real(dp), intent(in), optional :: tol

        real(dp) :: r_tol

        integer :: i

        PetscReal :: t1,t2
        PetscInt :: nconv
        PetscScalar :: eig_i
        Vec :: vec_i

        if (present(tol)) then
            r_tol = tol
        else
            r_tol = default_EPS_tol
        end if

        ! Create EPS object and set options
        PetscCallA(EPSCreate(PETSC_COMM_WORLD,eps,ierr))
        PetscCallA(EPSSetOperators(eps,A,B,ierr))
        PetscCallA(EPSSetProblemType(eps,EPS_PGNHEP,ierr))
        PetscCallA(EPSSetTolerances(eps,r_tol,PETSC_DETERMINE_INTEGER,ierr))
        PetscCallA(EPSSetDimensions(eps,2,PETSC_DETERMINE_INTEGER,PETSC_DETERMINE_INTEGER,ierr))

        PetscCallA(EPSSetTarget(eps,(-0.5d0,1.0d-5),ierr))
        PetscCallA(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE,ierr))

        ! Get the spectral transformation object
        PetscCallA(EPSGetST(eps,st,ierr))
        PetscCallA(STSetType(st,STSINVERT,ierr))

        ! Get and set the KSP object
        PetscCallA(STGetKSP(st,ksp,ierr))
        PetscCallA(KSPSetType(ksp,KSPGMRES,ierr))

        ! Get and set the PC object
        PetscCallA(KSPGetPC(ksp,pc,ierr))
        PetscCallA(PCSetType(pc,PCCHOLESKY,ierr))
        PetscCallA(PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr))

        ! Set MUMPS options
        PetscCallA(STGetOperator(st,PETSC_NULL_MAT,ierr))
        PetscCallA(PCFactorSetUpMatSolverType(pc,ierr))
        PetscCallA(PCFactorGetMatrix(pc,F,ierr))
        PetscCallA(MatMumpsSetIcntl(F,35,1,ierr))
        PetscCallA(MatMumpsSetCntl(F,7,1.0d-5,ierr))
        PetscCallA(MatMumpsSetIcntl(F,37,1,ierr))
        PetscCallA(STRestoreOperator(st,PETSC_NULL_MAT,ierr))

        ! Final KSP setup
        PetscCallA(KSPSetTolerances(ksp,default_KSP_tol,PETSC_DETERMINE_REAL,PETSC_DETERMINE_REAL,PETSC_DETERMINE_INTEGER,ierr))
        PetscCallA(KSPSetFromOptions(ksp,ierr))

        ! Final ESP setup
        PetscCallA(EPSSetFromOptions(eps,ierr))
        PetscCallA(EPSSetup(eps,ierr))

        ! Solve the eigenvalue problem and print info
        if (rank == 0) write(6,*) 'Solving eigenproblem ...'
        t1 = MPI_Wtime()
        PetscCallA(EPSSolve(eps,ierr))
        t2 = MPI_Wtime()
        if (rank == 0) write(6,*) "Done! Walltime for eigenproblem (s): ",t2-t1
        ! PetscCallA(EPSValuesView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr))
        ! PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL,ierr))
        PetscCallA(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD,ierr))
        PetscCallA(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD,ierr))

        PetscCallA(EPSGetConverged(eps,nconv,ierr))
        PetscCallA(MatCreateVecs(A,vec_i,PETSC_NULL_VEC,ierr))
        do i  = 0,nconv-1
            PetscCallA(EPSGetEigenpair(eps,i,eig_i,PETSC_NULL_SCALAR,vec_i,PETSC_NULL_VEC,ierr))
            if (rank == 0) write(6,*) i,eig_i
        end do

        ! Free work space.  All PETSc objects should be destroyed when they
        ! are no longer needed.
        PetscCallA(MatDestroy(A,ierr))
        PetscCallA(MatDestroy(B,ierr))
        PetscCallA(EPSDestroy(eps,ierr))
        PetscCallA(VecDestroy(vec_i,ierr))
    end subroutine setup_SLEPc_EPS

end module SLEPc_tools