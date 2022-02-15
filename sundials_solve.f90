module sundials_solve

    use, intrinsic :: iso_c_binding
    use fcvode_mod ! CVODE
    use fsundials_nvector_mod ! Generic vector
    use fnvector_serial_mod ! Serial vector
    use fsundials_linearsolver_mod ! Linear system Ax=b 
    use fsundials_matrix_mod ! Matrix
    use fsunmatrix_sparse_mod ! Sparse matrix
    ! use fsunlinsol_klu_mod ! External KLU sparse matrix direct solver methods
    use fsunnonlinsol_fixedpoint_mod ! Non linear fixed point solver
    use fsunnonlinsol_newton_mod ! Non linear fixed point solver
    use fsundials_nonlinearsolver_mod ! Non linear solver
    use fsunmatrix_dense_mod       ! Fortran interface to dense SUNMatrix
    use fsunlinsol_dense_mod       ! Fortran interface to dense SUNLinearSolver
    use fsunlinsol_spbcgs_mod ! Fortran interface to Bi-CGStab iterative linear solver

    ! Gradient calculating ode
    use ODEContainer
    implicit none

    type sundials_solve_class
        type(SUNLinearSolver), pointer :: sunlinsol ! sundials linear solver
        type(SUNNonLinearSolver), pointer :: sunnonlinsol ! sundials non-linear solver
        type(c_ptr) :: cvode_mem ! Memory location for CVODE 
        type(N_Vector), pointer :: solution_vector ! Holds a sundials solution vector
        real(c_double), allocatable :: c_solution_vector(:)
        type(SUNMatrix), pointer :: A ! Generic SUNDIALS matrix class, will be sparse, dense or banded
    contains
        procedure, public, pass :: create_cvode_environment
        procedure, public, pass :: create_sparse_matrix_linear_solver
        procedure, public, pass :: create_dense_matrix_linear_solver
        procedure, public, pass :: create_linear_solver
        procedure, public, pass :: create_non_linear_solver
        procedure, public, pass :: create_root_finding_function
        procedure, public, pass :: reinitialise
        procedure, public, pass :: associate_user_data
        procedure, public, pass :: set_min_max_step
        procedure, public, pass :: solve_timestep
        procedure, public, pass :: get_solution
        procedure, public, pass :: set_solution
        procedure, public, pass :: print_diagnostics
        procedure, public, pass :: destroy => destroy_sundials_solver
    endtype

    private check_association_sunmat
    contains

    subroutine create_cvode_environment(this, method, initial_condition_vector, tstart, rtol, atol, mxsteps)
        !! Allocates memory and creates solver that will use a given method, of a given type, with specified initial conditions and with specified relative and absolute error tolerances
        class(sundials_solve_class), intent(inout) :: this
        character(*), intent(in) :: method ! bdf (backwards differentiation - multi-step fully implicit method) or adams (multi-step explicit method)
        real(c_double), intent(in) :: initial_condition_vector(:) ! initial conditions
        real(c_double), intent(in) :: tstart ! start time
        real(c_double), intent(in), optional :: rtol ! relative error tolerance
        real(c_double), intent(in), optional :: atol ! absolute error tolerance
        integer, intent(in), optional :: mxsteps ! maximum number of steps to be taken by solver, default (500), a value <0 will disable max step test (not recommended, code may get stuck without reporting it)

        ! internal variables
        integer(c_int) :: ierr
        integer(c_long) :: neq
        integer(c_long) :: mxsteps_act
        real(c_double) :: rtol_act = 1.d-6
        real(c_double) :: atol_act = 1.d-6
 
        ! ensures that the solution vector is pointed to c-data allocated within this module so that pointed to data remains clearly visible.
        call this%set_solution(initial_condition_vector)

        neq = size(initial_condition_vector)

        ! create the SUNDIALS solution vector so that it points to the prescribed initial conditions
        this%solution_vector => FN_VMake_Serial(neq, this%c_solution_vector)

        ! create CVode memory
        if(method=='adams') then
            this%cvode_mem = FCVodeCreate(CV_ADAMS)
        elseif(method=='bdf') then
            this%cvode_mem = FCVodeCreate(CV_BDF)
        else
            print *, 'Error in FCVodeCreate, solver "method" must be "adams" or "bdf"'
        endif
        if (.not. c_associated(this%cvode_mem)) then
            print *, 'ERROR: cvode_mem = NULL'
            stop 1
        end if

        ! associate the template rhs function
        ierr = FCVodeInit(this%cvode_mem, c_funloc(RhsFn), tstart, this%solution_vector)
        if (ierr /= 0) then
            print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! if optional tolerances are not specified then assign them
        if(present(rtol)) rtol_act = rtol 
        if(present(atol)) atol_act = atol

        ! set problem toleratnces
        ierr = FCVodeSStolerances(this%cvode_mem, rtol_act, atol_act)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! set the number of maximum steps if argument has been provided
        if(present(mxsteps)) then
            mxsteps_act = mxsteps
            ierr = FCVodeSetMaxNumSteps(this%cvode_mem, mxsteps_act)
            if(ierr /= 0) then
                print*, 'Error in FCVodeSetMaxNumSteps, ierr = ', ierr, '; halting'
                stop 1
            endif
        endif
    end subroutine

    subroutine create_sparse_matrix_linear_solver(this, neq, total_non_zeros, comp_row_col)
        !! The problem can be formed as A y = dy where the matrix A is the linear operator used to produce gradient from y.
        !! We are not required to actually form matrix A, but instead decide on it's type and set it's size for the solver.
        !! This subroutine creates a sparse matrix using either sparse column or row.
        class(sundials_solve_class), intent(inout) :: this
        integer, intent(in) :: neq ! total number of odes that are going to be solved
        integer, intent(in) :: total_non_zeros ! total number of non-zero values in matrix
        character(*), intent(in), optional :: comp_row_col ! either "compressed_row" or "compressed_column", default is to use compressed row.

        integer(c_long) :: neqn
        integer(c_long) :: nnz

        neqn = neq; nnz = total_non_zeros

        call check_association_sunmat(this%A, .false.)

        if(present(comp_row_col)) then
            if(comp_row_col == 'compressed_column') then
                this%A => FSUNSparseMatrix(neqn, neqn, nnz, CSC_MAT)
                call check_association_sunmat(this%A, .true.)
                return
            endif
        endif
        this%A => FSUNSparseMatrix(neqn, neqn, nnz, CSR_MAT)
        call check_association_sunmat(this%A, .true.)
    end subroutine

    subroutine create_root_finding_function(this, num_roots, ode_container)
        !! Initialises the functionality to find roots in the system
        class(sundials_solve_class), intent(inout) :: this
        integer, intent(in) :: num_roots
        type(odesContainerClass), intent(inout), target :: ode_container

        integer(c_int) :: ierr

        ierr = FCVodeRootInit(this%cvode_mem, num_roots, c_funloc(RootFn))
        if (ierr /= 0) then
            print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
            stop 1
        end if
    end subroutine

    subroutine create_dense_matrix_linear_solver(this, neq)
        !! Creates a dense linear matrix, only optimal if neq < 100 approximately, otherwise use banded or sparse.
        !! N.B. sparse requires user supplied Jacobian, dense and banded can estimate Jacobian through difference
        !! quotients (alternatively use iterative linear solvers such as GMRES, Bi-CGSTAG, etc...)
        class(sundials_solve_class), intent(inout) :: this
        integer, intent(in) :: neq

        integer(c_long) :: neqn

        call check_association_sunmat(this%A, .false.)

        neqn = neq

        ! create a dense matrix
        this%A => FSUNDenseMatrix(neqn, neqn)
        if (.not. associated(this%A)) then
            print *, 'ERROR: sunmat = NULL'
            stop 1
        end if
        call check_association_sunmat(this%A, .true.)
    end subroutine

    subroutine create_linear_solver(this, solver_type)
        !! Create the linear solver for solving linear system of odes, depends on the matrix type
        !! that has been specified. E.g. Sparse matrix storage for linear system requires the use of the KLU direct sparse linear solver. Iterative Krylov solvers and
        !! direct dense matrix solvers are also available.
        class(sundials_solve_class), intent(inout) :: this
        character(*), intent(in) :: solver_type
        
        integer(c_int) :: ierr

        
        if(solver_type=='sparse_linear') then
            ! must specify a direct sparse KLU solver to go with sparse matrix storage scheme, this also requires jacobian evaluating methods since not dense or banded matrix.
            call check_association_sunmat(this%A, .true.)
            ! this%sunlinsol => FSUNLinSol_KLU(this%solution_vector, this%A)
        endif

        if(solver_type=='dense') then
            ! Dense linear solver is to be created
            call check_association_sunmat(this%A, .true.)
            this%sunlinsol => FSUNDenseLinearSolver(this%solution_vector, this%A)
        endif

        if(solver_type=='bicgstab') then
            ! A iterative linear solver is specified that doesn't require an A matrix
            this%sunlinsol => FSUNLinSol_SPBCGS(this%solution_vector, PREC_NONE, 5)
        endif

        if (.not. associated(this%sunlinsol)) then
            print *, 'ERROR: this%sunlinsol = NULL'
            stop 1
        end if

        ! attach linear solver
        ierr = FCVodeSetLinearSolver(this%cvode_mem, this%sunlinsol, this%A);
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetLinearSolver, ierr = ', ierr, '; halting'
            stop 1
        end if 
    end subroutine

    subroutine create_non_linear_solver(this, solver_type)
        !! Create a solver that is non-linear. At present can create fixed-point iterative solver that does not require mass matrix or Jacobian.
        class(sundials_solve_class), intent(inout) :: this
        character(*), intent(in) :: solver_type

        integer(c_int) :: ierr

        if(solver_type == 'fixed_point') then
            this%sunnonlinsol => FSUNNonlinSol_FixedPoint(this%solution_vector, 0)
            if (.not. associated(this%sunnonlinsol)) then
                print *,'ERROR: this%sunnonlinsol = NULL'
                stop 001
            end if
        elseif(solver_type == 'newton') then
            ! call check_association_sunmat(this%A, .true.)
            this%sunnonlinsol => FSUNNonlinSol_Newton(this%solution_vector)
            if (.not. associated(this%sunnonlinsol)) then
                print *,'ERROR: this%sunnonlinsol = NULL'
                stop 003
            end if
        endif

        ! attach nonlinear solver object to CVode
        ierr = FCVodeSetNonlinearSolver(this%cvode_mem, this%sunnonlinsol)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetNonlinearSolver, ierr = ', ierr, '; halting'
            stop 002
        end if
    end subroutine

    subroutine reinitialise(this, t, vals)
        !! Associate the user data and gradient calculating function with the ode solver. The argument "ode_container" will vary in type
        !! depending on the problem
        class(sundials_solve_class), intent(inout) :: this
        real(8) :: t
        real(8), intent(in) :: vals(:)

        real(c_double) :: t_c
        integer(c_int) :: ierr
        integer(c_long) :: neq

        call FN_VDestroy(this%solution_vector)

        t_c = t
        ! Re-create nvector with new initialised values
        neq = size(vals)
        this%c_solution_vector = vals
        this%solution_vector => FN_VMake_Serial(neq, this%c_solution_vector)
        

        ierr = FCVodeReInit(this%cvode_mem , t_c, this%solution_vector)
    end subroutine

    subroutine associate_user_data(this, ode_container)
        !! Associate the user data and gradient calculating function with the ode solver. The argument "ode_container" will vary in type
        !! depending on the problem
        class(sundials_solve_class), intent(inout) :: this
        type(odesContainerClass), intent(inout), target :: ode_container

        integer(c_int) :: ierr

        ierr = FCVodeSetUserData(this%cvode_mem , c_loc(ode_container))
    end subroutine

    subroutine set_min_max_step(this, min_step, max_step)
        !! Sets upper and lower bounds for the time step size, useful for ensuring CFL condition is matched and we don't step too large.
        class(sundials_solve_class), intent(inout) :: this
        real(c_double), intent(in) :: min_step
        real(c_double), intent(in) :: max_step

        integer(c_int) :: ierr

        ierr = FCVodeSetMinStep(this%cvode_mem, min_step)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetMinStep, ierr = ', ierr, '; halting'
            stop 1
        end if

        ierr = FCVodeSetMaxStep(this%cvode_mem, max_step)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetMinStep, ierr = ', ierr, '; halting'
            stop 1
        end if
    end subroutine

    subroutine solve_timestep(this, curr_time, end_time, iflag)
        !! Performs a single time-stepped solve and returns certain error codes if unable to reach output time or other errors appear.
        !! Therefore, may not actually reach tout.
        class(sundials_solve_class), intent(inout) :: this
        real(8), intent(inout) :: curr_time
        real(8), intent(inout) :: end_time
        integer(c_int) :: iflag

        real(8) :: arr_curr_time(1)

        arr_curr_time(1) = curr_time

        ! convert the input vector to the sunvector type
        iflag = FCVode(this%cvode_mem, end_time, this%solution_vector, arr_curr_time, CV_NORMAL)
        ! if(iflag /= 0) then
        !     write(*,*) iflag
        !     ! stop 100
        ! endif
        ! if(iflag /= 0) then
        !     write(*,*) 'issue here', iflag
        ! endif

        curr_time = arr_curr_time(1)
    end subroutine

    subroutine get_solution(this, vals)
        !! Obtains the current solution vector that the solver is pointing to. Editing this vector will affect solver!
        class(sundials_solve_class), intent(inout) :: this
        real(8), intent(inout) :: vals(:)

        vals = this%c_solution_vector
    end subroutine

    subroutine set_solution(this, vals)
        !! Sets the current solution vector that the solver is pointing to. Editing this vector will affect solver!
        class(sundials_solve_class), intent(inout) :: this
        real(8), intent(in) :: vals(:)

        if(.not.allocated(this%c_solution_vector)) allocate(this%c_solution_vector(size(vals)))
        this%c_solution_vector = vals
    end subroutine

    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) result(ierr) bind(C,name='RhsFn')
        !! Compute the gradient 
        real(c_double), value :: tn
        type(N_Vector) :: sunvec_y
        type(N_Vector) :: sunvec_f
        type(c_ptr), value :: user_data
        type(odesContainerClass), pointer :: f_user_data => null()

        
        ! C-data types for interface to sundials (c_double is 64-bit, equivalent to fortran double precision so may be casted to one another).
        real(c_double), pointer :: yvec(:) => null()
        real(c_double), pointer :: fvec(:) => null()

        call C_F_POINTER(user_data, f_user_data)

        ! get sundials vectors as conventional data types
        yvec => FN_VGetArrayPointer(sunvec_y)
        fvec => FN_VGetArrayPointer(sunvec_f)

        ! Compute the RHS vector
        call f_user_data%calculate_gradient(tn, yvec(:), fvec(:))

        ! call this%set_solution(yvec)

        ! write(*,*) fvec(1:100)

        yvec => null()
        fvec => null()
        f_user_data => null()

        ierr = 0
        return
    end function

    integer(c_int) function RootFn(tn, sunvec_y, g_out, user_data) result(ierr) bind(C,name='RootFn')
        !! Compute the gradient 
        real(c_double), value :: tn
        type(N_Vector) :: sunvec_y
        real(8) :: g_out(7)
        type(c_ptr), value :: user_data
        type(odesContainerClass), pointer :: f_user_data => null()

        
        ! C-data types for interface to sundials (c_double is 64-bit, equivalent to fortran double precision so may be casted to one another).
        real(c_double), pointer :: yvec(:) => null()

        call C_F_POINTER(user_data, f_user_data)

        g_out = 0.d0

        ! get sundials vectors as conventional data types
        yvec => FN_VGetArrayPointer(sunvec_y)
        ! Compute the RHS vector
        ! write(*,*) g_out
        call f_user_data%root_evaluation(tn, yvec(:), g_out)
        ! write(*,*) g_out

        yvec => null()
        f_user_data => null()

        ierr = 0
        return
    end function

    subroutine destroy_sundials_solver(this)
        !! Free up all memory associated with the SUNDIALS solver
        class(sundials_solve_class), intent(inout) :: this
        
        integer(c_int) :: ierr

        ! clean up
        ! call FCVodeFree(this%cvode_mem)
        ! ierr = FSUNLinSolFree(this%sunlinsol)
        ! call FSUNMatDestroy(this%A)
        ! call FN_VDestroy(this%solution_vector)
        ! if(allocated(this%c_solution_vector)) deallocate(this%c_solution_vector)
    end subroutine

    subroutine print_diagnostics(this)
        !! Prints out statistics associated with the solve that has been performed.
        !! Print statements below clearly describe information that is to be printed.
        class(sundials_solve_class), intent(in) :: this

        integer(c_int)  :: ierr          ! error flag

        integer(c_long) :: nsteps(1)     ! num steps
        integer(c_long) :: nfevals(1)    ! num function evals
        integer(c_long) :: nlinsetups(1) ! num linear solver setups
        integer(c_long) :: netfails(1)   ! num error test fails

        integer(c_int)  :: qlast(1)      ! method order in last step
        integer(c_int)  :: qcur(1)       ! method order for next step

        real(c_double)  :: hinused(1)    ! initial step size
        real(c_double)  :: hlast(1)      ! last step size
        real(c_double)  :: hcur(1)       ! step size for next step
        real(c_double)  :: tcur(1)       ! internal time reached

        integer(c_long) :: nniters(1)    ! nonlinear solver iterations
        integer(c_long) :: nncfails(1)   ! nonlinear solver fails

        ! general solver statistics
        ierr = FCVodeGetIntegratorStats(this%cvode_mem, nsteps, nfevals, nlinsetups, &
        netfails, qlast, qcur, hinused, hlast, hcur, tcur)
        if (ierr /= 0) then
            print *, 'Error in FCVodeGetIntegratorStats, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! nonlinear solver statistics
        ierr = FCVodeGetNonlinSolvStats(this%cvode_mem, nniters, nncfails)
        if (ierr /= 0) then
            print *, 'Error in FCVodeGetNonlinSolvStats, ierr = ', ierr, '; halting'
            stop 1
        end if

        print *, ' '
        print *, ' General Solver Stats:'
        print '(4x,A,i9)'    ,'Total internal steps taken =',nsteps
        print '(4x,A,i9)'    ,'Total rhs function calls   =',nfevals
        print '(4x,A,i9)'    ,'Num lin solver setup calls =',nlinsetups
        print '(4x,A,i9)'    ,'Num error test failures    =',netfails
        print '(4x,A,i9)'    ,'Last method order          =',qlast
        print '(4x,A,i9)'    ,'Next method order          =',qcur
        print '(4x,A,es12.5)','First internal step size   =',hinused
        print '(4x,A,es12.5)','Last internal step size    =',hlast
        print '(4x,A,es12.5)','Next internal step size    =',hcur
        print '(4x,A,es12.5)','Current internal time      =',tcur
        print '(4x,A,i9)'    ,'Num nonlinear solver iters =',nniters
        print '(4x,A,i9)'    ,'Num nonlinear solver fails =',nncfails
        print *, ' '
    end subroutine

    !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !! INTERNAL MODULE PROCEDURES
    !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    subroutine check_association_sunmat(sunmat, logic_associated)
        !! Internal function to check the association of the sun matrix class
        type(SUNMatrix), pointer :: sunmat
        logical, intent(in) :: logic_associated ! supply false if want to ensure not associated and true if want to ensure associated

        if (.not.(associated(sunmat) .eqv. logic_associated)) then
            print *,'ERROR: sunmat = NULL'
            stop 1
        end if
    end subroutine
end module