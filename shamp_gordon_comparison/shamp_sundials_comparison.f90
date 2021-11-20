program sundials_test
	use ShampGordon
	use sundials_solve
	implicit none


	type(sundials_solve_class) :: solver
	type(ode_test) :: odes
	integer :: neqn = 3
	integer :: iflag
	real(8) :: curr_time, end_time
	real(8) :: initial_conditions(3)
	real(8) :: solution(3)
	real(8) :: relerr, abserr
	integer :: iwork(5)
	REAL(8), ALLOCATABLE, DIMENSION(:) :: work
	real(8) :: cpu_start, cpu_end
	real(8) :: sim_end

	initial_conditions = 1.d0
	neqn = 3
	curr_time = 0.1d0
	relerr= 1.d-5
	abserr = 1.d-5
	sim_end = 1.d2
	ALLOCATE( work( 100+21*neqn ) )

	write(*,*) 'hi'

	write(*,*) 'Creating cvode environment.'
	call solver%create_cvode_environment('adams', initial_conditions, curr_time, relerr, abserr)
	write(*,*) 'environment created.'

	call odes%initialise(neqn)

	! write(*,*) 'Creating sparse matrix.'
	! call solver%create_sparse_matrix_linear_solver(neq=neqn, total_non_zeros=9, comp_row_col='compressed_row')
	! write(*,*) 'Sparse matrix created.'
	! write(*,*) 'Creating linear solver.'
	! call solver%create_linear_solver('sparse_linear')
	! write(*,*) 'Linear solver created.'

	! write(*,*) 'Creating dense matrix.'
	! call solver%create_dense_matrix_linear_solver(neqn)
	! write(*,*) 'Created densse matrix.'
	! write(*,*) 'Creating linear solver.'
	! call solver%create_linear_solver('dense')
	! write(*,*) 'Linear solver created.'

	write(*,*) 'Creating linear solver.'
	call solver%create_linear_solver('bicgstab')
	write(*,*) 'Linear solver created.'

	write(*,*) 'Creating non-linear solver.'
	call solver%create_non_linear_solver('newton')
	write(*,*) 'Non-linear solver created.'

	write(*,*) 'Associate user data.'
	call solver%associate_user_data(odes)
	write(*,*) 'Associated.'

	write(*,*) 'Solving ODE.'
	! open(145, file='lorenz_out.dat')

	! Perform SUNDIALS solve
	call cpu_time(cpu_start)
	call solver%set_min_max_step(0.0d0, 1.d0)
	do while (curr_time < sim_end)
		end_time = min(curr_time + 1.d-2, sim_end)
		iflag = 0
		call solver%solve_timestep(curr_time, end_time, iflag)
		! call solver%get_solution(solution)
		! write(*,*) curr_time, solution, iflag
		! write(145,*) solution
	enddo
	call cpu_time(cpu_end)
	! Performed solve.
	write(*,*) 'Solved using SUNDIALS ODE in: ', cpu_end - cpu_start

	write(*,*) 'Solving problem using Shampine Gordon'
	call cpu_time(cpu_start)
	! Perform Shampine Gordon Solve
	curr_time = 0.1d0
	initial_conditions = 1.d0
	relerr= 1.d-5
	abserr = 1.d-5
	do while (curr_time < sim_end)
		end_time = min(curr_time + 1.d-5, sim_end)
		iflag = 1
		call shampODE ( odes, 3, initial_conditions, curr_time, end_time, relerr, abserr, iflag, work, iwork )
		write(*,*) curr_time, initial_conditions, iflag
	enddo
	call cpu_time(cpu_end)
	! Performed Solve.
	write(*,*) 'Solved problem using Shampine Gordon in: ', cpu_end - cpu_start
	! close(145)
	! call solver%print_diagnostics()
	
end program