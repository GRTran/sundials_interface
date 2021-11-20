program sundials_test
	use sundials_solve
	implicit none


	type(sundials_solve_class) :: solver
	type(ode_test) :: odes
	integer :: neqn = 3
	integer :: iflag
	real(8) :: curr_time, end_time
	real(8) :: initial_conditions(3)
	real(8) :: solution(3)

	initial_conditions = 1.d0
	neqn = 3
	curr_time = 0.1d0

	write(*,*) 'hi'

	write(*,*) 'Creating cvode environment.'
	call solver%create_cvode_environment('adams', initial_conditions, curr_time, 1.d-5, 1.d-5)
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

	! write(*,*) 'Creating linear solver.'
	! call solver%create_linear_solver('bicgstab')
	! write(*,*) 'Linear solver created.'

	write(*,*) 'Creating non-linear solver.'
	call solver%create_non_linear_solver('fixed_point')
	write(*,*) 'Non-linear solver created.'

	write(*,*) 'Associate user data.'
	call solver%associate_user_data(odes)
	write(*,*) 'Associated.'

	write(*,*) 'Solving ODE.'
	open(145, file='lorenz_out.dat')
	call solver%set_min_max_step(0.0d0, 1.d0)
	do while (curr_time < 1000.d0)
		end_time = min(curr_time + 1.d-2, 1000.d0)
		iflag = 0
		call solver%solve_timestep(curr_time, end_time, iflag)
		call solver%get_solution(solution)
		! write(*,*) curr_time, solution, iflag
		write(145,*) solution
	enddo
	close(145)
	write(*,*) 'Solved ODE'
	call solver%print_diagnostics()
	
end program