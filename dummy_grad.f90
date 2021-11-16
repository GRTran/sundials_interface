module dummy_grad
implicit none

	type :: ode_test
		integer :: neqn
	contains
		procedure, public, pass :: initialise => initialise_ode
		procedure, public, pass :: calculate_gradient
	endtype

	contains

	subroutine initialise_ode(this, neqn)
		class(ode_test), intent(inout) :: this
		integer, intent(in) :: neqn
		this%neqn = neqn
	end subroutine

	subroutine calculate_gradient(this, tn, y, dy)	
		class(ode_test), intent(inout) :: this
		real(8), intent(inout) :: tn
		real(8), intent(inout) :: y(this%neqn)
		real(8), intent(inout) :: dy(this%neqn)

		real(8) :: rho, sigma, beta

		rho = 20.d0; sigma = 10.d0; beta = 8.d0/3.d0;

		! write(*,*) rho, sigma, beta
		dy(1) = sigma * (y(2) - y(1))
		dy(2) = y(1) * (rho - y(3)) - y(2)
		dy(3) = y(1) * y(2) - beta * y(3)
		! write(*,*) dy, y
	end subroutine


end module