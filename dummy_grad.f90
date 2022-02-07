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

		! Lorenz attractor
		! rho = 15.d0; sigma = 10.d0; beta = 8.d0/3.d0;
		! dy(1) = sigma * (y(2) - y(1))
		! dy(2) = y(1) * (rho - y(3)) - y(2)
		! dy(3) = y(1) * y(2) - beta * y(3)
		! write(*,*) dy, y

		! Chemical balance equations
		dy(1) = -0.04d0 * y(1) + 1.d4 * y(2) * y(3)
		dy(2) = 0.04 * y(1) - 1.d4 * y(2) * y(3) - 3.d7 * y(2)**2
		dy(3) = 3.d7 * y(2) **2

		write(*,*) 'dy', dy
	end subroutine


end module