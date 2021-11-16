# sundials_interface

Interface is created between custom configure and build of the SUNDIALS CVODE solver.

At present the solution of non-linear ODEs using BDF or Adams-Moulton fixed-point multi-step methods are allowed.

Example problem has been created for the solution of the Lorenz system of equations.

To solve a new ODE it is simply a matter of creating new derived type with matching functions as the the module "dummy_grad.f90" type "ode_test". Updates to module "sundials_solve.f90" functions "RhsFn" is required to specify the new type as the Fortran user data argument "f_user_data" so that the C_F_POINTER intrinsic function operates correctly. Corresponding changes may be required to the use statement for the module containing the ode type at the top of "sundials_solve.f90". The subroutine "sundials_solve_class::associate_user_data" should also read in the new derived type used to calculate the ode gradients.
