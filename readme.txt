This program calculates the relative residual norm of the solution of the linear system Ax = b. 
The relative residual norm is visualised in a graph over the number of iterations. 

To run it follow these steps:

1) in terminal run: ./main -a offshore_petsc -b offshore_b_petsc -ksp_monitor ascii:convergence.txt -ksp_monitor   
Explanation:	-a offshore_petsc takes the file offshore_petsc as the matrix A
		-b offshore_b_petsc takes the file offshore_b_petsc as the vector b
		-ksp_monitor ascii:convergence.txt stores the relative residual norms in a file convergence.txt
		-ksp_monitor prints the relative residual norms to the console

2) in terminal run: bash create_plot
Explanation: This creates the file residual_norms.png with gnu plot
