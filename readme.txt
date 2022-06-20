This program solves the linear system Ax = b. 

You can monitor the absolute residual norm, or calculate the runtime of the solver in milliseconds.


To run it follow these steps:

1) in terminal run: 
mpirun -n <number_of_processes> ./main -a <matrixfilename> -b <vectorfilename>

- Explanation
The <number_of_processes> is a number between 1 and 48 (on Harris)
<matrixfilename> is the filename of the matrix in the "vectors_and_matrices" directory.
<vectorfilename> is the filename of the right hand side vector in the "vectors_and_matrices" directory.

If the flag -b is omitted then the right hand side vector will be created with random values between 0 and 1 (uniform distribution)

By default GMRES with GAMG preconditioning is used to solve the system.
The methods for the solver and the preconditioning can be changed with the flags -ksp_type <ksp_type> and -pc_type <pc_type> respectively. 
Values for <ksp_type> can be found here: https://petsc.org/main/docs/manualpages/KSP/KSPType/
Values for <pc_type> can be found here: https://petsc.org/main/docs/manualpages/PC/PCType/

By default the program monitors the residual norm in each iteration. This value is printed to the console and also stored in a file called "rnorm_<matrixname>_<ksp_type>_<pc_type>.txt". 

By default the program also prints the termination reason to the console. This is an integer number, a description of the meanings can be found here: https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/

If the flag -time is given, then the program does not do any monitoring (and also does not print the termination reason), but measures the runtime of the solver and prints it to the console. (In Milliseconds) The runtime is measured as average over 5 runs.

It is also important to set the tolerances of the ksp solver. By default they are rtol=1e-5, atol=1e-50, dtol=1e5, and maxits=1e4.
You can change the tolerances with -ksp_rtol <value>,  -ksp_atol <value>,  -ksp_divtol <value> and -ksp_max_it <value>.




- Example 
Solving the offshore matrix with the given right hand side, the cg method and Jacobi preconditioner on 24 processes with rtol=1e-50:

mpirun -n 24 ./main -a offshore_petsc -b offshore_b_petsc -ksp_type cg -pc_type jacobi -ksp_rtol 1e-50

Measuring the runtime of the Bump_2911 matrix with random right hand side, gmres method and Jacobi preconditioner on 16 processes:

mpirun -n 16 ./main -a Bump_2911_petsc -ksp_type gmres -pc_type jacobi -time



2) To create a plot run 
bash create_plot

Make sure that in the bash script the right .txt file is chosen. Also make sure that the for loop iterates over the correct amount of runs. E.G. if you only want to create a plot of one output from the program then set the for loop to one iteration.
