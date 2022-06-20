
static char help[] = "Reads a PETSc matrix and vector from a file.\n\
  -a <input_file> : first file containing the matrix\n\
  -b <input_file> : first file containing the right hand side\n\n";

/*
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h    - vectors
     petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers
*/
//#include <petscmat.h>
#include "petscksp.h" //includes vec,mat,dm,pc also enables ksp
#include<iostream>
#include<string>
#include<chrono>

void prepare_file(KSPType type_ksp, PCType type_pc, char* filename);
PetscErrorCode solve_linear_system(Mat* A, Vec* b, Vec* x, char* matrix_name, PetscMPIInt rank, bool flg_runtime, std::chrono::duration<double, std::milli>* runtime);
PetscErrorCode monitor(KSP,PetscInt,PetscReal,void* ptr);
PetscErrorCode read_matrix_from_file(char* filepath_A, Mat* A);
PetscErrorCode read_vector_from_file(char* filepath_b, Vec* b);
PetscErrorCode create_random_vector(Vec* b, PetscInt size);

PetscErrorCode read_matrix_from_file(char* filepath_A, Mat* A){
    PetscErrorCode ierr = 0;
    PetscViewer view_in_A;
    PetscInt block_size = 1;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filepath_A,FILE_MODE_READ,&view_in_A);CHKERRQ(ierr);
    //Load the matrix; then destroy the viewer.
    ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
    ierr = MatSetType(*A,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatLoad(*A,view_in_A);CHKERRQ(ierr);

    ierr = MatSetBlockSize(*A, block_size);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_in_A);CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode read_vector_from_file(char* filepath_b, Vec* b){
    PetscErrorCode ierr = 0;
    PetscViewer view_in_b;   
    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filepath_b,FILE_MODE_READ,&view_in_b);CHKERRQ(ierr);
    // lead the vecotr, then destroy the viewer
    ierr = VecCreate(PETSC_COMM_WORLD, b);CHKERRQ(ierr);
    ierr = VecSetType(*b,VECMPI);CHKERRQ(ierr);
    ierr = VecLoad(*b,view_in_b);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&view_in_b);CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode create_random_vector(Vec* b, PetscInt size){
    PetscErrorCode ierr = 0;
    
    PetscRandom rand;
    ierr = VecCreate(PETSC_COMM_WORLD, b);CHKERRQ(ierr);
    ierr = VecSetType(*b,VECMPI);CHKERRQ(ierr);
    ierr = VecSetSizes(*b,PETSC_DECIDE,size);CHKERRQ(ierr);
    ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rand);CHKERRQ(ierr);
    
    ierr = VecSetRandom(*b,rand);CHKERRQ(ierr);

    ierr = PetscRandomDestroy(&rand);CHKERRQ(ierr);
    //VecView(*b,PETSC_VIEWER_STDOUT_WORLD);
    return ierr;
}


PetscErrorCode solve_linear_system(Mat* A, Vec* b, Vec* x, char* matrix_name, PetscMPIInt rank, bool flg_runtime, std::chrono::duration<double, std::milli>* runtime){
    PetscErrorCode ierr = 0;

    KSP ksp;
    PC pc;
    // default ksp and pc types. These types can be overwritten with the flag -ksp_type and -pc_type
    KSPType type_ksp = KSPGMRES; // {KSPRICHARDSON, KSPCG, KSPGCR, KSPGMRES, ...}
    PCType type_pc = PCGAMG; // {PCJACOBI,PCBJACOBI,PCILU,PCICC, PCGAMG...} PCILU does not work for MPI
    PetscInt m,n;
    KSPConvergedReason reason;
    PetscReal b_norm;

    // compute norm of RHS vector b
    ierr = VecNorm(*b, NORM_2,&b_norm);CHKERRQ(ierr);


    //#######################################################################################################################################
    //prepare solution vector
    ierr = MatGetSize(*A,&m,&n);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,x);CHKERRQ(ierr);
    ierr = VecSetType(*x,VECMPI);CHKERRQ(ierr);   // use same type as vector b
    ierr = PetscObjectSetName((PetscObject) *x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(*x,PETSC_DECIDE,n);CHKERRQ(ierr);

    //#######################################################################################################################################
    // prepare the KSP solver
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
    ierr = KSPSetType(ksp, type_ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr); // rtol, atol, dtol, maxit
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr); // assign preconditioner to ksp
    ierr = PCSetType(pc,type_pc);CHKERRQ(ierr); // set pc type
    ierr = KSPSetOperators(ksp, *A, *A);CHKERRQ(ierr); // The first matrix defines the linear system. The second matrix is used in constructing the preconditioner
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr); // allows setting flags like -ksp_type to specify the type from the terminal at execution
    
    // the actual types of ksp and pc after KSPSetFromOptions (to store them in the .txt file for later documentation)
    ierr = KSPGetType(ksp, &type_ksp);CHKERRQ(ierr);
    ierr = PCGetType(pc, &type_pc);CHKERRQ(ierr);

    // monitoring (writing residual norm to file) is only done by process 0 
    // and only if runtime is not measures, because monitoring influences the runtime
    if (rank == 0 && !flg_runtime){
        char filename[100];
        sprintf(filename, "rnorm_%s_%s_%s.txt", matrix_name, type_ksp, type_pc);
        prepare_file(type_ksp, type_pc, filename);
        ierr = KSPMonitorSet(ksp, monitor, filename,NULL);CHKERRQ(ierr);    
    }

    //measure the times before and after
    auto t1 = std::chrono::high_resolution_clock::now();

    // solving the linear system
    ierr = KSPSolve(ksp,*b,*x);CHKERRQ(ierr);

    //measure the times before and after
    auto t2 = std::chrono::high_resolution_clock::now();
    
    //get reason for termination of the solver (for meaning of this variable see https://petsc.org/main/docs/manualpages/KSP/KSPConvergedReason/)
    KSPGetConvergedReason(ksp, &reason);

    // rank 0 should either compute the runtime or print the termination reason
    if(rank == 0){
        if (flg_runtime){
            *runtime += t2-t1;
        }else{
            std::cout << "Termination Reason: " << reason << std::endl;
        }
    }
    //#######################################################################################################################################    
    // deallocate memory
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

    return ierr;
}


void prepare_file(KSPType type_ksp, PCType type_pc, char* filename){

    FILE* datafile;
    datafile = fopen(filename, "w");
    fprintf(datafile, "%s+%s\n", type_ksp, type_pc);
    fclose(datafile);

}


PetscErrorCode monitor(KSP kps,PetscInt it,PetscReal rnorm,void* ptr){

    PetscErrorCode ierr = 0;
    char* filename = (char*)ptr;
    FILE* datafile;
    datafile = fopen(filename, "a");
    std::cout << it << " " << rnorm << std::endl;
    fprintf(datafile, "%i %.12e\n", it, rnorm);
    fclose(datafile);
    return ierr;
}


int main(int argc,char **args)
{
    // variable definition
    PetscErrorCode ierr = 0;
    char file_A[PETSC_MAX_PATH_LEN], file_b[PETSC_MAX_PATH_LEN];           /* input file name */
    char filepath_A[PETSC_MAX_PATH_LEN] = "../vectors_and_matrices/";     
    char filepath_b[PETSC_MAX_PATH_LEN] = "../vectors_and_matrices/";     
    PetscBool flg_A, flg_b, flg_runtime;
    Mat A;
    Vec b,x;
    PetscMPIInt rank;
    std::chrono::duration<double, std::milli> runtime;

    // Initialization of petsc
    ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
    
    // read the strings which are given via the flags -a, -b
    // in case of -time, no string is read, just a boolean value "flg_runtime" is set to true
    ierr = PetscOptionsGetString(NULL,NULL,"-a",file_A,sizeof(file_A),&flg_A);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-b",file_b,sizeof(file_b),&flg_b);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-time",NULL,0,&flg_runtime);CHKERRQ(ierr);

    // concatenate file directory with filename
    sprintf(filepath_A, "%s%s", filepath_A, file_A );
    sprintf(filepath_b, "%s%s", filepath_b, file_b );
    
    if (flg_A){
        // read matrix from file if the filename is given via flag -a
        ierr = read_matrix_from_file(filepath_A, &A);CHKERRQ(ierr);
    }
    if(flg_b){
        // (optional flag) read RHS vector from file if the filename is given via flag -b 
        ierr = read_vector_from_file(filepath_b, &b);CHKERRQ(ierr);
    }else{
        // if -b is not given, create a RHS vector with random values between 0 and 1
        PetscInt m,n;
        ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
        //std::cout << m << " " << n << std::endl; 
        ierr = create_random_vector(&b, n);
    }

    // assign rank numbers to processes
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

    // either runtime is measured over 5 runs or the error is monitored (for only one run)
    if(flg_runtime){
        int no_runs = 5;
        int i;
        for(i=0; i<no_runs;i++){ierr = solve_linear_system(&A, &b, &x, file_A, rank, flg_runtime, &runtime);CHKERRQ(ierr);}
        if(rank == 0){std::cout << "Average Runtime: " << runtime.count()/no_runs << std::endl;}
    }else{
        ierr = solve_linear_system(&A, &b, &x, file_A, rank, flg_runtime, &runtime);CHKERRQ(ierr);
    }
    
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);
    //std::cout << type_ksp << std::endl;

    return 0;
}