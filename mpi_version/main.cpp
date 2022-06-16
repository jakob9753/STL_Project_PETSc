
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


void prepare_file(KSPType type_ksp, PCType type_pc, char* filename);
PetscErrorCode solve_linear_system(Mat* A, Vec* b, Vec* x, char* matrix_name);
PetscErrorCode monitor(KSP,PetscInt,PetscReal,void* ptr);
PetscErrorCode read_matrix_from_file(char* filepath_A, Mat* A);
PetscErrorCode read_vector_from_file(char* filepath_b, Vec* b);
PetscErrorCode create_random_vector(Vec* b, PetscInt size);

PetscErrorCode read_matrix_from_file(char* filepath_A, Mat* A){
    PetscErrorCode ierr = 0;
    PetscViewer view_in_A;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filepath_A,FILE_MODE_READ,&view_in_A);CHKERRQ(ierr);
    //Load the matrix; then destroy the viewer.
    ierr = MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
    ierr = MatSetType(*A,MATMPIAIJ);CHKERRQ(ierr);
    ierr = MatLoad(*A,view_in_A);CHKERRQ(ierr);
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


PetscErrorCode solve_linear_system(Mat* A, Vec* b, Vec* x, char* matrix_name){
    PetscErrorCode ierr = 0;

    KSP ksp;
    PC pc;
    // default ksp and pc types. These types can be overwritten with the flag -ksp_type and -pc_type
    KSPType type_ksp = KSPGMRES; // {KSPRICHARDSON, KSPCG, KSPGCR, KSPGMRES, ...}
    PCType type_pc = PCGAMG; // {PCJACOBI,PCBJACOBI,PCILU,PCICC, PCGAMG...} PCILU does not work for MPI
    PetscInt m,n;
    KSPConvergedReason reason;
    // assign ranks numbers to processes
    PetscMPIInt rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
    // print rank of each process
    ierr = PetscPrintf(PETSC_COMM_SELF,"Hello World from %d\n",rank);CHKERRQ(ierr);
    

    //#######################################################################################################################################
    //prepare solution vector
    ierr = MatGetSize(*A,&m,&n);CHKERRQ(ierr);

    ierr = VecCreate(PETSC_COMM_WORLD,x);CHKERRQ(ierr);
    ierr = VecSetType(*x,VECMPI);CHKERRQ(ierr);   // use same type as vector b
    ierr = PetscObjectSetName((PetscObject) *x, "Solution");CHKERRQ(ierr);
    ierr = VecSetSizes(*x,PETSC_DECIDE,n);CHKERRQ(ierr);

    //#######################################################################################################################################
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
    ierr = KSPSetType(ksp, type_ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr); // rtol, atol, dtol, maxit
    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr); // assign preconditioner to ksp
    ierr = PCSetType(pc,type_pc);CHKERRQ(ierr); // set pc type
    ierr = KSPSetOperators(ksp, *A, *A);CHKERRQ(ierr); // The first matrix defines the linear system. The second matrix is used in constructing the preconditioner
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    ierr = KSPGetType(ksp, &type_ksp);CHKERRQ(ierr);
    ierr = PCGetType(pc, &type_pc);CHKERRQ(ierr);

    // monitoring (writing residual norm to file) is only done by process 0
    if (rank == 0){
        char filename[100];
        sprintf(filename, "rnorm_%s_%s_%s.txt", matrix_name, type_ksp, type_pc);
        prepare_file(type_ksp, type_pc, filename);
        ierr = KSPMonitorSet(ksp, monitor, filename,NULL);CHKERRQ(ierr);    
    }
    
    ierr = KSPSolve(ksp,*b,*x);CHKERRQ(ierr);

    KSPGetConvergedReason(ksp, &reason);
    if(rank == 0){
        std::cout << "Termination Reason: " << reason << std::endl;

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
    PetscErrorCode ierr = 0;
    
    char file_A[PETSC_MAX_PATH_LEN], file_b[PETSC_MAX_PATH_LEN];           /* input file name */
    char filepath_A[PETSC_MAX_PATH_LEN] = "../vectors_and_matrices/";     
    char filepath_b[PETSC_MAX_PATH_LEN] = "../vectors_and_matrices/";     
    PetscBool flg_A, flg_b;
    Mat A;
    Vec b,x;
    // Initialization
    ierr = PetscInitialize(&argc,&args,(char*)0,help);CHKERRQ(ierr);
    
    ierr = PetscOptionsGetString(NULL,NULL,"-a",file_A,sizeof(file_A),&flg_A);CHKERRQ(ierr);
    ierr = PetscOptionsGetString(NULL,NULL,"-b",file_b,sizeof(file_b),&flg_b);CHKERRQ(ierr);
    
    sprintf(filepath_A, "%s%s", filepath_A, file_A );
    sprintf(filepath_b, "%s%s", filepath_b, file_b );
    if (flg_A){
        ierr = read_matrix_from_file(filepath_A, &A);CHKERRQ(ierr);
    }
    if(flg_b){
        ierr = read_vector_from_file(filepath_b, &b);CHKERRQ(ierr);
    }else{
        PetscInt m,n;
        ierr = MatGetSize(A,&m,&n);CHKERRQ(ierr);
        std::cout << m << " " << n << std::endl; 
        ierr = create_random_vector(&b, n);
    }
    

    ierr = solve_linear_system(&A, &b, &x, file_A);CHKERRQ(ierr);


    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);
    //std::cout << type_ksp << std::endl;

    return 0;
}