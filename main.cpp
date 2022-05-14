
static char help[] = "Reads a PETSc matrix and vector from a file.\n\
  -f0 <input_file> : first file to load (small system)\n\
  -f1 <input_file> : second file to load (larger system)\n\n";

/*
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h    - vectors
     petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers
*/
//#include <petscmat.h>
#include "petscksp.h" //includes vec,mat,dm,pc also enables ksp

int main(int argc,char **args)
{
    Mat A;                      /* matrix */
    Vec b, x, r; 
    PetscViewer view_in_A, view_in_b;          /* viewer */
    char file_A[PETSC_MAX_PATH_LEN];           /* input file name */
    char file_b[PETSC_MAX_PATH_LEN];
    char filepath_A[PETSC_MAX_PATH_LEN] = "vectors_and_matrices/";     
    char filepath_b[PETSC_MAX_PATH_LEN] = "vectors_and_matrices/";     
    PetscBool flg_A;
    PetscBool flg_b;
    KSP ksp;
    PC pc;
    MatType type_A = MATSEQAIJ; // // MATSEQAIJ,MATAIJ,MATMPIAIJ
    VecType type_b = VECSEQ; // VECMPI,VECSEQ,VECSTANDARD
    KSPType type_ksp = KSPCG; // {KSPRICHARDSON, KSPCG, KSPGCR, KSPGMRES, ...}
    PCType type_pc=PCILU; // {PCJACOBI,PCBJACOBI,PCILU,PCICC...}
    PetscInt m,n;
    //PetscDraw image;

    // Initialization
    PetscCall(PetscInitialize(&argc,&args,(char*)0,help));
    

    //#######################################################################################################################################
    //Determine files from which we read the matrix.

    PetscCall(PetscOptionsGetString(NULL,NULL,"-a",file_A,sizeof(file_A),&flg_A));
    PetscCall(PetscOptionsGetString(NULL,NULL,"-b",file_b,sizeof(file_b),&flg_b));
    // Error message if file is not passed with flag -a
    PetscCheck(flg_A,PETSC_COMM_WORLD,PETSC_ERR_USER,"Must indicate binary file with the -a option");
    // Error message if file is not passed with flag -b
    PetscCheck(flg_b,PETSC_COMM_WORLD,PETSC_ERR_USER,"Must indicate binary file with the -b option");
    
    strncat(filepath_A, file_A, PETSC_MAX_PATH_LEN);
    strncat(filepath_b, file_b, PETSC_MAX_PATH_LEN);
    


    //#######################################################################################################################################
    //Open binary file.  Note that we use FILE_MODE_READ to indicate reading from this file.
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filepath_A,FILE_MODE_READ,&view_in_A));


    //Load the matrix; then destroy the viewer.
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetType(A,type_A));
    PetscCall(MatLoad(A,view_in_A));
    PetscCall(PetscViewerDestroy(&view_in_A));
    //PetscCall(MatView(A,PETSC_VIEWER_STDOUT_WORLD));

    //#######################################################################################################################################
    //Open binary file.  Note that we use FILE_MODE_READ to indicate reading from this file.
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD,filepath_b,FILE_MODE_READ,&view_in_b));

    // lead the vecotr, then destroy the viewer
    PetscCall(VecCreate(PETSC_COMM_WORLD, &b));
    PetscCall(VecSetType(b,type_b));
    PetscCall(VecLoad(b,view_in_b));
    PetscCall(PetscViewerDestroy(&view_in_b));

    //#######################################################################################################################################
    //get the dimensions of matrix A
    PetscCall(MatGetSize(A,&m,&n));
    // prepare solution vector
    PetscCall(VecCreate(PETSC_COMM_WORLD,&x));
    PetscCall(VecSetType(x,type_b));    // use same type as vector b
    PetscCall(PetscObjectSetName((PetscObject) x, "Solution"));
    PetscCall(VecSetSizes(x,PETSC_DECIDE,n));

    // initialize the residual vector same as vector x
    PetscCall(VecDuplicate(x,&r));

    //#######################################################################################################################################
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
    PetscCall(KSPSetType(ksp, type_ksp));
    PetscCall(KSPSetTolerances(ksp, 1e-15,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT)); // rtol, atol, dtol, maxit
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPGetPC(ksp,&pc)); // assign preconditioner to ksp
    PetscCall(PCSetType(pc,type_pc)); // set pc type

    PetscCall(KSPSetOperators(ksp, A, A)); // The first matrix defines the linear system. The second matrix is used in constructing the preconditioner
    PetscCall(KSPSolve(ksp,b,x)); 

    //PetscCall(MatResidual(A,b,x,r));
    //PetscCall(VecView(r, PETSC_VIEWER_STDOUT_WORLD));
    //#######################################################################################################################################    
    // deallocate memory
    PetscCall(MatDestroy(&A));
    PetscCall(VecDestroy(&b));
    PetscCall(VecDestroy(&x));
    PetscCall(KSPDestroy(&ksp));
    // finalize petsc
    PetscCall(PetscFinalize());
    return 0;
}