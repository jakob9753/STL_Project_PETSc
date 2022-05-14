#! /usr/bin/env python3

import os
import sys

os.environ['PETSC_DIR'] = '/usr/lib/petsc'
sys.path.append('/usr/share/petsc/3.14/lib/petsc/bin')

import PetscBinaryIO
import scipy.io


matrix_name = sys.argv[1]
matrix = scipy.io.mmread(matrix_name + '.mtx')
with open(matrix_name + '_petsc', 'w') as petsc_file:
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(petsc_file, matrix)