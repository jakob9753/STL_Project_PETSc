#! /usr/bin/env python3

# How To Use: 
# this file must be in the same directory as petsc_conf.py and PetscBinaryIO.py 
# additinoally copy your matrixmarket (.mtx) file into the same directory, that you want to convert
# call this file from the shell with: python convert_MatrixMarket_to_PetSc.py <yourMatrixMarketFile.mtx> <type>  
# where <type> can be either 'matrix' or 'vector', accordingly to what you want to convert

import os
import sys

os.environ['PETSC_DIR'] = '/usr/lib/petsc'
sys.path.append('/usr/share/petsc/3.14/lib/petsc/bin')

import PetscBinaryIO
import scipy.io

file_name = sys.argv[1]
object_type = sys.argv[2]

matrix_or_vector = scipy.io.mmread(file_name + '.mtx')
with open(file_name + '_petsc', 'w') as petsc_file:
    if object_type == 'matrix':
        PetscBinaryIO.PetscBinaryIO().writeMatSciPy(petsc_file, matrix_or_vector)
    elif object_type == 'vector':
        PetscBinaryIO.PetscBinaryIO().writeVec(petsc_file, matrix_or_vector)