# convert_MatrixMarket_to_PetSc

This program converts a matrix or a vector in .mtx format to a pets binary file. 


To run it follow these steps:

1) in terminal run: python convert_MatrixMarket_to_PetSc.py <matrixMarketFile> <type>
Explanation:	
		<matrixMarketFile> is the name of the file you want to convert without the extension ".mtx". So if your filename is "matrix.mtx", just type "matrix".
		
		<type> is either "matrix" or "vector", depending on the type of object you want to convert