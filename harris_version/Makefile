include /usr/lib/petsc/lib/petsc/conf/variables
include /usr/lib/petsc/lib/petsc/conf/rules

all: main

main: main.o
	#g++ -o main main.o -L/usr/lib/petsc #linking
	${CLINKER} -o main main.o ${PETSC_LIB} #linking

main.o: main.cpp
	g++ -c -Wall -std=c++11 -I/usr/include/petsc -pedantic main.cpp #compiling
