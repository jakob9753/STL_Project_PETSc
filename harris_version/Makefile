all: main

main: main.o
	mpic++ -o main -L/usr/lib/petsc main.o #linking

main.o: main.cpp
	mpic++ -c -Wall -I/usr/include/petsc -std=c++11 -pedantic main.cpp #compiling
