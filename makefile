#!/bin/bash

#tested with g++-8
GCC     = g++ -Wall -mprefer-vector-width=512 -O3 -ffast-math -mfma -fno-trapping-math -march=native -mtune=native -fopenmp-simd -mavx
#ICC     = icpc -Ofast -qopenmp-simd -fma -xhost -qopt-report=5 -qopt-report-phase=vec -qopt-zmm-usage=high -fPIC -qopenmp
#DEBUG   = -g -fopt-info-optall -fverbose-asm
INCLUDE = -lm -lfftw3 -lblas

testinterface : testinterface.o 
	$(GCC) $(DEBUG) testinterface.o -o testinterface $(INCLUDE)

testinterface.o: ./test/testinterface.cpp ./include/*.hpp ./include/*/*.hpp 
	$(GCC) $(DEBUG) $(INCLUDE) -c ./test/testinterface.cpp -o testinterface.o
