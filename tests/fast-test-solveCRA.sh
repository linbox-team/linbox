#!/bin/bash

mpicxx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -mavx2 -mfma -DLinBoxTestOnly -O3 -Wall -DNDEBUG -U_LB_DEBUG -DNDEBUG -U_LB_DEBUG -I.. -fopenmp -fabi-version=6 -I/home/breusta/build/fflas-ffpack/include -I/home/breusta/build/open-blas/include -I/home/breusta/build/givaro/include        -DLinBoxTestOnly -O3 -Wall -DNDEBUG -U_LB_DEBUG -DNDEBUG -U_LB_DEBUG -I.. -fopenmp -fabi-version=6 -I/home/breusta/build/fflas-ffpack/include -I/home/breusta/build/open-blas/include -I/home/breusta/build/givaro/include      test-solveCRA.C -o test-solveCRA -L/home/breusta/build/open-blas/lib -L/home/breusta/build/givaro/lib -fopenmp -lopenblas -lgivaro -lgmp -lgmpxx

mpirun -n 3 ./test-solveCRA -n 100 -b 10 -q 1 -t 2
