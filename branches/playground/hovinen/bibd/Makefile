LINBOX_INCLUDE_DIR=$(HOME)/build/linbox/include
GMP_DIR=$(HOME)/work/linbox/build

CFLAGS=-g -Wall -I$(LINBOX_INCLUDE_DIR) -I$(GMP_DIR)/include
LIBS=-L$(GMP_DIR)/lib -lgmp

TARGETS=bibd-solve pretty-print stage1 test1

all: $(TARGETS)

test1: test1.o matrix.o sparse-compact-vector.o sparse-vector.o dense-vector.o domain.o
	g++ -g -o $@ $(LIBS) $^
stage1: stage1.o matrix.o sparse-vector.o sparse-compact-vector.o dense-vector.o domain.o
	g++ -g -o $@ $(LIBS) $^
pretty-print: pretty-print.o matrix.o sparse-vector.o sparse-compact-vector.o domain.o
	g++ -g -o $@ $(LIBS) $^
bibd-solve: bibd-solve.o bibd-matrix.o moore-penrose.o matrix.o sparse-compact-vector.o sparse-vector.o dense-vector.o domain.o
	g++ -g -o $@ $(LIBS) $^

dense-vector.o: dense-vector.cpp dense-vector.h vector.h domain.h
sparse-vector.o: sparse-vector.cpp sparse-vector.h vector.h domain.h
sparse-compact-vector.o: sparse-compact-vector.cpp sparse-compact-vector.h vector.h domain.h
matrix.o: matrix.cpp matrix.h vector.h domain.h
domain.o: domain.cpp domain.h
moore-penrose.o: moore-penrose.cpp matrix.h vector.h domain.h
bibd-matrix.o: bibd-matrix.cpp matrix.h sparse-compact-vector.h vector.h domain.h

test1.o: test1.cpp matrix.h sparse-compact-vector.h sparse-vector.h dense-vector.h vector.h domain.h
stage1.o: stage1.cpp matrix.h sparse-vector.h dense-vector.h vector.h domain.h
pretty-print.o: pretty-print.cpp matrix.h sparse-vector.h vector.h domain.h
bibd-solve.o: bibd-solve.cpp bibd-matrix.h moore-penrose.h matrix.h

clean:
	rm -f *.o *~ $(TARGETS)

%.o: %.cpp
	g++ -c $(CFLAGS) $< -o $@
