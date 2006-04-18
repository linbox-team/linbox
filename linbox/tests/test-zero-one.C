#include <iostream>
#include <vector>
#include <utility>

#include "linbox/blackbox/zero-one.h"
#include "linbox/field/modular.h" 

#include "test-common.h"
#include "test-generic.h"

int main(int argc, char **argv) {
  bool pass = true;
  uint32 prime = 31337;
  size_t *rows, *cols, i;
  static size_t n = 1000, iter = 1;

  static Argument args[] = 
    {{ 'n', "-n N", "Set dimension of test matrix to NxN (default 100)", TYPE_INT, &n }, { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 31337)", TYPE_INT, &prime }, { 'i', "-i I", "Perform each test for I iterations (default 1)", TYPE_INT, &iter}};

  parseArguments(argc, argv, args);

  typedef LinBox::ZeroOne<LinBox::Modular<LinBox::uint32> > Matrix;

  LinBox::Modular<LinBox::uint32> afield(prime);

  rows = new size_t[3 * n];
  cols = new size_t[3 * n];

  // "arrow" matrix
  for(i = 0; i < n; i++) { rows[i] = 0; cols[i] = i; } // first row
  for(i = 0; i < n; i++) { rows[n+i] = i; cols[i] = 0; } // first col
  for(i = 0; i < n; i++) { rows[2*n+i] = i; cols[2*n+i] = i; } // diag
  Matrix testMatrix(afield, rows, cols, n, n, 3*n - 2);

  /*
  for(i = 0; i < n; i++) { rows[i] = i; cols[i] = i; } // diag
  Matrix testMatrix(afield, rows, cols, n, n, n);
  */
/*
  Matrix testMatrix(afield);
  ifstream mat_in("data/n4c6.b9.186558x198895.sms");
  testMatrix.read(mat_in);
  std::cout << testMatrix.rowdim() << " " << testMatrix.coldim() << " " << testMatrix.nnz() << std::endl;
*/


  std::cout << std::endl << "ZeroOne matrix blackbox test suite" << std::endl;

  pass = pass && testBlackbox(afield, testMatrix);
  
  delete [] rows;
  delete [] cols;

  return pass ? 0 : -1;
}
