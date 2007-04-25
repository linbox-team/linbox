#include <fstream>
#include <iostream>
#include <vector>
#include <utility>

#include "linbox/blackbox/zo.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/field/modular.h" 
#include "linbox/util/commentator.h"


#include "test-common.h"
#include "test-generic.h"

int main(int argc, char **argv) {
  bool pass = true;
  uint32 prime = 31337;
  static size_t n = 100, iter = 5;

  static Argument args[] = 
    {{ 'n', "-n N", "Set dimension of test matrix to NxN (default 100)", TYPE_INT, &n }, 
     { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 31337)", TYPE_INT, &prime }, 
     { 'i', "-i I", "Perform each test for I iterations (default 5)", TYPE_INT, &iter}};

  parseArguments(argc, argv, args);

  typedef LinBox::Modular<uint32> Field;
  //typedef LinBox::Modular<LinBox::uint32> Field;
  typedef LinBox::ZeroOne<Field> Matrix;

  Field afield(prime);

#if 1
  size_t *rows, *cols, i;
  rows = new size_t[3 * n + 1 - 2];
  cols = new size_t[3 * n + 1 - 2];

  // "arrow" matrix
  for(i = 0; i < n; i++) { rows[i] = 0; cols[i] = i; } // first row
  for(i = 0; i < n - 1; i++) 
    { rows[n+2*i] = i + 1; cols[n+2*i] = 0; rows[n+2*i+1] = i + 1; cols[n+2*i+1] = i + 1; } // first col and the diag
  Matrix testMatrix(afield, rows, cols, n, n, 3 * n - 2);
#else
  Matrix testMatrix(afield);
  //ifstream mat_in("data/m133.b3.200200x200200.sms");
  ifstream mat_in("data/n4c6.b9.186558x198895.sms");
  //ifstream mat_in("data/small21x21.sms");
  testMatrix.read(mat_in);
  //LinBox::Transpose<Matrix> testMat(testMatrix);
#endif

  //print out the dimensions and the number of non-zero entries of the matrix
  std::cout << testMatrix.rowdim() << " " << testMatrix.coldim() << " " << testMatrix.nnz() << std::endl;


  std::cout << std::endl << "ZeroOne matrix blackbox test suite" << std::endl;

  pass = pass && testBlackbox(afield, testMatrix);
  //bool pass2 = testBlackbox(afield, testMat);
  
  //delete [] rows;
  //delete [] cols;

  //return pass&&pass2 ? 0 : -1;
  return pass ? 0 : -1;
}
