#include <iostream>
#include <vector>
#include <utility>

#include "linbox/blackbox/zero-one.h"
#include "linbox/field/modular.h" 

#include "test-common.h"
#include "test-generic.h"

using std::cout;
using std::endl;
using std::vector;
using LinBox::Modular;
using LinBox::ZeroOne;
using LinBox::uint32;

int main(int argc, char **argv) {
  bool pass = true;
  uint32 prime = 31337;
  size_t *rows, *cols, i;
  static size_t n = 100, iter = 1;

  static Argument args[] = 
    {{ 'n', "-n N", "Set dimension of test matrix to NxN (default 100)", TYPE_INT, &n }, { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 31337)", TYPE_INT, &prime }, { 'i', "-i I", "Perform each test for I iterations (default 1)", TYPE_INT, &iter}};

  parseArguments(argc, argv, args);


  Modular<uint32> afield(prime);

  rows = new size_t[2 * n];
  cols = new size_t[2 * n];

  for(i = 0; i < n; i++) {
    rows[2 * i] = cols[2 * i] = i; // The diagonal
    
    rows[2 * i + 1] = i;         // The reverse diagonal
    cols[2 * i + 1] = 99 - i;
    
  }
  
  ZeroOne<Modular<uint32> > testMatrix(afield, rows, cols, n, n, 2 * n);

  cout << endl << "ZeroOne matrix blackbox test suite" << endl;

  pass = pass && testBlackbox(afield, testMatrix);
  
  delete [] rows;
  delete [] cols;

  return pass ? 0 : -1;
}
