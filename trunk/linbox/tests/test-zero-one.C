#include <iostream>
#include <vector>
#include <fstream>
#include <utility>

#include "linbox/util/commentator.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/field/modular.h" 
#include "linbox/randiter/modular.h"
#include "test-common.h"
#include "test-generic.h"

using std::ostream;
using std::endl;
using std::vector;
using LinBox::Modular;
using LinBox::uint32;
using LinBox::Commentator;
using LinBox::commentator;

template<class Field, class Vector>
static bool testZeroOneBlackBox(Field &F, const Vector &x, size_t n, size_t iter)
{
  typedef typename Field::Element Element;
  Element blank;
  bool res = true;
  int i, j, k, h;
  size_t *rows, *cols;
  Vector y1, y2;
  typename Vector::iterator y1i, y2i;
  LinBox::ModularRandIter<Element> gen(F);

  ostream &report = commentator.report(Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

  report << "ZeroOne blackbox suite." << endl;
  
  report << "Form Test. . .";
  report << "Generating a " << n << "x" << n << " zero-one Matrix with 1s on the main diagonal, and" << endl;
  report << "the reverse diagonal." << endl;
  report << "Now creating Matrix." << endl;
  
  rows = new size_t[2 * n];
  cols = new size_t[2 * n];

  for(i = 0; i < n; i++) {
    rows[2 * i] = cols[2 * i] = i; // The diagonal
    
    rows[2 * i + 1] = i;         // The reverse diagonal
    cols[2 * i + 1] = 99 - i;
    
  }
  
  LinBox::ZeroOne<Field, Vector> testMatrix(F, rows, cols, n, n, 2 * n);

  for(h = 0; h < iter; h++) {
    report << "Now running test iteration: " << h + 1 << endl;

    report << "Testing rowdim: " << endl;
    if( testMatrix.rowdim() != n) {
      report << "rowdim FAIL!" << endl;
      report << "Was expecting " << n << ", recieved " << testMatrix.rowdim() << "." << endl; 
      res = false;
    }
  
    report << "Testing coldim: " << endl;
    if( testMatrix.coldim() != n) {
      report << "coldim FAIL!" <<  endl;
      report << "Was expecting " << n << ", recieved " << testMatrix.coldim() << "." << endl;
      res = false;
    }


    report << "Testing apply:" << endl;
    F.init(blank);
    y1.assign(100, blank); // ensure proper size of y values
    y2.assign(100,blank); 
  
    report << "Performing a manual application of Ax." << endl;
    for(i = 0; i < 100; i++) {
      for(j = 0; j < 100; j++) {
	if( i == j) {
	  F.addin( y1[i], x[j] );
	}
	else if(j == 99 - i) {
	  F.addin( y1[i], x[j]  );
	}
      }
    }
    // Now perform the apply
    report << "Performing application Ax using blackbox." << endl;
    testMatrix.apply(y2, x);
  
    // Now compare the two
    k = 0;
    for( i = 0; i < 100; i++ ) {
      if( y1[i] != y2[i]) {
	report << "Apply test FAIL!" << endl;
	report << "Was expecting " << y1[i] << ", but recieved " << y2[i] << "." << endl;
	report << "At location: " << i << endl;
	res = false;
	k++;
      }
    }
    if(k > 0) report << "Total number of erroneous vector entries: " << k << endl;

    report << "Testing applyTranspose" << endl;
    // 0's out y for the applyTranspose test
    y1.assign(100,blank);
    for(i = 0; i < 100; i++) {
      for(j = 0; j < 100; j++) {
	
	if(i == j) {
	  F.addin(y1[i], x[j]);
	}
      
	else if(j == 99 - i) {
	F.addin(y1[i], x[j]);
	}
      }
    }
  
    // Now perform the apply transpose
    testMatrix.applyTranspose(y2, x);
  
    k = 0;
    // Now go through and compare the two
    for( i = 0; i < 100; i++) {
      if( y1[i] != y2[i] ) {
	report << "applyTranspose test FAIL!" << endl;
	report << "Was expecting " << y1[i] << ", but recived " << y2[i] << " instead." << endl;
      res = false;
      report << "At location: " << i << endl; 
      k++;
      }
    }
    
    if( k > 0) report << "Total numberof erroneous vector entries: " << k << endl;

    report << "Checking Raw and Index Iterators." << endl;
    report << "Running pre-increment test." << endl;
    LinBox::ZeroOne<Field, Vector>::RawIterator ra_i = testMatrix.rawBegin();
    F.init(blank, 1); 
    for(i = 0; ra_i != testMatrix.rawEnd(); ++i, ++ra_i) {
      if(blank != *ra_i) {
	report << "Error, was expecting a 1, but got " << *ra_i << " instead." << endl;
      }
    }
    if( i != 2 * n) {
      report << "ERROR, should have processed " << 2 * n << " elements, processed " << i << " instead." << endl;
    }
  
    report << "Now running post-increment test." << endl;
    ra_i = testMatrix.rawBegin();
    for(i = 0; ra_i != testMatrix.rawEnd(); ++i, ra_i++) {
      if( blank != *ra_i) {
	report << "Error, was expecting a 1, but got " << *ra_i << " instead." << endl;
      }
    }
    if( i != 2 * n) {
      cout << "Error, should have processed " << 2 * n << " elements, processed " << i << " instead." << endl;
    }   

    // Now tests the RawIndexIterators using the same style of
    // linear tests as above.  This time checks the index pair
    // returned by the RawIndexIterator against the two row &
    // column initialization arrays used above
    report << "RawIterator test finished.  Now testing RawIndexIterator." << endl;
    LinBox::ZeroOne<Field, Vector>::RawIndexIterator ri_i = testMatrix.indexBegin();
    pair<size_t, size_t> testPair;
    report << "First running test using pre-increment." << endl;
    for(k = 0; ri_i != testMatrix.indexEnd(); ++k, ++ri_i) {
      testPair = *ri_i;
      if( testPair.first != rows[k] ) {
	res = false;
        report << "ERROR: expected row " << k << " of RawIndexIterator to be " << rows[k] << ", recieved " << testPair.first << "." << endl;
      }
      if( testPair.second != cols[k] ) {
	res = false;
	report << "ERROR: expect cols " << k << " of RawIndexIterator to be " << cols[k] << ", recieved " << testPair.second << "." << endl;
      }
    }
  
    if( k != 2 * n) {
      res = false;
      report << "ERROR: Supposed to process " << 2 * n << " indices.  Only processed " << k + 1 << "." << endl;
    }
    
    // Now repeats the test using, you guessed it, postincrement
    report << "Now running Index test using postincrement." <<endl;
    ri_i = testMatrix.indexBegin();
    for(k = 0; ri_i != testMatrix.indexEnd(); k++, ri_i++) {
      testPair = *ri_i;
      if( testPair.first != rows[k] ) {
	res = false;
	report << "RawIndexIterator test FAIL." << endl;
	report << "Expected row " << k << " of RawIndexIterator to be " << rows[k] << ", recieved " << testPair.first << "." << endl;
      }
      if( testPair.second != cols[k] ) {
	res = false;
	report << "RawIndexIterator test FAIL." << endl;
	report << "Expected col " << k << " of RawIndexIterator to be " << cols[k] << ", recieved " << testPair.second << "." << endl;
      }
    }
  
    if( k != 2 * n) {
      res = false;
      report << "ERROR: Supposed to process " << 2 * n << " indice pairs.  Only processed " << k << "." << endl;
    }
  

  }
  

  if(res == false) {
    report << "Form Test. . .FAILED." << endl;
  }
  else {
    report << "Form Test. . .PASSED." << endl;
  }

  delete [] rows;
  delete [] cols;
 
  return res;
}
  

int main(int argc, char **argv) {
  bool pass = true;
  uint32 prime = 31337;
  vector<uint32> x;
  vector<uint32>::iterator xp;
  size_t i, n = 100, p, iter = 1;
  int pb;

  static Argument args[] = 
    {{ 'n', "-n N", "Set dimension of test matrix to NxN (default 100)", TYPE_INT, &n }, { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 31337)", TYPE_INT, &prime }, { 'i', "-i I", "Perform each test for I iterations (default 1)", TYPE_INT, &iter}};

  parseArguments(argc, argv, args);


  Modular<uint32> afield(prime);
  LinBox::ModularRandIter<uint32> gen(afield);

  x.resize(n);

  for(xp = x.begin(); xp != x.end(); ++xp) {
    gen.random(*xp);
  }

  commentator.getMessageClass( INTERNAL_DESCRIPTION).setMaxDepth(3);
  commentator.getMessageClass( INTERNAL_DESCRIPTION).setMaxDetailLevel( Commentator::LEVEL_UNIMPORTANT);


  if( testZeroOneBlackBox(afield, x, n, iter) ) {
    pb = 0;
  }
  else
    pb = -1;

  return pb;
}
