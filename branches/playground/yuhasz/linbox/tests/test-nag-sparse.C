/* 	written by Rich Seagraves <Seagrave@mail.eecis.udel.edu>
	Modified by Z. Wan
*/
#include <linbox/blackbox/nag-sparse.h>
#include <linbox/field/modular.h>
#include <iostream>
#include <fstream>
#include <test-common.h>
#include <test-common.C>

/* Linbox testing routine.  Tests the class to ensure that it functions
 * properly.  Uses the standard diagnostic model used by all Linbox BlackBox's
 * Runs 3 tests, a form test to ensure that results produced are proper, and a
 * function test to ensure that all BlackBoxArchetype API works properly.
 */

template<class Vector,class Field>
bool test(std::ostream & report, const Field &F)
{
	typedef typename Field::Element Element;
	typedef size_t Index;
	Element blank;
	bool res = true;
	int i, j, k;
	cout << "NAGSparse blackbox suite." << endl;
	report << "NAGSparse blackbox suite." << endl;

	cout << "Form Test. . .";
	report << "Form Test. . .";
	report << "A:: Values: 1 3 5  7 9 11 13 15" << endl;
	report << "   columns: 0 0 0  1 2  3  3  3" << endl;
	report << "and   rows: 0 1 4  1 2  0  3  4" << endl;
	report << "A:: mxn = 4x5.  NNZ = 8.  " << endl;
	report << "Now creating Matrix." << endl;

	Element val[8];
	for(k = 0; k < 8; k++)
		val[k] = F.init(blank, 2 * k - 1);

	Index cols[8] = {0,0,0,1,2,3,3,3};
	Index rows[8] = {0,1,4,1,2,0,3,4};
	LinBox::NAGSparse<Field> testNAG(F, val, cols, rows, 4U, 5U, 8U);

	int y_check[4][5];
	Vector x, y;
	report<< "Now checking each row:" << endl;
	report << "Create x as a m element vector w/ a 1 to look at each row" << endl;

	for(i = 0; i < 4; ++i)
		for(j = 0; j < 5; ++j)
			y_check[i][j] = 0;

	F.init(blank);
	x.assign(5,blank);
	y.assign(4,blank); 
	y_check[0][0] = 1;  y_check[0][1] = 3; y_check[0][4] = 4; y_check[1][1] = 7; y_check[2][2] = 9; y_check[3][0] = 11; y_check[3][3] = 13; y_check[3][4] = 15;

	for(k = 0; (unsigned int) k < testNAG.coldim(); ++k) {
		report << "Checking row " << k + 1 << endl;
		if(k > 0) x[k-1] = 0;
		x[k] = 1;
		testNAG.apply(y,x);

		// Now checks y against y_check
		for( i = 0; i < 5; ++i)
			if((int) y[i] != y_check[k][i]) {
				res = false;
				report << "ERROR:  expected " << y_check[k][i] << " at position (" << k << "," << i << "), got "<<y[i] << ".  Not cool." << std::endl;
			}
	}

	/* Now tests Raw and Index Iterators */
	report << "Checking Raw and Index Iterators." << endl;
	report << "Running pre-increment test." << endl;

	// First does a linear run through and checks if the value
	// returned by the RawIterator is the same as that passed into
	// testNAG at initialization.  Runs the checks with 
	// pre-incriments
	typename LinBox::NAGSparse<Field>::RawIterator r_i = testNAG.rawBegin();
	for(k = 0; r_i != testNAG.rawEnd(); ++k, ++r_i) {
		if( *r_i != val[k] ) {
			res = false;
			report << "ERROR: expected element " << k << " of RawIterator to be " << val[k] << ", recieved " << *r_i << "." << endl;
		}
	}
	// Now checks to ensure that all 8 elements were processed
	if(k != 8) {
		res = false;
		report << "ERROR: expected 8 elements to be processed in test, processed " << k + 1 << ".  Something is screwy w/ rawEnd()" << endl;
	}

	report << "Finished with pre-increment test, now running postincrement test." << endl;
	// Now runs the exact same test, but uses post-incrementing
	r_i = testNAG.rawBegin();
	for( k = 0; r_i != testNAG.rawEnd(); k++, r_i++) {
		if( *r_i != val[k] ) {
			res = false;
			report << "ERROR: expected element " << k << " of RawIterator to be " << val[k] << ", recieved " << *r_i << "." << endl;
		}
	}
	if(k != 8) {
		res = false;
		report << "ERROR: expected 8 elements to be processed in test, processed " << k + 1 << ".  Something is screwy w/ rawEnd()" << endl;
	}


	// Now tests the RawIndexIterators using the same style of
	// linear tests as above.  This time checks the index pair
	// returned by the RawIndexIterator against the two row &
	// column initialization arrays used above
	report << "RawIterator test finished.  Now testing RawIndexIterator." << endl;
	typename LinBox::NAGSparse<Field>::RawIndexIterator ri_i = testNAG.indexBegin();
	std::pair<size_t, size_t> testPair;
	report << "First running test using pre-increment." << endl;
	for(k = 0; ri_i != testNAG.indexEnd(); ++k, ++ri_i) {
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
		
	if( k != 8) {
		res = false;
		report << "ERROR: Supposed to process 8 indices.  Only processed " << k + 1 << "." << endl;
	}

	// Now repeats the test using, you guessed it, postincrement
	report << "Now running Index test using postincrement." <<endl;
	for(k = 0; ri_i != testNAG.indexEnd(); k++, ri_i++) {
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
		
	if( k != 8) {
		res = false;
		report << "ERROR: Supposed to process 8 indices.  Only processed " << k + 1 << "." << endl;
	}


	if(res == false) {
		cout << "FAILED." << endl;
		report << "Form Test. . .FAILED." << endl;
	} else {
		cout << "PASSED." << endl;
		report << "Form Test. . .PASSED." << endl;
	}

	// Now on to the functional test
	cout << "FUNCTIONAL TEST. . .";
	report << "FUNCTIONAL TEST. . ." << endl;
	report << "Uses BBGeneralTest(A,F,report)" << endl;

	/*
	if( BBGeneralTest(testNAG, F, report) == false) {
		res = false;
		cout << "FAILED " << endl;
		report << "FUNCTIONAL TEST. . .FAILED" << endl;
	} else {
		cout << "PASSED." << endl;
		report << "FUNCTIONAL TEST. . .PASSED." << endl;
	}
	*/
	return res;
}

template<class Field>
bool test(std::ostream& report, const Field& F) {
	return test<std::vector<typename Field::Element>, Field>(report, F);
}

int main(int argc, char** argv) {
	LinBox::commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);
 	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
  
  	bool pass = true;
  
  
	static size_t n = 10000;
	static long q = 2147483647;
	static int iterations = 1;
	
  static Argument args[] = {
	  { 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n},
	  { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INT, &q},
	  { 'i', "-i I", "Perform each test for I iterations (default 10)",          TYPE_INT,     &iterations },
  };

	 parseArguments (argc, argv, args);
	LinBox::Modular<LinBox::uint32> F(101);
	pass = test(report,F);
	if(pass)
		cout<<"passed\n";
	else
		cout<<"failed\n";
	return 0;
}
