/* Copyright (C) LinBox
 *
 * Copyright (C) 2002 Austin Lobo, B. David Saunders
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-ntl-toeplitz.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc
 */



#include <linbox/linbox-config.h>
#include <iostream>
#include <fstream>

#include "linbox/ring/ntl.h"

#include "linbox/integer.h"
#include "linbox/vector/blas-vector.h"
//#include "linbox/blackbox/ntl-toeplitz.h"
#include "linbox/blackbox/toeplitz.h"


#include "test-generic.h"


using namespace std;


int main(int argc, char* argv[])
{
	LinBox::commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);
	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	bool pass = true;

	static size_t n = 1000;
	static int64_t q = 65521;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	LinBox::parseArguments (argc, argv, args);


	//------ Read q and construct F(q)
	//NTL::ZZ modulus; 	// prime modulus
	//modulus = q;


	//std::cout << std::endl << "Enter a prime number for the modulus of the field: ";
	//std::cin >> modulus;
	//report <<  "The modulus is " << modulus << std::endl;
	report <<  "The modulus is " << q << std::endl;
	//NTL::ZZ_p::init(modulus); // NOTE: This is essential for using NTL

	LinBox::commentator().start("Toeplitz black box test suite", "Toeplitz");
	report  <<"     \tDimension= " << n << "\t modulus= " << q << endl;

	typedef LinBox::NTL_ZZ_p Field;
	typedef LinBox::NTL_ZZ_pX PolyRing;
	// typedef Field::Element element;
	typedef LinBox::BlasVector<Field> Vector;

	// Now we are using the NTL wrapper as the field, call the instance F
	Field F(q);

	// Use the default constructor to create a matrix
	LinBox::Toeplitz<Field, PolyRing> T(F);

	// Use a special constructor to construct a matrix of dim TSIZE
	size_t TSIZE = 2*(n)-1;
	Vector tdata(F,TSIZE);
	report << "The random vector is:" << std::endl;
	for (unsigned int i=0; i < tdata.size(); i++) {
		tdata[i] = NTL::random_ZZ_p() ;
		report << tdata[i] << " ";
	}
	report << std::endl;

	LinBox::Toeplitz<Field, PolyRing> TT(F,tdata);
	report << "The matrix is: " << std::endl;
	TT.write(report);

	// Create an interesting input vector called idata
	Vector idata(F,(TSIZE+2)/2), odata(F,(TSIZE+2)/2);
	report << "A random col vector:\t" << std::endl;
	for (unsigned int i=0; i < idata.size(); i++) {
		idata[i] = NTL::random_ZZ_p() ;
		report << idata[i] << " ";
	}
	report << std::endl;

	// Apply the matrix to the vector just created
	// Testing the apply function when both input and output are over ZZ_p
	TT.applyTranspose(odata, idata);
	report << "Testing apply Transpose:----------------- \nResult is[";
	for (unsigned int i = 0; i < odata.size(); i++)
		report << odata[i] << " ";
	report << "]\n";



	TT.apply(odata, idata);
	report << "\n\nTesting  apply :--------------------- \nResult is[";
	for (unsigned int i = 0; i < odata.size(); i++)
		report << odata[i] << " ";
	report << "]\n";

	// TT.setToUniModLT();
	//  report << "Setting the matrix to UniGivaro::Modular Upper Triangular";
	//  TT.setToUniModUT();
	//  report << "The Upper Triangular matrix is: " << std::endl;
	//TT.write();

	pass = testBlackboxNoRW(TT);
	LinBox::commentator().stop("Toeplitz black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
