/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 *    test-ntl-sylvester.C
 *    Copyright (C) 2003 Austin Lobo, B. David Saunders
 *    Copyright (C) LinBox
 *
 *    Tests for  Sylvester matrix specification with ntl Arithmetic,
 *    for 2 polynomials in one variable.
 *    LinBox version 2003
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 * Function definitions for block Lanczos iteration
 *-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/

/*! @file  tests/test-ntl-sylvester.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc
 */


#include <iostream>
#include <fstream>
#include <vector>
#include "linbox/field/ntl.h"

#include "linbox/field/ntl.h"
#include "linbox/integer.h"
#include "linbox/blackbox/ntl-sylvester.h"
#include "test-generic.h"



using namespace std;

int main(int argc, char* argv[])
{
	LinBox::commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);
	ostream &report = LinBox::commentator.report(
						     LinBox::Commentator::LEVEL_IMPORTANT,
						     INTERNAL_DESCRIPTION );
	bool pass = true;

	static size_t n = 1000;
	static long q = 134217689;
	//   q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT, &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",
			TYPE_INT, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",
			TYPE_INT, &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);




	//------ Read q and construct F(q)
	NTL::ZZ modulus; 	// prime modulus
	modulus = q;
	size_t  m = n-3;
	// std::cout << std::endl << "Enter a prime number for the modulus of the field: ";
	// std::cin >> modulus;
	report <<  "The modulus is " << modulus << std::endl;
	report <<  "Dimension (m+n) is " << m+n << std::endl;
	NTL::ZZ_p::init(modulus); // NOTE: This is essential for using NTL

	LinBox::commentator.start("Sylvester black box test suite", "Sylvester");
	report <<"Dimension(m+n)= " << m+n << "\t modulus= " << q << endl;

	// typedef LinBox::UnparametricField<NTL::ZZ_p> Field;
	typedef LinBox::NTL_ZZ_p Field;
	typedef Field::Element element;
	typedef std::vector<element> Vector;

	// Now we are using the NTL wrapper as the field, call the instance F
	Field F(q); // XXX same bug ?
	element zero;
	F.init(zero, 0);

	// Use the default constructor to create a matrix
	LinBox::Sylvester<Field> T;

	// Use a special constructor to construct a matrix of dim TSIZE

	Vector pdata(n), qdata(m);

	report << "\n\tpx:=";

	for (size_t i=pdata.size()-1; i > 0; i-- ) {
		pdata[i] = NTL::random_ZZ_p() ;
		report << pdata[i] << "*X^" << i << " + ";
	}
	pdata[0] = NTL::random_ZZ_p() ;
	report << pdata[0];

	report << std::endl;
	report << "\nqx is: \n\t";

	for (size_t i=qdata.size()-1; i > 0; i-- ) {
		qdata[i] = NTL::random_ZZ_p() ;
		report << qdata[i] << "*X^" << i << " + ";
	}

	qdata[0] = NTL::random_ZZ_p() ;
	report << qdata[0];
	report << std::endl;


	LinBox::Sylvester<Field> TT(F,pdata,qdata);
	report << "The matrix is: " << std::endl;
	//  TT.printcp( "cpout.txt");
	//  TT.print(report);

	report << std::endl;


	// Create an interesting input vector called idata
	Vector idata( TT.sysdim() ), odata( TT.sysdim() );
	report << "A random col vector:\npx:=[" << std::endl;

	for (unsigned int i=0; i < idata.size(); i++) {
		idata[i] = NTL::random_ZZ_p() ;

		if (i!= idata.size()-1)     report << idata[i] << ",";
	}
	report << "]\n";

	TT.apply(odata, idata);
	report << "\n\nTesting  apply :--------------------- \nResult is[";

	for (unsigned int i = 0; i < odata.size(); i++)
		report << odata[i] << " ";

	report << "]\n";

	// Apply the matrix to the vector just created
	// Testing the apply function when both input and output are over ZZ_p

	report << "Testing apply Transpose:----------------- \nResult is[";

	TT.applyTranspose(odata, idata);

	for (unsigned int i = 0; i < odata.size(); i++)
		report << odata[i] << " ";

	pass = testBlackbox(TT);
	report <<"<====\tDone Sylvester matrix black box test suite" << endl;


	LinBox::commentator.stop("Sylvester black box test suite");
	return pass ? 0 : -1;

}

