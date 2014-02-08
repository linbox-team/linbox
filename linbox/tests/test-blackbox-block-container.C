
/* tests/test-blackbox-block-container.C
 *
 * Written by bds
 *
 * -----------------------------------------------------
 *
 * Copyright (c) LinBox
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========


 */


/*! @file  tests/test-blackbox-block-container.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */
#if 0


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/algorithms/blackbox-block-container.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"

#include "test-common.h"

using namespace LinBox;


int main (int argc, char **argv)
{

	bool pass = true;

	static size_t n = 40;
	static size_t k = 4;
	static integer q = 65519U;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'k', "-k K", "Set dimension of blocks matrices to KxK.", TYPE_INT,     &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("block container test", "bbbc");
	ostream& report = commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "over Modular<double>" << std::endl;
	typedef Modular<double> Field;
	typedef BlasMatrix<Field> Block;
	typedef BlasMatrix<Field> Blackbox;
	Modular<double> F (q);

	Blackbox A(F, n, n);
	Block U(F, k, n);
	Block V(F, n, k);
	Block W(F, k, k);

	// Solved a segfault in init.
	BlackboxBlockContainer<Field, Blackbox> B(&A, F, U, V);

	// A more thorough test should be constructed.
#endif

 /*
 * Written by George Yuhasz <yuhasz@gmail.com>
 *
 * --------------------------------------------------------
 *
 */
#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;
// using namespace std;

template<class Blackbox>
bool testContainer (const Blackbox& A, size_t r, size_t c);

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 4;
	static size_t r = 2;
	static size_t c = 2;
	static size_t q = 5;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'r', "-r R", "Set rowdim of blocks R.", TYPE_INT,     &r },
		{ 'c', "-c C", "Set coldim of blocks to C.", TYPE_INT,     &c },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &q },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;

	parseArguments (argc, argv, args);
	Field F ((uint32_t)q);

	commentator().start("bmseq test suite", "BlasMatrix");

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	// ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	commentator().start("SparseMatrix test");
	SparseMatrix<Field> A(F, n, n);
	/*
	if (n > 2) {
		A.setEntry(0,1,F.one);
		A.setEntry(1,2,F.one);
		A.setEntry(2,0,F.one);
	}
	for(size_t i=(n > 2 ? 3 : 0); i<n;i++)
	*/
	for(size_t i=0; i<n;i++)
			A.setEntry(i,n-1-i,F.one);
 	pass = pass and	testContainer(A, r, c);
	commentator().stop("SparseMatrix test");

#if 0 // BlackboxBlockContainer<BlasMatrix<..> > is not working.
	commentator().start("BlasMatrix<Modular<int> > test");
	BlasMatrix<Field> B(F, n, n);
	for(size_t i=0; i<n;i++)
			B.setEntry(i,i,F.one);
	 	pass = pass and testContainer(B, r, c);
	commentator().stop("BlasMatrix<Modular<int> > test");

	commentator().start("BlasMatrix<Modular<double> > test");
	Modular<double> G(q);
	BlasMatrix<Modular<double> > C(G, n, n);
	for(size_t i=0; i<n;i++)
			C.setEntry(i,i,G.one);
	 	pass = pass and testContainer(C, r, c);
	commentator().stop("BlasMatrix<Modular<double> > test");
#endif

	// A more thorough test should be constructed.
	if (pass) commentator().stop("block container test pass");
	else commentator().stop("block container test FAIL");
	return pass ? 0 : -1;
}

template<class Blackbox>
bool testContainer (const Blackbox& A, size_t r, size_t c) {
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	bool pass = true;
	typedef typename Blackbox::Field Field;
	MatrixDomain<Field> MD(A.field());
	size_t n = A.rowdim(); // = A.coldim()
	BlasMatrix<Field> U(A.field(),r,n);
	BlasMatrix<Field> V(A.field(),n,c);
	BlasMatrix<Field> AV(A.field(),n,c);
	BlasMatrix<Field> UAV(A.field(),r,c);
	typename Field::RandIter rand(A.field());
	/*
	for(size_t i=0; i<r;i++)
		for(size_t j=0; j<n; j++)
			U.setEntry(i,j, A.field().zero);
	for(size_t i=0; i<r;i++)
		U.setEntry(i,i, A.field().one);
	*/
	for(size_t i=0; i<r;i++)
		for(size_t j=0; j<n; j++)
			rand.random(U.refEntry(i,j));

	/*
	for(size_t i=0; i<n;i++)
		for(size_t j=0; j<c; j++)
			V.setEntry(i,j, A.field().zero);
	for(size_t i=0; i<c;i++)
		V.setEntry(i,i, A.field().one);
		*/
	for(size_t i=0; i<n;i++)
		for(size_t j=0; j<c; j++){
			rand.random(V.refEntry(i,j));
		}
	MD.copy(AV, V);
	report << std::endl << "A" << std::endl;
	A.write(report);
	report << std::endl << "U" << std::endl;
	U.write(report);
	report << std::endl << "V" << std::endl;
	V.write(report);
	report << std::endl << "AV" << std::endl;
	AV.write(report);
	BlackboxBlockContainer<Field, Blackbox > blockseq(&A,A.field(),U,V);
	MD.mul(UAV,U,AV);
	typename BlackboxBlockContainer<Field, Blackbox >::const_iterator contiter(blockseq.begin());
	report << std::endl << "container size is " << blockseq.size() << std::endl;
	report << std::endl;
	bool pass1 = MD.areEqual(UAV, *contiter);
	if (not pass1) report << "sequences differ at index 0" << std::endl;
	else report << "sequences agree at index 0" << std::endl;
	report << "My UA^0V";
	UAV.write(report) << std::endl;
	report << "Container UA^0V";
	(*contiter).write(report ) << std::endl << std::endl;
	for (size_t i=1; i<10; i++){
		MD.leftMulin(A,AV);
		MD.mul(UAV,U,AV);
		++contiter;
		pass1 = MD.areEqual(UAV, *contiter);
		if (not pass1) report << "sequences differ at index " << i << std::endl;
		else report << "sequences agree at index " << i << std::endl;
		report << "My UA^" << i << "V ";
		UAV.write(report) << std::endl;
		report << "Container UA^" << i << "V ";
		(*contiter).write(report) << std::endl << std::endl;
		pass = pass and pass1;
	}
	return pass;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

