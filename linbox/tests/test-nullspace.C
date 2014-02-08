/* Copyright (C) 2009 LinBox
 * Written by <boyer.brice@gmail.com>
 * Inspired and adapted from test-ffpack.C
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
 *
 * Function definitions for block Lanczos iteration
 */

/** \file tests/test-nullspace.C
 * \brief Tests the dense nullspace functions for Zp
 * @ingroup tests
 * @test dense nullspace
 *  @todo test non dense nullspace
 *  @todo test for submatrices
 *  @todo make sure this is faster than FFPACK ?
 */

#include "linbox/linbox-config.h"
#include <iostream>
#include "linbox/integer.h"
#include "linbox/matrix/matrix-domain.h"
//#include "linbox/field/givaro-zpz.h"
#include "linbox/field/modular.h"
//#include "fflas-ffpack/ffpack/ffpack.h"
#include "linbox/algorithms/dense-nullspace.h"

#include "./test-common.h"
#include "fflas-ffpack/utils/Matio.h"
#include "linbox/linbox-tags.h"

#include "test-matrix-utils.h"

using namespace LinBox;




/*!
 * @brief Tests the NullSpace routines.
 * @param F field
 * @param m row
 * @param n col
 * @param rank \p n-rank is the size of the NullSpace
 * @param iterations number of its
 * @param a_droite \p true if.. \p false if on the left
 * @return \p true hopefully if test's passed!
 */
template <class Field >
static bool testNullSpaceBasis (const Field& F, size_t m, size_t n, size_t rank, int iterations, bool a_droite)
{

	//Commentator commentator;
	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().start ("Testing NullSpace Decomposition","testNullSpace",(unsigned int)iterations);

	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
	}
	for (int k=0; k<iterations; ++k) {

		commentator().progress(k);
		BlasMatrix<Field> A(F,m,n+5);
		BlasSubmatrix<BlasMatrix<Field> > Aref(A,0,0,m,n);
		RandomMatrixWithRank(F,A.getWritePointer(),m,n,n+5,rank);

		//tests on a submatrix (more general/prone to errors)
		BlasSubmatrix<BlasMatrix<Field> > Abis(A,0,0,m,n);
		BlasMatrix<Field> Kern(F);
		size_t ker_dim;
		if (a_droite) {
			NullSpaceBasis (Tag::Side::Right,Abis,Kern,ker_dim);
			if (ker_dim != (Abis.coldim() - rank)) {
				ret = false;
				cout << "wrong: (1) bad dim : " << ker_dim << " != " << (Abis.coldim() - rank) << endl;
				break ;
			}
		}
		else {
			NullSpaceBasis ( Tag::Side::Left,Abis,Kern,ker_dim);
			if (ker_dim != (Abis.rowdim() - rank) ) {
				ret = false;
				cout << "wrong : (1) bad dim " << ker_dim << " != " << (Abis.rowdim() - rank)  << endl;
				break ;
			}
		}

		assert(CheckRank(F,Kern,ker_dim));
		BlasMatrix<Field> NullMat(F,(a_droite)?m:ker_dim,(a_droite)?ker_dim:n);
		BlasMatrixDomain<Field> BMD(F);

		if ( a_droite){
			BMD.mul(NullMat,Aref,Kern);
		}
		else{
			BMD.mul(NullMat,Kern,Aref);
		}

		if (!BMD.isZero(NullMat)) {
			std::cout << "wrong (3) NullMat non zero :" << NullMat << std::endl;
			ret = false;
			break;
		}

	}

	commentator().stop(MSG_STATUS (ret), (const char *) 0, "testNullSpace");
	return ret;
}

int main(int argc, char** argv)
{
	//-----------------------------------------------------------------------
	// Choice of the finite field representation
	//typedef GivaroZpz<Std32> Field;
	typedef Modular<double> Field;
	//typedef Modular<float> Field;
	//typedef Modular<uint32_t> Field;
	//------------------------------------------------------------------------

	bool pass = true;

	static size_t n = 15;
	static size_t m = 8;
	static size_t r = 3;
	static integer q = 101;
	static int iterations =4;

	static Argument args[] = {
		{ 'n', "-n N", "Set width of test matrices.",			TYPE_INT,     &n },
		{ 'm', "-m M", "Set hight of test matrices.",			TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices.",			TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].",		TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",           TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	Field F (q);

	//srand(time (NULL));

	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().start("NullSpace test suite", "nullspace");

	std::ostream& report = commentator().report();


	TESTE("left kernel");
	if (!testNullSpaceBasis (F, m,n,r, iterations, false))
		pass=false;
	RAPPORT("left kernel");

	TESTE("left kernel");
	if (!testNullSpaceBasis (F, n,m,r, iterations, false))
		pass=false;
	RAPPORT("left kernel");

	TESTE("left kernel");
	if (!testNullSpaceBasis (F, m,n,0, iterations, false))
		pass=false;
	RAPPORT("left kernel");


	TESTE("left kernel");
	if (!testNullSpaceBasis (F, n,m,0, iterations, false))
		pass=false;
	RAPPORT("left kernel");

	TESTE("left kernel");
	if (!testNullSpaceBasis (F, m,n,std::min(m,n), iterations, false))
		pass=false;
	RAPPORT("left kernel");





	TESTE("right kernel");
	if (!testNullSpaceBasis (F, m,n,r, iterations, true))
		pass=false;
	RAPPORT("right kernel");

	TESTE("right kernel");
	if (!testNullSpaceBasis (F, n,m,r, iterations, true))
		pass=false;
	RAPPORT("right kernel");

	TESTE("right kernel");
	if (!testNullSpaceBasis (F, m,n,0, iterations, true))
		pass=false;
	RAPPORT("right kernel");

	TESTE("right kernel");
	if (!testNullSpaceBasis (F, n,m,0, iterations, true))
		pass=false;
	RAPPORT("right kernel");

	TESTE("right kernel");
	if (!testNullSpaceBasis (F, n,m,std::min(m,n), iterations, true))
		pass=false;
	RAPPORT("right kernel");



	// if we are here, no RAPPORT exited
	report << "\033[1;32m +++ ALL MY TESTS PASSED +++\033[0;m" << endl;


	commentator().stop(MSG_STATUS (pass),"NullSpace test suite");
	return (pass ? 0 : -1);
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
