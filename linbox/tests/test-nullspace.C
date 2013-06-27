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
	typedef typename Field::Element			Element;

	//Commentator commentator;
	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	//commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().start ("Testing NullSpace Decomposition","testNullSpace",(unsigned int)iterations);
	// typename Field::Element one,zero;
	// F.init(one,1UL);
	// F.init(zero,0UL);

	bool ret = true;
	{
		size_t min = std::min(m,n);
		if (rank > min ) {
			rank = min; // rank <= min(m,n)...
		}
	}
	for (int k=0; k<iterations; ++k) {

		commentator().progress(k);
		Element * A = new Element[m*n];
		size_t ld_a =  n ;
		size_t wd_a =  m ;
		RandomMatrixWithRank(F,A,m,n,rank);

		Element * Abis = new Element[m*n]; // copie de A
		for (size_t i=0; i<m*n; ++i)
			*(Abis+i) = *(A+i);
		size_t ker_dim = 0 ; // or coker_dim
		Element * Kern  = NULL;
		size_t ld_k = 0 ;
		if (a_droite) {
			NullSpaceBasis (F, LinBoxTag::Side::Right,m,n,A,ld_a,Kern,ld_k,ker_dim);
			if (ker_dim != (ld_a - rank)) {
				ret = false;
				cout << "faux : (1) mauvaises dim : " << ker_dim << " != " << (ld_a - rank) << endl;
				delete[] Kern;
				delete[] A;
				delete[] Abis;
				break ;
			}
		}
		else {
			NullSpaceBasis (F, LinBoxTag::Side::Left,m,n,A,ld_a,Kern,ld_k,ker_dim);
			if (ker_dim != (wd_a - rank) ) {
				ret = false;
				cout << "faux : (1) mauvaises dim " << ker_dim << " != " << (wd_a - rank)  << endl;
				delete[] Kern;
				delete[] A;
				delete[] Abis;
				break ;
			}
		}
		size_t ld_ker = (a_droite)?ker_dim:m ;
		size_t wd_ker = (a_droite)?n:ker_dim ;
		assert(ld_ker == ld_k) ;
		size_t ld_n = (a_droite)?ker_dim:ld_a;
		size_t wd_n = (a_droite)?wd_a:ker_dim;
		assert(CheckRank(F,Kern,wd_ker,ld_ker,ld_ker,ker_dim)); // ...il est bien de rang plein...
		Element * NullMat = new Element[ld_n*wd_n] ;// ...et on s'attend à ce que ça soit nul !

		if ( a_droite){
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, wd_a, ld_ker, ld_a,
				     F.one, Abis, ld_a, Kern, ld_ker , F.zero, NullMat, ld_n);
		}
		else{
			FFLAS::fgemm(F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, wd_ker,ld_a, ld_ker,
				     F.one,  Kern, ld_ker , Abis, ld_a, F.zero, NullMat, ld_n);
		}

		// write_field (F, std::cout<<"final: NullMat"<<std::endl, NullMat, (int)wd_n, (int)ld_n, (int)ld_n, true);
		//write_field (F, std::cout<<"A="<<endl, A, m, n, n,true);
		//write_field (F, std::cout<<"Abis="<<endl, Abis, m, n, n, true);
		delete[] Abis ;
		delete[] A ;
		delete[] Kern ;
#if 1
		for (size_t i = 0 ; i < wd_n ; ++i ){
			for (size_t j = 0 ; j < ld_n ; ++j ) {
				if (!F.isZero(*(NullMat + j+i*ld_n)) ){
					    	write_field (F, std::cout<<"faux : (3) NullMat pas nulle. "<<std::endl, NullMat, (int)wd_n, (int)ld_n, (int)ld_n, true);
					ret = false;
					break;
				}
			}
			if (!ret)
				break;
		}
#endif
		delete[] NullMat;
		if (!ret)  break;


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

	TESTE("right kernel");
	if (!testNullSpaceBasis (F, m,n,r, iterations, true))
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
