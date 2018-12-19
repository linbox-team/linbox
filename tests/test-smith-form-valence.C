/* Copyright (C) LinBox
 *
 *  Author: JGD
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

/*! @file tests/test-smith-form-valence.h
 * @ingroup tests
 * @brief testing Msith Form computations via the Valence algorithm
 */

#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#include "test-smith-form.h"
#include <linbox/algorithms/smith-form-valence.h>


using namespace LinBox;
typedef Givaro::ZRing<Integer> PIR;
typedef SparseMatrix<PIR>  Blackbox;

static bool testValenceSmith(const char * name) {
    const std::string filename("data/sms.matrix");
	std::ifstream input (filename);
	PIR ZZ;
	MatrixStream< PIR > ms( ZZ, input );
	Blackbox A (ms);
	input.close();
    std::vector<Givaro::Integer> SmithDiagonal;

    PAR_BLOCK {
        smithValence(SmithDiagonal, A, filename);
    }
    
	const size_t k = std::min(A.rowdim(),A.coldim());
	BlasVector<PIR> x(ZZ,k);
    smithForm(x,A);
    
    BlasVector<PIR> sdz(ZZ, SmithDiagonal);

    return checkSNFExample(sdz,x);
}


int main(int argc, char** argv)
{
    bool pass(true);
    
    pass &= testValenceSmith("data/sms.matrix");
    pass &= testValenceSmith("data/fib25.sms");

	return pass;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
