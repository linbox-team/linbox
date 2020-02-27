/* tests/test-rank-ex.C
 * Time-stamp: <18 Dec 19 13:32:34 Jean-Guillaume.Dumas@imag.fr>
 * -----------------------------------------------------
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


/*! @file  tests/test-rank-ex.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */


#include <linbox/linbox-config.h>
#include <iostream>
#include <sstream>
#include <givaro/givrational.h>

#include <linbox/ring/modular.h>
#include <linbox/field/gf2.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/zero-one.h>
#include <linbox/solutions/rank.h>
#include <linbox/util/matrix-stream.h>

#define SP_STOR SparseMatrixFormat::SparseSeq

using namespace LinBox;

bool writing=false;


/// rank or rank mod p
bool runrank (const char * name, size_t irank, int32_t mod)
{

	std::ifstream input (name);
	if (!input) {
		std::cerr << "Error opening matrix file: " << name << std::endl;
		return false;
	}

	size_t r;

	Givaro::QField<Givaro::Rational> QQ;
	LinBox::Timer tim ; tim.clear() ; tim.start();
	MatrixStream<Givaro::QField<Givaro::Rational>> ms (QQ, input);
	SparseMatrix<Givaro::QField<Givaro::Rational>, SP_STOR> A ( ms );

	tim.stop();
	if (writing) std::clog << "matrix is " << A.rowdim() << " by " << A.coldim() << " (" << tim << ")" << std::endl;
	tim.clear() ; tim.start();
	
	if (mod < 0) { // rank over the rational numbers.
		/* We could pick a random prime and work mod that prime, But
		 * the point here is that the rank function in solutions/
		 * handles that issue.  Our matrix here is an integer or
		 * rational matrix and our concept is that we are getting the
		 * rank of that matrix by some blackbox magic inside linbox.
		 */
		LinBox::rank (r, A);
	} else { // rank mod a prime
		typedef Givaro::Modular<double> Field;

		Field F(mod);
		SparseMatrix<Field, SP_STOR > B (F, A.rowdim(), A.coldim());// modular image of A
		MatrixHom::map(B, A);
        tim.stop();
		if (writing) std::clog << "matrix is " << B.rowdim() << " by " << B.coldim() <<" (time for map: "<< tim << ")" << std::endl;

        tim.clear();tim.start();
		//if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::clog) << std::endl;

		// Using the adaptive LinBox Solution
		LinBox::rank(r,B);
	}
	tim.stop();

    bool pass = (r == irank);

    if (pass) {
        std::clog << "[PASS]: " << name << " rank " << r;
        if (mod >0) std::clog << " mod " << mod;
        std::clog << std::endl;
    } else {
        std::cerr << "[FAILED]: Rank found " << std::flush;
        std::cerr << r << std::flush;
        if (mod >0) std::cerr << " mod " << mod;
        std::cerr << " where " << irank << " expected.";
        std::cerr << " (" << tim << " )" << std::endl;
    }
    
    return pass;
}


int main (int argc, char **argv)
{

        // text is written to clog/cerr/cout iff a command line argument is present.
    if (argc > 1) writing = true;

    bool pass(true);

    pass &= runrank("data/rk9_7_10.sms", 10, -1);
    pass &= runrank("data/rk9_7_10.sms", 9, 7);

    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
