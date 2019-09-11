/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by :
 *          Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *
 * --------------------------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#include "linbox/linbox-config.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>

#include "linbox/integer.h"
#include "givaro/zring.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/dense-matrix.h"

#include "test-common.h"
#include "test-blackbox.h"

using namespace LinBox ;

//! @bug remove NoRW
template<class Matrix>
bool testMatrix(const typename Matrix::Field & F, size_t m, size_t n, bool rw = true) {
	bool pass = true;
	Matrix A(F, m, n);
	A.random();
	pass = pass && testBlackbox(A,rw);

    if (std::min(m,n)>1) {
        BlasSubmatrix<Matrix> B(A,1,1,m/2,n/2);
        pass = pass && testBlackboxNoRW(B);
	}

    BlasSubmatrix<Matrix> C(A,1,1, m-1, n-1);
	if (std::min(m,n)>1) {
		BlasSubmatrix<BlasSubmatrix<Matrix> > D(C,1,1,m/2,n/2);
		pass = pass && testBlackboxNoRW(D);
	}

	return pass ;
}

int main (int argc, char **argv)
{
	// ofstream report;

	bool pass = true;

	static size_t m = 4;
	static size_t n = 10;
	// static size_t nnz = 0;
	static integer q = 101 ;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrix to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension of test matrix to N.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));

	commentator().start("BlasMatrix black box test suite", "triplesbb");

	{ /* Givaro::Modular<float> */
		//Field
		typedef Givaro::Modular<float> Field;

		Field F (q);
		commentator().start("Givaro::Modular<float>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);

		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<float>");
	}

	{ /* Givaro::ModularBalanced<float> */
		//Field
		typedef Givaro::ModularBalanced<float> Field;

		Field F (q);
		commentator().start("Givaro::ModularBalanced<float>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ModularBalanced<float>");
	}

	{ /* Givaro::Modular<double> */
		//Field
		typedef Givaro::Modular<double> Field;

		Field F (q);
		commentator().start("Givaro::Modular<double>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);

		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<double>");
	}

	{ /* Givaro::ModularBalanced<double> */
		//Field
		typedef Givaro::ModularBalanced<double> Field;

		Field F (q);
		commentator().start("Givaro::ModularBalanced<double>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ModularBalanced<double>");
	}

	{ /* Givaro::Modular<int64_t> */
		//Field
		typedef Givaro::Modular<int64_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<int64_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<int64_t>");
	}

	{ /* Givaro::ModularBalanced<int64_t> */
		//Field
		typedef Givaro::ModularBalanced<int64_t> Field;

		Field F (q);
		commentator().start("Givaro::ModularBalanced<int64_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ModularBalanced<int64_t>");
	}
	{ /* Givaro::Modular<uint64_t> */
#if 0
		//Field
		typedef Givaro::Modular<uint64_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<uint64_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<uint64_t>");
#endif
	}
	{ /* Givaro::Modular<int32_t> */
		//Field
		typedef Givaro::Modular<int32_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<int32_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<int32_t>");
	}

	{ /* Givaro::ModularBalanced<int32_t> */
		//Field
		typedef Givaro::ModularBalanced<int32_t> Field;

		Field F (q);
		commentator().start("Givaro::ModularBalanced<int32_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ModularBalanced<int32_t>");
	}

	{ /* Givaro::Modular<uint32_t> */
		//Field
		typedef Givaro::Modular<uint32_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<uint32_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<uint32_t>");
	}

	{ /* Givaro::Modular<int16_t> */
		//Field
		typedef Givaro::Modular<int16_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<int16_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<int16_t>");
	}

	{ /* Givaro::ModularBalanced<int16_t> */
#if 0 /* not working -> PG : this is normal it does not exist anymore in Givaro */
		//Field
		typedef Givaro::ModularBalanced<int16_t> Field;

		Field F (q);
		commentator().start("Givaro::ModularBalanced<int16_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ModularBalanced<int16_t>");
#endif
	}

	{ /* Givaro::Modular<uint16_t> */
		//Field
		typedef Givaro::Modular<uint16_t> Field;

		Field F (q);
		commentator().start("Givaro::Modular<uint16_t>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<uint16_t>");
	}

	{ /* Givaro::Modular<uint32_t> */
		//Field
		typedef Givaro::Modular<uint32_t> Field;

		Field F (q);

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
	}

	{ /* Givaro::Modular<integer> */
		//Field
		typedef Givaro::Modular<integer> Field;

		Field F (integer("123456789124"));
		commentator().start("Givaro::Modular<integer>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Modular<integer>");
	}

	{ /* Givaro::ZRing<Integer> */
		//Field
		typedef Givaro::ZRing<Integer> Field;

		Field F ;
		commentator().start("Givaro::ZRing<Integer>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::ZRing<Integer>");
	}

	{ /* Givaro::Extension<> */
		//Field
		typedef Givaro::Extension<> Field;

		Field F(103,4) ;
		commentator().start("Givaro::Extension<>");

		typedef	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testMatrix<Matrix>(F,m,n, false);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Givaro::Extension<>");
	}

	commentator().stop(MSG_STATUS(pass),(const char *) 0,"BlasMatrix BB test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
