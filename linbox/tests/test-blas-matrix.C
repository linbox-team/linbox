
/*
 * Copyright (C) 2014 the LinBox group
 *
 * Written by :
 *          BB <bbboyer@ncsu.edu>
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
#include "linbox/integer.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/modular.h"
#include "linbox/field/givaro.h"

#include "linbox/matrix/dense-matrix.h"

#include "test-common.h"
#include "test-blackbox.h"

using namespace LinBox ;

//! @bug remove rw
template<class Matrix>
bool testField(const typename Matrix::Field & F, size_t m, size_t n, bool rw = true) {
	bool pass = true;
	Matrix A(F, m, n);
	A.random();
	pass = pass && testBlackbox(A,rw);

	if (std::min(m,n)>1) {
		BlasSubmatrix<Matrix> B(A,1,1,m/2,n/2);
		pass = pass && testBlackboxNoRW(B);
	}

	return pass ;
}

int main (int argc, char **argv)
{
	// ofstream report;

	bool pass = true;

	static size_t m = 4;
	static size_t n = 20;
	// static size_t nnz = 0;
	static integer q = 2147483647U;
	q = 101 ;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrix to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension of test matrix to N.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));

	commentator().start("BlasMatrix black box test suite", "triplesbb");

	{
		//Field
		typedef Modular<double> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
	}

	{
		//Field
		typedef Modular<int64_t> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
	}

#if 0 /*  bug somewhere */
	{
		//Field
		typedef GivaroZpz<Givaro::Unsigned32> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
	}
#endif

	{
		//Field
		typedef GivaroZpz<integer> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
	}

	{
		//Field
		typedef PID_integer Field;

		Field F ;

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
	}

	{
		//Field
		typedef GivaroExtension<> Field;

		Field F(103,4) ;

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n, false);
	}

	commentator().stop(MSG_STATUS(pass),(const char *) 0,"TriplesBB black box test suite");
	return pass ? 0 : -1;
}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
