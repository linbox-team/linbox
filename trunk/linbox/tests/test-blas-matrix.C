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

	{ /* Modular<float> */
		//Field
		typedef Modular<float> Field;

		Field F (q);
		commentator().start("Modular<float>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);

		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<float>");
	}

	{ /* ModularBalanced<float> */
		//Field
		typedef ModularBalanced<float> Field;

		Field F (q);
		commentator().start("ModularBalanced<float>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<float>");
	}

	{ /* Modular<double> */
		//Field
		typedef Modular<double> Field;

		Field F (q);
		commentator().start("Modular<double>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);

		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<double>");
	}

	{ /* ModularBalanced<double> */
		//Field
		typedef ModularBalanced<double> Field;

		Field F (q);
		commentator().start("ModularBalanced<double>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<double>");
	}

	{ /* Modular<int64_t> */
#if 0 /*  bug somewhere */
		//Field
		typedef Modular<int64_t> Field;

		Field F (q);
		commentator().start("Modular<int64_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<int64_t>");
#endif
	}

	{ /* ModularBalanced<int64_t> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<int64_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<int64_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<int64_t>");
#endif
	}

	{ /* Modular<uint64_t> */
#if 0 /* not working */
		//Field
		typedef Modular<uint64_t> Field;

		Field F (q);
		commentator().start("Modular<uint64_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<uint64_t>");
#endif
	}

	{ /* ModularBalanced<uint64_t> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<uint64_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<uint64_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<uint64_t>");
#endif
	}

	{ /* Modular<int32_t> */
		//Field
		typedef Modular<int32_t> Field;

		Field F (q);
		commentator().start("Modular<int32_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<int32_t>");
	}

	{ /* ModularBalanced<int32_t> */
		//Field
		typedef ModularBalanced<int32_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<int32_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<int32_t>");
	}

	{ /* Modular<uint32_t> */
#if 0 /*  bug somewhere */
		//Field
		typedef Modular<uint32_t> Field;

		Field F (q);
		commentator().start("Modular<uint32_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<uint32_t>");
#endif
	}

	{ /* ModularBalanced<uint32_t> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<uint32_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<uint32_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<uint32_t>");
#endif
	}

	{ /* Modular<int16_t> */
#if 0 /* not working */
		//Field
		typedef Modular<int16_t> Field;

		Field F (q);
		commentator().start("Modular<int16_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<int16_t>");
#endif
	}

	{ /* ModularBalanced<int16_t> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<int16_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<int16_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<int16_t>");
#endif
	}

	{ /* Modular<uint16_t> */
#if 0 /* not working */
		//Field
		typedef Modular<uint16_t> Field;

		Field F (q);
		commentator().start("Modular<uint16_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"Modular<uint16_t>");
#endif
	}

	{ /* ModularBalanced<uint16_t> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<uint16_t> Field;

		Field F (q);
		commentator().start("ModularBalanced<uint16_t>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"ModularBalanced<uint16_t>");
#endif
	}

	{ /* Modular<char> */
#if 0 /* not working */
		//Field
		typedef Modular<char> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
#endif
	}

	{ /* ModularBalanced<char> */
#if 0 /* not working */
		//Field
		typedef ModularBalanced<char> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
#endif
	}

	{ /* GivaroZpz<uint32_t> */
#if 0 /*  bug somewhere */
		//Field
		typedef GivaroZpz<uint32_t> Field;

		Field F (q);

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
#endif
	}

	{ /* GivaroZpz<integer> */
		//Field
		typedef GivaroZpz<integer> Field;

		Field F (123456789124);
		commentator().start("GivaroZpz<integer>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"GivaroZpz<integer>");
	}

	{ /* PID_integer */
		//Field
		typedef PID_integer Field;

		Field F ;
		commentator().start("PID_integer");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"PID_integer");
	}

	{ /* GivaroExtension<> */
		//Field
		typedef GivaroExtension<> Field;

		Field F(103,4) ;
		commentator().start("GivaroExtension<>");

		typedef 	BlasMatrix<Field,Vector<Field>::Dense>  Matrix ;

		pass = pass && testField<Matrix>(F,m,n, false);
		commentator().stop(MSG_STATUS (pass), (const char *) 0,"GivaroExtension<>");
	}

	commentator().stop(MSG_STATUS(pass),(const char *) 0,"BlasMatrix BB test suite");
	return pass ? 0 : -1;
}


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
