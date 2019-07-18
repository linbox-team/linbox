/*
 * examples/smith.C
 *
 * Copyright (C) 2005, 2010  D. Saunders, Z. Wang, J-G Dumas
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

/** \file examples/smith.C
 * @example  examples/smith.C
 \brief mod m Smith form by elmination
 \ingroup examples

 \author bds & zw

 Various Smith form algorithms may be used for matrices over the
 integers or over Z_m.  Moduli greater than 2^32 are not supported here.
 Several types of example matrices may be constructed or the matrix be read from a file.
 Run the program with no arguments for a synopsis of the command line parameters.

 For the "adaptive" method, the matrix must be over the integers.
 This is expected to work best for large matrices.

 For the "2local" method, the computation is done mod 2^32.

 For the "local" method, the modulus must be a prime power.

 For the "ilio" method, the modulus may be arbitrary composite.
 If the modulus is a multiple of the integer determinant, the integer Smith form is obtained.  
 Determinant plus ilio may be best for smaller matrices.

 This example was used during the design process of the adaptive algorithm.
 The matrix example generation code that was here is now in matrices.C.
*/

#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;


#include <linbox/ring/modular.h>

#include <linbox/util/timer.h>

#include <linbox/ring/local2_32.h>
// place A: Edit here and at place B for ring change
//#include <linbox/ring/pir-modular-int32.h>
#include <linbox/ring/pir-ntl-zz_p.h>
#include <linbox/algorithms/smith-form-local.h>
#include <linbox/algorithms/smith-form-iliopoulos.h>
#include <linbox/algorithms/smith-form-adaptive.h>

#include <linbox/algorithms/matrix-hom.h>

using namespace LinBox;

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);
template <class PIR> void Mat(DenseMatrix<PIR>& M, string src);

int main(int argc, char* argv[])
{
// place B: Edit here and at place A for ring change
	//typedef PIRModular<int32_t> PIR;
	typedef PIR_ntl_ZZ_p PIR;

	if (argc < 3 or argc > 4) {

		cout << "\nUsage: " << argv[0] << " alg file [m]\n"  << endl;

		cout << " alg = `adaptive', `ilio', `local', or `2local'," << endl
			 << " Modulus m is needed for `local' and `ilio'" << endl 
			 << " m must be a prime power for `local', arbitrary composite for `ilio'." << endl
			 << " Integer smith form is obtained by `ilio' if m is a multiple of the " 
			 << " largest invariant, eg. det." << endl
             << " The matrix is read from file (from cin if file is `-').\n" << endl;
		cout << " Regardless of file format, internal matrix rep is dense." << endl 
			<< " For algoritms using sparse matrix rep, see the examples"
			<< " smithvalence.C, power_rank.C, and poweroftwo_ranks.C." << endl;

		return 0;
	}

	string algo = argv[1];

	string src = argv[2];

	uint64_t m = 1; if (argc == 4) m = atoi(argv[3]);

	UserTimer T;

	if (algo == "adaptive")
	{
		typedef Givaro::ZRing<Integer> Ints;
		Ints Z;
		DenseMatrix<Ints> M(Z);

		Mat(M, src);

		DenseVector<Givaro::ZRing<Integer> > v(Z,M.coldim());
		T.start();
		SmithFormAdaptive::smithForm(v, M);
		T.stop();
		list<pair<integer, size_t> > p;

		distinct(v.begin(), v.end(), p);

		//cout << "#";

		cout << "Integer Smith Form using adaptive alg :\n";
		display(p.begin(), p.end());

		//cout << "# adaptive, Ints, n = " << M.coldim() << endl;

		cout << "T" << M.coldim() << "adaptive(Ints)" << m << " := ";

	}
	else if (algo == "ilio") {

        PIR R( (int32_t)m);

		DenseMatrix<PIR> M(R);

		Mat(M, src);

		T.start();

		SmithFormIliopoulos::smithFormIn (M);
		//PIR::Element d;
		//IliopoulosDomain<PIR> ID(R);
		//ID.smithFormIn (M,d);

		T.stop();

		typedef list< PIR::Element > List;

		List L;

		for (size_t i = 0; i < M.rowdim(); ++i)
			L.push_back(M[(size_t)i][(size_t)i]);

		list<pair<PIR::Element, size_t> > p;

		distinct(L.begin(), L.end(), p);

		//cout << "#";

		cout << "Modular Smith Form using ilio alg :\n";
		display(p.begin(), p.end());

		//cout << "# ilio, PIR-Modular-int32_t(" << m << "), n = " << M.coldim() << endl;

		cout << "T" << M.coldim() << "ilio(PIR-Modular-int32_t)" << m << " := ";
	}

	else if (algo == "local") { // m must be a prime power

		PIR R( (int32_t)m);

		DenseMatrix<PIR> M(R);

		Mat(M, src);

		typedef list< PIR::Element > List;

		List L;

		SmithFormLocal<PIR> SmithForm;

		T.start();

		SmithForm( L, M, R );

		T.stop();

		list<pair<PIR::Element, size_t> > p;

		distinct(L.begin(), L.end(), p);

		//cout << "#";

		PIR::Element x = p.back().first, y;
		R.neg(y,x);
		R.gcdin(x,y);
		if (not R.areEqual(p.back().first, x)) 
			R.write(R.write (cerr << "x ", x) << ", back ", p.back().first) << endl;;
		p.back().first = x;

		cout << "Local Smith Form :\n";
		display(p.begin(), p.end());

		//cout << "# local, PIR-Modular-int32_t(" << m << "), n = " << M.coldim() << endl;

		cout << "T" << M.coldim() << "local(PIR-Modular-int32_t)" << m << " := ";
	}

	else if (algo == "2local") {

		Local2_32 R;

		DenseMatrix<Local2_32> M(R);
		Mat(M, src);

		typedef list< Local2_32::Element > List;

		List L;

		SmithFormLocal<Local2_32> SmithForm;

		T.start();

		SmithForm( L, M, R );

		T.stop();

		list<pair<Local2_32::Element, size_t> > p;

		distinct(L.begin(), L.end(), p);

		//cout << "#";

		cout << "2-Local Smith Form :\n";
		display(p.begin(), p.end());

		//cout << "# 2local, Local2_32, n = " << M.coldim() << endl;

		cout << "T" << M.coldim() << "local2_32 := ";
	}

	else 

		printf ("Unknown algorithm ");

	T.print(cout); cout << /*";" << */ endl;

	return 0 ;
}
template <class PIR> void Mat(DenseMatrix<PIR>& M, string src) {
	if (src[0]=='-') M.read(cin);
	else {
		ifstream in(src);
		M.read(in);
	}
}

//! @bug this already exists elsewhere
template<class I1, class Lp>
void distinct (I1 a, I1 b, Lp& c)
{
	typename iterator_traits<I1>::value_type e;
	size_t count = 0;
	if (a != b) {e = *a; ++a; count = 1;}
	else return;
	while (a != b)
	{  if (*a == e) ++count;
    else
    { c.push_back(typename Lp::value_type(e, count));
    e = *a; count = 1;
    }
    ++a;
	}
	c.push_back(typename Lp::value_type(e, count));
	return;
}

template <class I>
void display(I b, I e)
{ cout << "(";
 for (I p = b; p != e; ++p) cout << "[" << p->first << "," << p->second << "] ";
 cout << ")" << endl;
}
//@}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
