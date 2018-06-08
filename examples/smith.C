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
 The matrix example generation code that was here is now in mats.C.
*/

#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;


#include <linbox/ring/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>

#include <linbox/util/timer.h>

#include <linbox/ring/local2_32.h>
// place A: Edit here and at place B for ring change
//#include <linbox/ring/pir-modular-int32.h>
#include <linbox/ring/pir-ntl-zz_p.h>
#include <linbox/algorithms/smith-form-local.h>
#include <linbox/algorithms/smith-form-iliopoulos.h>
#include <linbox/algorithms/smith-form-adaptive.h>

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

		cout << "usage: " << argv[0] << " alg file [m]"  << endl;

		cout << " alg = `adaptive', `ilio', `local', or `2local'," << endl
			 << " Modulus m is needed for `local' and `ilio'" << endl 
			 << " m must be a prime power for `local', arbitrary composite for `ilio'." << endl
			 << " Integer smith form is obtained by `ilio' if m is a multiple of largest invariant, eg. det." << endl
             << " Matrix is read from file, from cin if file is `-'." << endl;

		return 0;
	}

	string algo = argv[1];

	string src = argv[2];

	unsigned long m = 1; if (argc == 4) m = atoi(argv[3]);

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

		display(p.begin(), p.end());

		//cout << "# ilio, PIR-Modular-int32_t(" << m << "), n = " << M.coldim() << endl;

		cout << "T" << M.coldim() << "ilio(PIR-Modular-int32_t)" << m << " := ";
	}

	else if (algo == "local") { // m must be a prime power

#if 0
		if (format == "sparse" ) {
			typedef Givaro::Modular<int32_t> Field;
			Field F(m);
			std::ifstream input (argv[4]);
			if (!input) { std::cerr << "Error opening matrix file: " << argv[1] << std::endl; return -1; }

			MatrixStream<Field> ms( F, input );
			SparseMatrix<Field, SparseMatrixFormat::SparseSeq > B (ms);
			std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
			if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::cout) << std::endl;



			Integer p(m), im(m);
                // Should better ask user to give the prime !!!
            Givaro::IntPrimeDom IPD;
			for(unsigned int k = 2; ( ( ! IPD.isprime(p) ) && (p > 1) ); ++k)
                Givaro::root( p, im, k );

                // using Sparse Elimination
			LinBox::PowerGaussDomain< Field > PGD( F );
			std::vector<std::pair<size_t,Field::Element> > local;
            LinBox::Permutation<Field> Q(F,B.coldim());

			PGD(local, B, Q, (int32_t)m, (int32_t)p);

			typedef list< Field::Element > List;
			List L;
			for ( auto p_it = local.begin(); p_it != local.end(); ++p_it) {
				for(size_t i = 0; i < (size_t) p_it->first; ++i)
					L.push_back((Field::Element)p_it->second);
			}
			size_t M = (B.rowdim() > B.coldim() ? B.coldim() : B.rowdim());
			for (size_t i = L.size(); i < M; ++i)
				L.push_back(0);

			list<pair<Field::Element, size_t> > pl;

			distinct(L.begin(), L.end(), pl);

			std::cout << "#";

                //display(local.begin(), local.end());
			display(pl.begin(), pl.end());
			cout << "# local, PowerGaussDomain<int32_t>(" << M << "), n = " << n << endl;

		}
		else {
#endif

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
 for (I p = b; p != e; ++p) cout << p->first << " " << p->second << ", ";
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
