/*
 * examples/smithsparse.C
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

/** \file examples/smithsparse.C
 * @example  examples/smithsparse.C
 \brief Smith form by sparse elmination, Integer Smith by valence method, or local at a prime power. 
 \ingroup examples

 See smithvalence for valence method with more options and information.

 For local smith, moduli greater than 2^32 are not supported here (easily changed).
 matrix is read from file or generated from a small selection of example families.
 Run the program with no arguments for a synopsis of the
 command line parameters.

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

#include <linbox/ring/pir-modular-int32.h>
#define NOT_USING_OMP
#include "smithvalence.h"
#undef NOT_USING_OMP

using namespace LinBox;

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

int main(int argc, char* argv[])
{
	Givaro::Timer chrono; 

	if (argc < 2 or 3 < argc) {
		cout << "Usage: " << argv[0] << " file [m]"  << endl;
		cout << "where file contains the input matrix and m is the modulus, a prime power." << endl; 
		cout << "If file is "-", cin is used" << endl;
		cout << "If m is given, a sparse elimination over Z_m is performed." << endl;
		cout << "if m is absent, the valence method is used.  See smithvalence for more options and info." << endl;
		cout << "See mats.C for some examples that have been used in smith form algorithm testing" << endl;
		return 0;
	}

	chrono.start();

	if (argc > 2) { // so over Z_m
		unsigned long m = atoi(argv[2]);

		typedef Givaro::Modular<int32_t> SPIR;
		SPIR R(m);

		SparseMatrix<SPIR, SparseMatrixFormat::SparseSeq > B (R);
		if (argv[1] == "-") 
			B.read(std::cin);
		else {
			std::ifstream input (argv[1]);
			if (!input) { std::cerr << "Error opening matrix file: " << argv[2] << std::endl; return -1; }
			B.read(input);
		}
		std::cout << "matrix is " << B.rowdim() << " by " << B.coldim() << std::endl;
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::cout) << std::endl;

		Integer p(m), im(m);
                // Should better ask user to give the prime !!!
    	Givaro::IntPrimeDom IPD;
		for(unsigned int k = 2; ( ( ! IPD.isprime(p) ) && (p > 1) ); ++k)
       		Givaro::root( p, im, k );

    	// using Sparse Elimination
		LinBox::PowerGaussDomain< SPIR > PGD( R );
		std::vector<std::pair<size_t,SPIR::Element> > vec;
    	LinBox::Permutation<SPIR> Q(R,B.coldim());

		PGD(vec, B, Q, (int32_t)m, (int32_t)p);

		typedef list< SPIR::Element > List;
		List L;
		for ( auto p_it = vec.begin(); p_it != vec.end(); ++p_it) {
			for(size_t i = 0; i < (size_t) p_it->first; ++i)
				L.push_back((SPIR::Element)p_it->second);
		}
		size_t M = (B.rowdim() > B.coldim() ? B.coldim() : B.rowdim());
		size_t Min = (B.rowdim() < B.coldim() ? B.coldim() : B.rowdim());
		for (size_t i = L.size(); i < M; ++i)
			L.push_back(0);

		list<pair<SPIR::Element, size_t> > pl;

		distinct(L.begin(), L.end(), pl);

		//std::cout << "#";

         //display(local.begin(), local.end());
		display(pl.begin(), pl.end());
		//cout << "# local, PowerGaussDomain<int32_t>(" << M << "), n = " << Min << endl;


		chrono.stop();
		cout << "T" << M << "local(PowerGaussDomain<int32_t>)" << m << " := "
			<< chrono << std::endl;

	} else {// argc is 2, so use valence method over ZZ

		Givaro::ZRing<Integer> ZZ;
		typedef SparseMatrix<Givaro::ZRing<Integer> >  Blackbox;
		Blackbox A (ZZ);
	
		std::ifstream input (argv[1]);
		if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }
		A.read(input);
		input.close();
	
		std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
	
		Givaro::ZRing<Integer>::Element val_A;
	
		chrono.start();
		Transpose<Blackbox> T(&A);
		if (A.rowdim() > A.coldim()) {//ata
			Compose< Transpose<Blackbox>, Blackbox > C (&T, &A);
		//	std::cout << "A^T A is " << C.rowdim() << " by " << C.coldim() << std::endl;
			valence(val_A, C);
		}
		else if (A.rowdim() < A.coldim()) {//aat
			Compose< Blackbox, Transpose<Blackbox> > C (&A, &T);
		//	std::cout << "A A^T is " << C.rowdim() << " by " << C.coldim() << std::endl;
			valence(val_A, C);
		}
		else { // square, just a
			valence(val_A, A);
		}
	
		//std::cout << "Valence is " << val_A << std::endl;
	
		std::vector<integer> Moduli;
		std::vector<size_t> exponents;
	 	Givaro::IntFactorDom<> FTD;
	
		typedef std::pair<integer,unsigned long> PairIntRk;
		std::vector< PairIntRk > smith;
	
	
		integer coprimeV=2;
		while ( gcd(val_A,coprimeV) > 1 ) {
			FTD.nextprimein(coprimeV);
		}
	
		//std::cout << "integer rank: " << std::endl;
	
		unsigned long coprimeR; LRank(coprimeR, argv[1], coprimeV);
		smith.push_back(PairIntRk(coprimeV, coprimeR));
		//         std::cerr << "Rank mod " << coprimeV << " is " << coprimeR << std::endl;
	
		//std::cout << "Some factors (5000 factoring loop bound): ";
		FTD.set(Moduli, exponents, val_A, 5000);
		std::vector<size_t>::const_iterator eit=exponents.begin();
		//for(std::vector<integer>::const_iterator mit=Moduli.begin();
		//    mit != Moduli.end(); ++mit,++eit)
		//	std::cout << *mit << '^' << *eit << ' ';
		//std::cout << std::endl;
	
		std::vector<integer> SmithDiagonal(coprimeR,integer(1));
	
		for(std::vector<integer>::const_iterator mit=Moduli.begin();
		    mit != Moduli.end(); ++mit) {
			unsigned long r; LRank(r, argv[1], *mit);
			//             std::cerr << "Rank mod " << *mit << " is " << r << std::endl;
			smith.push_back(PairIntRk(*mit, r));
			for(size_t i=r; i < coprimeR; ++i)
				SmithDiagonal[i] *= *mit;
		}
	
	
		eit=exponents.begin();
		std::vector<PairIntRk>::const_iterator sit=smith.begin();
		for( ++sit; sit != smith.end(); ++sit, ++eit) {
			if (sit->second != coprimeR) {
				std::vector<size_t> ranks;
				ranks.push_back(sit->second);
	            size_t effexp;
				if (*eit > 1) {
					PRank(ranks, effexp, argv[1], sit->first, *eit, coprimeR);
				}
				else {
					PRank(ranks, effexp, argv[1], sit->first, 2, coprimeR);
				}
				if (ranks.size() == 1) ranks.push_back(coprimeR);
	
	            if (effexp < *eit) {
	                for(size_t expo = effexp<<1; ranks.back() < coprimeR; expo<<=1) {
	                    PRankInteger(ranks, argv[1], sit->first, expo, coprimeR);
	                }
	            } else {
	
	                for(size_t expo = (*eit)<<1; ranks.back() < coprimeR; expo<<=1) {
	                    PRank(ranks, effexp, argv[1], sit->first, expo, coprimeR);
	                    if (ranks.size() < expo) {
	     //                   std::cerr << "It seems we need a larger prime power, it will take longer ..." << std::endl;
	                            // break;
	                        PRankInteger(ranks, argv[1], sit->first, expo, coprimeR);
	                    }
	                }
	            }
	
				std::vector<size_t>::const_iterator rit=ranks.begin();
				// unsigned long modrank = *rit;
				for(++rit; rit!= ranks.end(); ++rit) {
					if ((*rit)>= coprimeR) break;
					for(size_t i=(*rit); i < coprimeR; ++i)
						SmithDiagonal[i] *= sit->first;
					// modrank = *rit;
				}
			}
		}
	
		integer si=1;
		size_t num=0;
		std::cout << '(';
		for( std::vector<integer>::const_iterator dit=SmithDiagonal.begin();
		     dit != SmithDiagonal.end(); ++dit) {
			if (*dit == si) ++num;
			else {
				if (num > 0)
					std::cout << si << ' ' << num << ", ";
				num=1;
				si = *dit;
			}
		}
		std::cout << si << ' ' << num << ") " << std::endl;
		chrono.stop();
		cout << "T" << A.rowdim() << "smithvalence(ZRing<Integer>):= "
			<< chrono << std::endl;
	
	}
	return 0;
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


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
