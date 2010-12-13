/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox
 *
 *  Written by Brice Boyer <brice.boyer@imag.fr>
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*! @file  tests/test-cra.C
 * @ingroup tests
 * @ingroup CRA
 * @brief We test the various CRA algorithms here.
 */

#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-early-single.h"

using namespace LinBox ;

struct DummyIteration {
	double x_ ;
	DummyIteration(double x) :
		x_(x)
	{}
	~DummyIteration() {}
} ;

std::ostream & operator<<(std::ostream & o, const std::vector<double> & v)
{
	for (size_t i = 0 ; i < v.size() ; ++i )
		o << v[i] << ' ';
	return o;
}

int test_early_single(int ac, char ** av)
{

	size_t PrimeSize = 20 ; //! @todo why ???

	size_t n = 10 ;
	typedef std::vector<double> Vect ;
	Vect primes(n) ;
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < n ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}


	/*  very probably all coprime... */

	Vect residues(n) ;
	for (size_t i = 0 ; i < n ; ++i)
		residues[i] = Integer::random(PrimeSize-1);

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif
	typedef LinBox::Modular<double> ModularField ;

	Vect::iterator genprime = primes.begin()  ; // prime iterator
	Vect::iterator residu = residues.begin()  ; // prime iterator

	EarlySingleCRA<ModularField> cra( 4UL ) ; //!\todo what is \p EARLY ?
	Integer res = 0;
	typedef ModularField::Element Element;
	Element residue ;
	{ // init
		ModularField F(*genprime);
		F.init(residue,*residu);
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() && !cra.terminated() ) { // progress
		if (cra.noncoprime((integer)*genprime))
		{
#if 0
			std::cout << "bad luck, you picked twice the same prime..." <<std::endl;
			std::cout << *genprime << std::endl;
#endif
			return EXIT_SUCCESS ; // pas la faute à cra...
		}
		ModularField F(*genprime);
		F.init(residue,*residu);
		cra.progress(F,residue);
		++genprime;
		++residu ;
	}

	cra.getResidue(res);

	for (size_t i = 0 ; i < n ; ++i){
		ModularField F(primes[i]);
		Element tmp1,tmp2 ;
		F.init(tmp1,res);
		F.init(tmp2,residues[i]);
		if(!F.areEqual(tmp1,tmp2)){
#if 0
			std::cout << "computed res : "   << res         << std::endl;
			std::cout << "but " << tmp1 << " not equal to " << residues[i] << " mod "    << primes[i] << std::endl;
#endif
			return EXIT_FAILURE ;
		}
	}

	std::cout << "early_cra_single passed" << std::endl;
	return EXIT_SUCCESS ;
}

int main(int ac, char ** av)
{
	if (test_early_single(ac,av)) return EXIT_FAILURE ;
	return EXIT_SUCCESS ;
}
