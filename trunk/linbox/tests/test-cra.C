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
 * @bug this test fails sometimes
 */

#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-early-single.h"
#include "linbox/algorithms/cra-early-multip.h"

#define _LB_ITERS 50

#define _LB_REPEAT(command) \
do { for (size_t i = 0 ; i < _LB_ITERS ; ++i) {  command } } while(0)

using namespace LinBox ;

struct DummyIteration {
	double x_ ;
	DummyIteration(double x) :
		x_(x)
	{}
	~DummyIteration() {}
} ;

template< class T >
std::ostream & operator<<(std::ostream & o, const std::vector<T> & v)
{
	for (size_t i = 0 ; i < v.size() ; ++i )
		o << v[i] << ' ';
	return o;
}

template< class T >
int test_early_single(size_t PrimeSize, size_t Size)
{
	// std::cout << ">" ;

	typedef typename std::vector<T> Vect ;
	typedef typename Vect::iterator Iterator;
	Vect primes(Size) ;
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}


	/*  probably not all coprime... */

	Vect residues(Size) ;
	for (size_t i = 0 ; i < Size ; ++i)
		residues[i] = Integer::random(PrimeSize-1);

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif
	typedef LinBox::Modular<double> ModularField ;

	Iterator genprime = primes.begin()  ; // prime iterator
	Iterator residu = residues.begin()  ; // residu iterator

	EarlySingleCRA<ModularField> cra( 4UL ) ;
	Integer res = 0; // the result
	typedef ModularField::Element Element;
	Element residue ; // temporary
	{ /* init */
		ModularField F(*genprime);
		F.init(residue,*residu);
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() && !cra.terminated() )
	{ /* progress */
		if (cra.noncoprime((integer)*genprime))
		{
#if 0
			std::cout << "bad luck, you picked twice the same prime..." <<std::endl;
			std::cout << *genprime << std::endl;
			for (size_t i = 0 ; i < Size ; ++i) {
				if (*genprime == primes[i])
					std::cout << "was n°" << i << std::endl;
			}
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

	for (size_t i = 0 ; i < Size ; ++i){
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

	// std::cout << "early_cra_single passed" << std::endl;
	return EXIT_SUCCESS ;
}

template< class T >
int test_early_multip(size_t PrimeSize, size_t Taille, size_t Size)
{
	// std::cout << ">" ;

	typedef typename std::vector<T> Vect ;
	typedef typename Vect::iterator Iterator;
	typedef typename std::vector<Vect>::iterator VectIterator;

	/*  primes */
	Vect primes(Size) ;
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}

	/*  residues */
	std::vector<Vect> residues(Size) ;
	for (size_t i = 0 ; i < Size ; ++i){
		residues[i].resize(Taille);
		for (size_t j = 0 ; j < Taille ; ++j)
			residues[i][j] = Integer::random(PrimeSize-1);
	}

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif
	typedef LinBox::Modular<double> ModularField ;

	Iterator   genprime = primes.begin()    ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	EarlyMultipCRA<ModularField> cra( 4UL ) ;
	std::vector<Integer> result (Taille); // the result
	typedef ModularField::Element Element;
	std::vector<Element> residue(Taille) ; // temporary
	{ /* init */
		ModularField F(*genprime);
		for (size_t i = 0 ; i < Taille ; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() && !cra.terminated() )
	{ /* progress */
		if (cra.noncoprime((integer)*genprime))
		{
#if 0
			std::cout << "bad luck, you picked twice the same prime..." <<std::endl;
			// std::cout << *genprime << std::endl;
			// for (size_t i = 0 ; i < Size ; ++i) {
				// if (*genprime == primes[i])
					// std::cout << "was n°" << i << std::endl;
			// }
#endif
			return EXIT_SUCCESS ; // pas la faute à cra...
		}
		ModularField F(*genprime);
		for (size_t i = 0 ; i < Taille ; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.progress(F,residue);
		++genprime;
		++residu ;
	}

	cra.result(result);

	for (size_t i = 0 ; i < Size ; ++i){
			ModularField F(primes[i]);
		for (size_t j = 0 ; j < Taille ; ++j){
			Element tmp1,tmp2 ;
			F.init(tmp1,result[j]);
			F.init(tmp2,residues[i][j]);
			if(!F.areEqual(tmp1,tmp2)){
#if 0
				std::cout << "computed res : "   << result[i]         << std::endl;
				std::cout << "but " << tmp1 << " not equal to " << residues[i][j] << " mod "    << primes[i] << std::endl;
#endif
				return EXIT_FAILURE ;
			}
		}
		}

	// std::cout << "early_cra_multip passed" << std::endl;
	return EXIT_SUCCESS ;
}


int main(int ac, char ** av)
{
	// srand(time(NULL));
	size_t PrimeSize = 22; // size of the residues/primes
	size_t Size      = 50 ; // nb of residues/primes
	size_t Taille    = 2*Size ; // nb of vectors of residues
	if (ac > 1) Size = atoi(av[1]);
	//if (ac > 2) PrimeSize = atoi(av[2]);
	// std::cout << Size << ',' << PrimeSize << std::endl;
	/* EARLY SINGLE */
	_LB_REPEAT( if (test_early_single<double>(22,Size))  return EXIT_FAILURE ;   );
	_LB_REPEAT( if (test_early_single<integer>(PrimeSize,Size))  return EXIT_FAILURE ;  ) ;

	/* EARLY MULTIPLE */
	_LB_REPEAT( if (test_early_multip<double>(22,Taille,Size))  return EXIT_FAILURE ;   );
	_LB_REPEAT( if (test_early_multip<integer>(PrimeSize,Taille,Size))  return EXIT_FAILURE ;  ) ;

	_LB_REPEAT( if (test_early_multip<double>(22,Taille/4,Size))  return EXIT_FAILURE ;   );
	_LB_REPEAT( if (test_early_multip<integer>(PrimeSize,Taille/4,Size))  return EXIT_FAILURE ;  ) ;


	return EXIT_SUCCESS ;
}
