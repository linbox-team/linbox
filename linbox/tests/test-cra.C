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
 * @bug this test fails  for LinBox::FullMultipFixedCRA and LinBox::FullMultipBlasMatCRA
 * @test cra algorithms
 */

#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-early-single.h"
#include "linbox/algorithms/cra-early-multip.h"

#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/cra-full-multip-fixed.h"

#define _LB_ITERS 5

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

// testing EarlySingleCRA
template< class T >
int test_early_single(size_t PrimeSize, size_t Size)
{
	// std::cout << ">" ;

	typedef typename std::vector<T> Vect ;
	typedef typename Vect::iterator Iterator;
	Vect primes(Size) ;
	RandomPrimeIterator RP(PrimeSize);
	/*  primes, probably not all coprime... */
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}

	/*  residues */
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

	cra.result(res);

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

// testing EarlyMultipCRA
template< class T >
int test_early_multip(size_t PrimeSize, size_t Taille, size_t Size)
{
	// std::cout << ">" ;

	typedef typename std::vector<T>                     Vect ;
	typedef typename Vect::iterator                  Iterator;
	typedef typename std::vector<Vect>::iterator VectIterator;
	typedef LinBox::Modular<double>             ModularField ;
	typedef ModularField::Element                     Element;
	typedef std::vector<Integer>                      IntVect;
	typedef std::vector<Element>                        pVect;

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

	Iterator   genprime = primes.begin()    ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	EarlyMultipCRA<ModularField> cra( 4UL ) ;
	IntVect result (Taille); // the result
	pVect residue(Taille) ; // temporary
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
				std::cout << "computed res : "   << result[j]         << std::endl;
				std::cout << "but " << tmp1 << " not equal to " << residues[i][j] << " mod "    << primes[i] << std::endl;
#endif
				return EXIT_FAILURE ;
			}
		}
	}

	// std::cout << "early_cra_multip passed" << std::endl;
	return EXIT_SUCCESS ;
}

// testing FullMultipBlasMatCRA
#if 0
template< class T>
int test_full_multip_matrix(size_t PrimeSize, size_t Size, std::pair<size_t, size_t> dims)
{
	// std::cout << ">" ;

	typedef typename std::vector<T>                    Vect ;
	typedef typename LinBox::BlasMatrix<T>           Matrix ;
	typedef typename std::vector<Matrix>            MatVect ;
	typedef typename Vect::iterator                 Iterator;
	typedef typename MatVect::iterator           MatIterator;
	typedef typename LinBox::BlasMatrix<Integer>  IntMatrix ;

	typedef LinBox::Modular<double>            ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename LinBox::BlasMatrix<Element>    pMatrix ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}

	/*  residues */
	const Matrix Zero (dims.first,dims.second);
	MatVect residues(Size) ;
	for (size_t k = 0 ; k < Size ; ++k) {
		residues[k] = Zero ;
		for (size_t i = 0 ; i < Zero.rowdim() ; ++i)
			for (size_t j = 0 ; j < Zero.coldim() ; ++j)
				residues[k].setEntry( i,j,Integer::random(PrimeSize-1) );
	}

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif

	Iterator  genprime =   primes.begin()  ; // prime iterator
	MatIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = std::log(PrimeSize)*std::log(2)+std::log(Size) ;

	std::pair<size_t,double> my_pair(Size,LogIntSize)  ;

	FullMultipBlasMatCRA<ModularField> cra( my_pair ) ;
	IntMatrix result(dims.first,dims.second); // the result
	pMatrix residue(dims.first,dims.second) ; // temporary
	{ /* init */
		ModularField F(*genprime);
		for (size_t i = 0 ; i < residue.rowdim(); ++i)
			for (size_t j = 0 ; j < residue.coldim(); ++j)
				F.init(residue.refEntry(i,j),(*residu).getEntry(i,j));
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() /*  && !cra.terminated() */ )
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
		for (size_t i = 0 ; i < residue.rowdim(); ++i)
			for (size_t j = 0 ; j < residue.coldim(); ++j)
				F.init(residue.refEntry(i,j),(*residu).getEntry(i,j));

		cra.progress(F,residue);
		++genprime;
		++residu ;
	}

	cra.result(result);

	for (size_t i = 0 ; i < Size ; ++i){
		ModularField F(primes[i]);
		for (size_t l = 0 ; l < dims.first ; ++l)
			for (size_t m = 0 ; m < dims.second ; ++m) {
				Element tmp1,tmp2 ;
				F.init(tmp1,result.getEntry(l,m));
				F.init(tmp2,residues[i].getEntry(l,m));
				if(!F.areEqual(tmp1,tmp2)){
					std::cout << "at " << i << ',' << l << ',' << m << ':' << std::endl;
#if 1
					std::cout << "computed res : "   << result.getEntry(l,m)  << " ( " << tmp1 << ')' << std::endl;
					std::cout << "but " << tmp1 << " not equal to " << residues[i].getEntry(l,m) << " ( " << tmp2 << ')'  << " mod "    << primes[i] << std::endl;
#endif
					return EXIT_FAILURE ;
				}
			}
	}

	std::cout << "full_multip_matrix passed" << std::endl;
	return EXIT_SUCCESS ;
}
#endif

// testing FullMultipCRA
template< class T>
int test_full_multip(size_t PrimeSize, size_t Size, size_t Taille)
{
	// std::cout << ">" ;

	typedef typename std::vector<T>                    Vect ;
	typedef typename std::vector<Vect>             VectVect ;
	typedef std::vector<Integer>                    IntVect ;
	typedef typename Vect::iterator                 Iterator;
	typedef typename VectVect::iterator         VectIterator;

	typedef LinBox::Modular<double >           ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename std::vector<Element>             pVect ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}

	/*  residues */
	VectVect residues(Size) ;
	for (size_t k = 0 ; k < Size ; ++k) {
		residues[k].resize(Taille) ;
		for (size_t i = 0 ; i < Taille ; ++i)
			residues[k][i] = Integer::random(PrimeSize-1) ;
	}

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif

	Iterator   genprime =   primes.begin()  ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = std::log(PrimeSize)*std::log(2)+std::log(Size) ;

	FullMultipCRA<ModularField> cra( LogIntSize ) ;
	IntVect result(Taille) ; // the result
	pVect  residue(Taille) ; // temporary
	{ /* init */
		ModularField F(*genprime);
		for (size_t i = 0 ; i < Taille; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() /* && !cra.terminated()*/ )
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
		for (size_t i = 0 ; i < Taille; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.progress(F,residue);
		++genprime;
		++residu ;
	}

	cra.result(result);

	for (size_t i = 0 ; i < Size ; ++i){
		ModularField F(primes[i]);
		for (size_t j = 0 ; j < Taille ; ++j) {
			Element tmp1,tmp2 ;
			F.init(tmp1,result[j]);
			F.init(tmp2,residues[i][j]);
			if(!F.areEqual(tmp1,tmp2)){
#if 0
				std::cout << "at " << i << ',' << j << ':' << std::endl;
				std::cout << "computed res : "   << result[j]  << " ( " << tmp1 << ')' << std::endl;
				std::cout << "but " << tmp1 << " not equal to " << residues[i][j] << " ( " << tmp2 << ')'  << " mod "    << primes[i] << std::endl;
#endif
				return EXIT_FAILURE ;
			}
		}
	}

	// std::cout << "full_multip passed" << std::endl;
	return EXIT_SUCCESS ;
}


// testing FullMultipFixedCRA
#if 0
template< class T>
int test_full_multip_fixed(size_t PrimeSize, size_t Size, size_t Taille)
{
	// std::cout << ">" ;

	typedef typename std::vector<T>                    Vect ;
	typedef typename std::vector<Vect>             VectVect ;
	typedef std::vector<Integer>                    IntVect ;
	typedef IntVect::iterator                IntVectIterator;
	typedef typename Vect::iterator                 Iterator;
	typedef typename VectVect::iterator         VectIterator;

	typedef LinBox::Modular<double >           ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename std::vector<Element>             pVect ;
	typedef typename pVect::iterator          pVectIterator ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	RandomPrimeIterator RP(PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = RP.randomPrime() ;
		++RP ;
	}

	/*  residues */
	VectVect residues(Size) ;
	for (size_t k = 0 ; k < Size ; ++k) {
		residues[k].resize(Taille) ;
		for (size_t i = 0 ; i < Taille ; ++i)
			residues[k][i] = Integer::random(PrimeSize-1) ;
	}

#if 0
	std::cout << "primes : " << std::endl ;
	std::cout << primes << std::endl;
	std::cout << "residues : " << std::endl ;
	std::cout << residues << std::endl;
#endif

	Iterator   genprime =   primes.begin()  ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = std::log(PrimeSize)*std::log(2)+std::log(Size) ;
	std::pair<size_t,double> my_pair(Size,LogIntSize)  ;

	FullMultipFixedCRA<ModularField> cra( my_pair ) ;
	IntVect result(Taille) ; // the result
	pVect  residue(Taille) ; // temporary
	pVectIterator residue_it = residue.begin();
	{ /* init */
		ModularField F(*genprime);
		for (size_t i = 0 ; i < Taille; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.initialize(F,residue_it);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() /* && !cra.terminated()*/ )
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
		pVectIterator residue_it = residue.begin();
		for (size_t i = 0 ; i < Taille; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.progress(F,residue_it);
		++genprime;
		++residu ;
	}

	IntVectIterator result_it = result.begin();

	cra.result(result_it);

	for (size_t i = 0 ; i < Size ; ++i){
		ModularField F(primes[i]);
		for (size_t j = 0 ; j < Taille ; ++j) {
			Element tmp1,tmp2 ;
			F.init(tmp1,result[j]);
			F.init(tmp2,residues[i][j]);
			if(!F.areEqual(tmp1,tmp2)){
#if 1
				std::cout << "at " << i << ',' << j << ':' << std::endl;
				std::cout << "computed res : "   << result[j]  << " ( " << tmp1 << ')' << std::endl;
				std::cout << "but " << tmp1 << " not equal to " << residues[i][j] << " ( " << tmp2 << ')'  << " mod "    << primes[i] << std::endl;
#endif
				return EXIT_FAILURE ;
			}
		}
	}

	std::cout << "full_multip passed" << std::endl;
	return EXIT_SUCCESS ;
}
#endif

// launching tests
int main(int ac, char ** av)
{

	typedef std::pair<size_t,size_t> Pair ;
	// srand(time(NULL));
	size_t PrimeSize = 22; // size of the residues/primes
	size_t Size      = 50 ; // nb of residues/primes
	size_t Taille    = 2*Size ; // nb of vectors of residues
	if (ac > 1) Size = atoi(av[1]);
	//if (ac > 2) PrimeSize = atoi(av[2]);
	// std::cout << Size << ',' << PrimeSize << std::endl;
	/* EARLY SINGLE */
	_LB_REPEAT( if (test_early_single<double>(22,Size))                   return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_early_single<integer>(PrimeSize,Size))           return EXIT_FAILURE ;  ) ;

	/* EARLY MULTIPLE */
	_LB_REPEAT( if (test_early_multip<double>(22,Taille*2,Size))          return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_early_multip<integer>(PrimeSize,Taille*2,Size))  return EXIT_FAILURE ;  ) ;

	_LB_REPEAT( if (test_early_multip<double>(22,Taille/4,Size))          return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_early_multip<integer>(PrimeSize,Taille/4,Size))  return EXIT_FAILURE ;  ) ;

	/* FULL MULTIPLE */
	_LB_REPEAT( if (test_full_multip<double>(22,Size,Taille))             return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip<integer>(PrimeSize,Size,Taille))     return EXIT_FAILURE ;  ) ;

	_LB_REPEAT( if (test_full_multip<double>(22,Size,Taille/4))           return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip<integer>(PrimeSize,Size,Taille/4))   return EXIT_FAILURE ;  ) ;

#if 0 /* FULL MULTIPLE FIXED */
	_LB_REPEAT( if (test_full_multip_fixed<double>(22,Size,Taille))           return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip_fixed<integer>(PrimeSize,Size,Taille))   return EXIT_FAILURE ;  ) ;

	_LB_REPEAT( if (test_full_multip_fixed<double>(22,Size,Taille/4))         return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip_fixed<integer>(PrimeSize,Size,Taille/4)) return EXIT_FAILURE ;  ) ;
#endif


#if 0
	/* FULL MULTIPLE MATRIX */
	Pair p(Taille,2*Taille);
	_LB_REPEAT( if (test_full_multip_matrix<double>(22,Size,p))           return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip_matrix<integer>(PrimeSize,Size,p))   return EXIT_FAILURE ;  ) ;

	Pair q(Taille*2,Taille);
	_LB_REPEAT( if (test_full_multip_matrix<double>(22,Size,q))           return EXIT_FAILURE ;  ) ;
	_LB_REPEAT( if (test_full_multip_matrix<integer>(PrimeSize,Size,q))   return EXIT_FAILURE ;  ) ;

#endif
	return EXIT_SUCCESS ;
}
