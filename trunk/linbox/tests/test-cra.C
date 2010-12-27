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

#define _LB_REPEAT(command) \
do { for (size_t i = 0 ; pass && i < iters ; ++i) {  command } } while(0)

using namespace LinBox ;

template< class T >
std::ostream & operator<<(std::ostream & o, const std::vector<T> & v)
{
	for (size_t i = 0 ; i < v.size() ; ++i )
		o << v[i] << ' ';
	return o;
}

// testing EarlySingleCRA
template< class T >
int test_early_single(std::ostream & report, size_t PrimeSize, size_t Size)
{

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

	typedef LinBox::Modular<double> ModularField ;

	Iterator genprime = primes.begin()  ; // prime iterator
	Iterator residu = residues.begin()  ; // residu iterator

	report << "EarlySingleCRA (" <<  4UL << ')' << std::endl;
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
		if (cra.noncoprime((integer)*genprime)) {
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "EarlySingleCRA exiting successfully." << std::endl;
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
			report << tmp1 << "!=" << tmp2 << std::endl;
			report << " *** EarlySingleCRA failed. ***" << std::endl;
			return EXIT_FAILURE ;
		}
	}

	report << "EarlySingleCRA exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}

// testing EarlyMultipCRA
template< class T >
int test_early_multip(std::ostream & report, size_t PrimeSize, size_t Taille, size_t Size)
{

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


	Iterator   genprime = primes.begin()    ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	report << "EarlyMultpCRA (" <<  4UL << ')' << std::endl;
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
		if (cra.noncoprime((integer)*genprime)) {
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "EarlyMultipCRA exiting successfully." << std::endl;
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
				report << " *** EarlyMultipCRA failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "EarlyMultipCRA exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}


#if 1 /* testing FullMultipBlasMatCRA */
template< class T>
int test_full_multip_matrix(std::ostream & report, size_t PrimeSize, size_t Size, std::pair<size_t, size_t> dims)
{

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

	Iterator  genprime =   primes.begin()  ; // prime iterator
	MatIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = (PrimeSize+1)*std::log(2)+std::log(Size)+1 ;

	std::pair<size_t,double> my_pair(dims.first*dims.second,LogIntSize)  ;

	report << "FullMultipBlasMatCRA (" <<  my_pair.first << ", " << my_pair.second << ')' << std::endl;
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
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "FullMultipBlasMatCRA exiting successfully." << std::endl;
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
					report << result.getEntry(l,m) << ';' << residues[i].getEntry(l,m) << '@' << primes[i] << std::endl;
					report << i << ':' << l << ',' << m << "> " << tmp1 << "!=" << tmp2 << std::endl;
					report << " *** FullMultipBlasMatCRA failed. ***" << std::endl;
					return EXIT_FAILURE ;
				}
			}
	}

	report << "FullMultipBlasMatCRA exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}
#endif

// testing FullMultipCRA
template< class T>
int test_full_multip(std::ostream & report, size_t PrimeSize, size_t Size, size_t Taille)
{

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


	Iterator   genprime =   primes.begin()  ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = PrimeSize*std::log(2)+std::log(Size)+1 ;

	report << "FullMultipCRA (" <<  LogIntSize << ')' << std::endl;
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
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "FullMultipCRA exiting successfully." << std::endl;
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
				report << " *** FullMultipCRA failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "FullMultipCRA exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}



#if 1 /* testing FullMultipFixedCRA */
template< class T>
int test_full_multip_fixed(std::ostream & report, size_t PrimeSize, size_t Size, size_t Taille)
{


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

	Iterator   genprime =   primes.begin()  ; // prime iterator
	VectIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = PrimeSize*std::log(2)+std::log(Size) ;

	std::pair<size_t,double> my_pair(Taille,LogIntSize)  ;

	report << "FullMultipFixedCRA (" <<  my_pair.first << ", " << my_pair.second << ')' << std::endl;

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
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "FullMultipFixedCRA exiting successfully." << std::endl;
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
				report << " *** FullMultipFixedCRA failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "FullMultipFixedCRA exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}
#endif

bool test_CRA_algos(size_t PrimeSize, size_t Size, size_t Taille, size_t iters)
{
	bool pass = true ;
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT,
							   INTERNAL_DESCRIPTION);



	typedef std::pair<size_t,size_t> Pair ;

	/* EARLY SINGLE */
	_LB_REPEAT( if (test_early_single<double>(report,22,Size))                       pass = false ;  ) ;
	_LB_REPEAT( if (test_early_single<integer>(report,PrimeSize,Size))               pass = false ;  ) ;

	/* EARLY MULTIPLE */
	_LB_REPEAT( if (test_early_multip<double>(report,22,Taille*2,Size))              pass = false ;  ) ;
	_LB_REPEAT( if (test_early_multip<integer>(report,PrimeSize,Taille*2,Size))      pass = false ;  ) ;

	_LB_REPEAT( if (test_early_multip<double>(report,22,Taille/4,Size))              pass = false ;  ) ;
	_LB_REPEAT( if (test_early_multip<integer>(report,PrimeSize,Taille/4,Size))      pass = false ;  ) ;

	/* FULL MULTIPLE */
	_LB_REPEAT( if (test_full_multip<double>(report,22,Size,Taille))                 pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip<integer>(report,PrimeSize,Size,Taille))         pass = false ;  ) ;

	_LB_REPEAT( if (test_full_multip<double>(report,22,Size,Taille/4))               pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip<integer>(report,PrimeSize,Size,Taille/4))       pass = false ;  ) ;

#if 1 /* FULL MULTIPLE FIXED */
	_LB_REPEAT( if (test_full_multip_fixed<double>(report,22,Size,Taille))           pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip_fixed<integer>(report,PrimeSize,Size,Taille))   pass = false ;  ) ;

	_LB_REPEAT( if (test_full_multip_fixed<double>(report,22,Size,Taille/4))         pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip_fixed<integer>(report,PrimeSize,Size,Taille/4)) pass = false ;  ) ;
#endif


#if 1 /* FULL MULTIPLE MATRIX */
	Taille = 15 ;
	Pair q(Taille,2*Taille);
	_LB_REPEAT( if (test_full_multip_matrix<double>(report,22,Size,q))               pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip_matrix<integer>(report,PrimeSize,Size,q))       pass = false ;  ) ;

	Pair s(Taille*2,Taille);
	_LB_REPEAT( if (test_full_multip_matrix<double>(report,22,Size,s))               pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip_matrix<integer>(report,PrimeSize,Size,s))       pass = false ;  ) ;

#endif

	return pass ;

}

#include "test-common.h"
#include "linbox/util/timer.h"

// launching tests
int main(int ac, char ** av)
{

	/*  Argument parsing/setting */

	static size_t       n = 50;    /*  Taille */
	static size_t       p = 22;    /*  PrimeSize */
	// static size_t    seed =  0;    /*  ! unused */
	static size_t   iters = 20;    /* _LB_REPEAT */

        static Argument as[] = {
                { 'n', "-n N", "Set number of primes.", TYPE_INT , &n },
                { 'p', "-p P", "Set size of test primes.", TYPE_INT , &p },
                { 'i', "-i I", "Perform each test for I iterations.",     TYPE_INT, &iters },
		END_OF_ARGUMENTS
        };

	parseArguments (ac, av, as);

	bool pass = true ;

	srand(time(NULL));             // seeding
	size_t PrimeSize   =  p;       // size of the residues/primes
	size_t Size        =  n ;      // nb of residues/primes
	size_t Taille      =  2*Size ; // nb of vectors of residues


	commentator.start("CRA-Algos test suite", "CRA-Algos");

	pass = test_CRA_algos(PrimeSize,Size,Taille,iters) ;

	commentator.stop(MSG_STATUS (pass), (const char *) 0,"CRA-Algos test suite");
	return !pass ;
}
