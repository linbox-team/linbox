/* Copyright (C) LinBox
 *
 *  Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-cra.C
 * @ingroup tests
 * @ingroup CRA
 * @brief We test the various CRA algorithms here.
 * @test cra algorithms
 */

#include "linbox/linbox-config.h"
#include <givaro/zring.h>
#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-builder-single.h"
#include "linbox/algorithms/cra-builder-early-multip.h"
#include "linbox/algorithms/rational-cra-builder-full-multip.h"

#include "linbox/matrix/dense-matrix.h"
#include "linbox/algorithms/cra-builder-full-multip.h"
#include "linbox/algorithms/cra-builder-full-multip-fixed.h"


#define _LB_REPEAT(command) \
do { for (size_t i = 0 ; pass && i < iters ; ++i) {  command } } while(0)

using namespace LinBox ;


// Need these call_* functions so it uses Integers with primes
// larger than 23 bits.
template <typename CRAType>
void call_initialize(CRAType& cra, const double p, const double r) {
	using ModularField = typename CRAType::Domain;
	using Element = typename ModularField::Element;
	ModularField F(p);
	Element residue;
	F.init(residue, r);
	cra.initialize(F, residue);
}

template <typename CRAType>
void call_initialize(CRAType& cra, const Integer& p, const Integer& r) {
	cra.initialize(p, r);
}

template <typename CRAType>
void call_progress(CRAType& cra, const double p, const double r) {
	using ModularField = typename CRAType::Domain;
	using Element = typename ModularField::Element;
	ModularField F(p);
	Element residue;
	F.init(residue, r);
	cra.progress(F, residue);
}

template <typename CRAType>
void call_progress(CRAType& cra, const Integer& p, const Integer& r) {
	cra.progress(p, r);
}

// testing CRABuilderEarlySingle
template< class T >
int test_early_single(std::ostream & report, size_t PrimeSize, size_t Size)
{

	typedef typename std::vector<T> Vect ;
	typedef typename Vect::iterator Iterator;
	Vect primes(Size) ;
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	/*  primes, probably not all coprime... */
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
		++RP ;
	}

	/*  residues */
	Vect residues(Size) ;
	for (size_t i = 0 ; i < Size ; ++i)
		residues[i] = Integer::random(PrimeSize-1);

	typedef Givaro::Modular<double> ModularField ;

	Iterator genprime = primes.begin()  ; // prime iterator
	Iterator residu = residues.begin()  ; // residu iterator

	report << "CRABuilderEarlySingle (" <<  LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD << ')' << std::endl;
	CRABuilderEarlySingle<ModularField> cra( LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD ) ;
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
			report << "CRABuilderEarlySingle exiting successfully." << std::endl;
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
			report << " *** CRABuilderEarlySingle failed. ***" << std::endl;
			return EXIT_FAILURE ;
		}
	}

	report << "CRABuilderEarlySingle exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}

// testing CRABuilderProbSingle
template< class T >
int test_prob_single(std::ostream & report, size_t PrimeSize, size_t Size)
{

	typedef typename std::vector<T> Vect ;
	typedef typename Vect::iterator Iterator;

        Integer pprod(1); // product of distinct primes
	Vect primes(Size) ;
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	/*  primes, probably not all coprime... */
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
		++RP ;
                if (pprod % primes[i]) {
			pprod *= primes[i];
                }
	}

	// true result
	size_t resbits = 1 + (random() % (pprod.bitsize() - 1));
	Integer actual = Integer::random(resbits);

	/*  residues */
	Vect residues(Size) ;
	for (size_t i = 0 ; i < Size ; ++i)
		residues[i] = actual % primes[i];

	typedef Givaro::Modular<double> ModularField ;

	Iterator genprime = primes.begin()  ; // prime iterator
	Iterator residu = residues.begin()  ; // residu iterator

	report << "CRABuilderProbSingle (" << pprod.bitsize()-1 << ")";
	report << " actual length " << actual.bitsize() << std::endl;
	CRABuilderProbSingle<ModularField> cra(pprod.bitsize()-1) ;
	Integer res = 0; // the result
	typedef ModularField::Element Element;
	{ /* init */
		call_initialize(cra, *genprime, *residu);
	}
	size_t itercount = 1;
	size_t skips = 0;
	while (genprime < primes.end() && !cra.terminated() )
	{ /* progress */
		if (cra.noncoprime((integer)*genprime)) {
			//report << "bad luck, you picked twice the same prime..." <<std::endl;
			++skips;
		}
		else {
			call_progress(cra, *genprime, *residu);
			++itercount;
		}
		++genprime;
		++residu ;
	}
	report << "  " << itercount << " iterations, " << itercount*(PrimeSize-1) << " bits "
		<< skips << " skips" << std::endl;

	cra.result(res);
	if (res != actual) {
		report << res << " != " << actual << std::endl;
		report << "pprod: " << pprod << "\n" << "pprod / actual: " << (pprod / actual) << "\n";
		report << " *** CRABuilderProbSingle failed. ***" << std::endl;
		return EXIT_FAILURE ;
	}

	for (size_t i = 0 ; i < Size ; ++i){
		ModularField F(primes[i]);
		Element tmp1,tmp2 ;
		F.init(tmp1,res);
		F.init(tmp2,residues[i]);
		if(!F.areEqual(tmp1,tmp2)){
			report << tmp1 << "!=" << tmp2 << std::endl;
			report << " *** CRABuilderProbSingle failed. ***" << std::endl;
			return EXIT_FAILURE ;
		}
	}

	report << "CRABuilderProbSingle exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}

// testing CRABuilderFullSingle
template< class T >
int test_full_single(std::ostream & report, size_t PrimeSize, size_t Size)
{
	// true result
        size_t maxbits = (PrimeSize-1) * Size;
	size_t resbits = 1 + (random() % maxbits);
	Integer actual = Integer::random(resbits);

	typedef Givaro::Modular<double> ModularField ;

        PrimeIterator<IteratorCategories::DeterministicTag> pgen(PrimeSize);

	report << "CRABuilderFullSingle (" << maxbits << ")";
	report << " actual length " << actual.bitsize() << std::endl;
	CRABuilderFullSingle<ModularField> cra(maxbits) ;
	Integer res = 0; // the result
	T residue;
	T prime;
	{ /* init */
		prime = *pgen;
		residue = actual % prime;
		call_initialize(cra, prime, residue);
		++pgen;
	}
	size_t itercount = 1;
	while (!cra.terminated())
	{ /* progress */
		prime = *pgen;
		residue = actual % prime;
		if (cra.noncoprime((integer)prime)) {
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "got a duplicate prime from deterministic prime gen" << std::endl;
			Integer mod;
			report << prime << ' ' << cra.getModulus(mod) << std::endl;
			report << " *** CRABuilderFullSingle failed. ***" << std::endl;
			return EXIT_FAILURE ;
		}
		call_progress(cra, prime, residue);
		++itercount;
		++pgen;
	}
	report << "  " << itercount << " iterations, " << itercount*(PrimeSize-1) << " bits "
		<< std::endl;

	cra.result(res);
	if (res != actual) {
		report << res << " != " << actual << std::endl;
		report << " *** CRABuilderFullSingle failed. ***" << std::endl;
		return EXIT_FAILURE ;
	}

	report << "CRABuilderFullSingle exiting successfully." << std::endl;
	return EXIT_SUCCESS ;
}

// testing CRABuilderEarlyMultip
template< class T >
int test_early_multip(std::ostream & report, size_t PrimeSize, size_t Taille, size_t Size)
{

	typedef typename std::vector<T>                     Vect ;
	typedef typename Vect::iterator                  Iterator;
	typedef typename std::vector<Vect>::iterator VectIterator;
	typedef Givaro::Modular<double>             ModularField ;
	typedef ModularField::Element                     Element;
	typedef std::vector<Integer>                      IntVect;
	typedef std::vector<Element>                        pVect;

	/*  primes */
	Vect primes(Size) ;
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
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

	report << "EarlyMultpCRA (" <<  LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD << ')' << std::endl;
	CRABuilderEarlyMultip<ModularField> cra( LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD ) ;
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
			report << "CRABuilderEarlyMultip exiting successfully." << std::endl;
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
				report << " *** CRABuilderEarlyMultip failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "CRABuilderEarlyMultip exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}


#if 1 /* testing CRABuilderFullMultipMatrix */
template< class T>
int test_full_multip_matrix(std::ostream & report, size_t PrimeSize,
			    size_t Size, std::pair<size_t, size_t> dims)
{

	typedef typename Givaro::ZRing<T>        Unparam ;
	typedef typename std::vector<T>                         Vect ;
	typedef typename LinBox::BlasMatrix<Unparam>          Matrix ;
	typedef typename std::vector<Matrix>                 MatVect ;
	typedef typename Vect::iterator                      Iterator;
	typedef typename MatVect::iterator                MatIterator;
	typedef typename LinBox::BlasMatrix<Givaro::ZRing<Integer> >  IntMatrix ;

	typedef Givaro::Modular<double>                        Field ;
	typedef Field::Element                                Element;
	typedef typename LinBox::BlasMatrix<Field>           pMatrix ;

	Givaro::ZRing<Integer> Z ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
		++RP ;
	}

	/*  residues */
	Unparam U ;
	const Matrix Zero (U,dims.first,dims.second);
	MatVect residues(Size,Zero) ;
	for (size_t k = 0 ; k < Size ; ++k) {
		residues[k] = Zero ;
		for (size_t i = 0 ; i < Zero.rowdim() ; ++i)
			for (size_t j = 0 ; j < Zero.coldim() ; ++j)
				residues[k].setEntry( i,j,Integer::random(PrimeSize-1) );
	}

	Iterator  genprime =   primes.begin()  ; // prime iterator
	MatIterator residu = residues.begin()  ; // residu iterator

	double LogIntSize = (double)(PrimeSize+1)*std::log(2.)+std::log((double)Size)+1 ;

	std::pair<size_t,double> my_pair(dims.first*dims.second,LogIntSize)  ;

	report << "CRABuilderFullMultipMatrix (" <<  my_pair.first << ", " << my_pair.second << ')' << std::endl;
	CRABuilderFullMultipMatrix<Field> cra( my_pair ) ;
	IntMatrix result(Z,dims.first,dims.second); // the result
	{ /* init */
		Field F(*genprime);
		pMatrix residue(F,dims.first,dims.second) ; // temporary
		//! @bug it is not possible to allocate some memory and use submatrices ?
		for (size_t i = 0 ; i < residue.rowdim(); ++i)
			for (size_t j = 0 ; j < residue.coldim(); ++j)
				F.init(residue.refEntry(i,j),(*residu).getEntry(i,j));
		cra.initialize(F,residue);
		++genprime;
		++residu;
	}
	while (genprime < primes.end() /*  && !cra.terminated() */ )
	{ /* progress */
		if (cra.noncoprime((integer)*genprime)) {
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "CRABuilderFullMultipMatrix exiting successfully." << std::endl;
			return EXIT_SUCCESS ; // pas la faute à cra...
		}
		Field F(*genprime);
		pMatrix residue(F,dims.first,dims.second) ; // temporary
		for (size_t i = 0 ; i < residue.rowdim(); ++i)
			for (size_t j = 0 ; j < residue.coldim(); ++j)
				F.init(residue.refEntry(i,j),(*residu).getEntry(i,j));

		cra.progress(F,residue);
		++genprime;
		++residu ;
	}

	cra.result(result);

	for (size_t i = 0 ; i < Size ; ++i){
		Field F(primes[i]);
		for (size_t l = 0 ; l < dims.first ; ++l)
			for (size_t m = 0 ; m < dims.second ; ++m) {
				Element tmp1,tmp2 ;
				F.init(tmp1,result.getEntry(l,m));
				F.init(tmp2,residues[i].getEntry(l,m));
				if(!F.areEqual(tmp1,tmp2)){
					report << result.getEntry(l,m) << ';' << residues[i].getEntry(l,m) << '@' << primes[i] << std::endl;
					report << i << ':' << l << ',' << m << "> " << tmp1 << "!=" << tmp2 << std::endl;
					report << " *** CRABuilderFullMultipMatrix failed. ***" << std::endl;
					return EXIT_FAILURE ;
				}
			}
	}

	report << "CRABuilderFullMultipMatrix exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}
#endif

// testing CRABuilderFullMultip
template< class T>
int test_full_multip(std::ostream & report, size_t PrimeSize, size_t Size, size_t Taille)
{

	typedef typename std::vector<T>                    Vect ;
	typedef typename std::vector<Vect>             VectVect ;
	typedef std::vector<Integer>                    IntVect ;
	typedef typename Vect::iterator                 Iterator;
	typedef typename VectVect::iterator         VectIterator;

	typedef Givaro::Modular<double >           ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename std::vector<Element>             pVect ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
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

	double LogIntSize = (double)PrimeSize*std::log(2.)+std::log((double)Size)+1 ;

	report << "CRABuilderFullMultip (" <<  LogIntSize << ')' << std::endl;
	CRABuilderFullMultip<ModularField> cra( LogIntSize ) ;
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
			report << "CRABuilderFullMultip exiting successfully." << std::endl;
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
				report << " *** CRABuilderFullMultip failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "CRABuilderFullMultip exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}

// testing RationalCRABuilderFullMultip
template< class T>
int test_full_multip_rat(std::ostream & report, size_t PrimeSize, size_t Size, size_t Taille)
{
	typedef typename std::vector<T>                    Vect ;
	typedef std::vector<Integer>                    IntVect ;

	typedef Givaro::Modular<double >           ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename std::vector<Element>             pVect ;

    static_assert(std::is_same<T,Element>::value, "Can only test modular<double> for now");

    /* true answer */
    size_t bitlen = (PrimeSize-1) * Size / 2;
    Integer act_den = Integer::random(bitlen);
    IntVect act_num(Taille);
    for (auto& num_elt : act_num) {
        num_elt = Integer::random<false>(bitlen);
    }

	Vect primes ;
    std::vector<ModularField> fields;
    Vect denom_imgs;
	/*  probably not all coprime... */
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
        fields.emplace_back(*RP);
        denom_imgs.emplace_back();
        while (fields.back().isZero(fields.back().init(denom_imgs.back(), act_den))) {
            // hit a denominator divisor, try again
            ++RP;
            fields.pop_back();
            fields.emplace_back(*RP);
        }
        primes.emplace_back(*RP);
        ++RP;
	}

	/*  residues */
    std::vector<pVect> residues;
    for (size_t i=0; i < Size; ++i) {
        residues.emplace_back(Taille);
        for (size_t j=0; j < Taille; ++j) {
            fields[i].init(residues.back()[j], act_num[j]);
            fields[i].divin(residues.back()[j], denom_imgs[i]);
        }
	}

	auto genprime = primes.begin()  ; // prime iterator
    auto field_it = fields.begin();
	auto residu = residues.begin()  ; // residu iterator

	double LogIntSize = (double)PrimeSize*std::log(2.)+std::log((double)Size)+1 ;

	report << "RationalCRABuilderFullMultip (" <<  LogIntSize << ')' << std::endl;
	RationalCRABuilderFullMultip<ModularField> cra( LogIntSize ) ;
	IntVect res_num(Taille) ; // the result
    Integer res_den;
	{ /* init */
		cra.initialize(*field_it, *residu);
		++genprime;
        ++field_it;
		++residu;
	}
	while (genprime != primes.end() /* && !cra.terminated()*/ )
	{ /* progress */
		if (cra.noncoprime((integer)*genprime))
		{
			report << "bad luck, you picked twice the same prime..." <<std::endl;
			report << "RationalCRABuilderFullMultip exiting successfully." << std::endl;
			return EXIT_SUCCESS ; // pas la faute à cra...
		}
		cra.progress(*field_it,*residu);
		++genprime;
        ++field_it;
		++residu ;
	}

	cra.result(res_num, res_den);

    if (act_den % res_den != 0) {
        report << " *** RationalCRABuilderFullMultip failed. ***" << std::endl;
        report << "denominator mismatch: " << res_den << " != " << act_den << std::endl;
        return EXIT_FAILURE ;
    }
    Integer dm = act_den / res_den;

    for (size_t i = 0; i < Taille; ++i) {
        if (act_num[i] != res_num[i]*dm) {
            report << " *** RationalCRABuilderFullMultip failed. ***" << std::endl;
            report << "numerator mismatch: " << res_num[i] << " != " << act_num[i] << std::endl;
            return EXIT_FAILURE ;
        }
    }

	for (size_t i = 0 ; i < Size ; ++i){
		for (size_t j = 0 ; j < Taille ; ++j) {
			Element tmp1, tmp2 ;
			fields[i].init(tmp1,res_num[j]);
            fields[i].init(tmp2, res_den);
            fields[i].mulin(tmp2, residues[i][j]);
			if(!fields[i].areEqual(tmp1,tmp2)){
				report << " *** RationalCRABuilderFullMultip failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "RationalCRABuilderFullMultip exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}



#if 1 /* testing CRABuilderFullMultipFixed */
template< class T>
int test_full_multip_fixed(std::ostream & report, size_t PrimeSize, size_t Size, size_t Taille)
{


	typedef typename std::vector<T>                    Vect ;
	typedef typename std::vector<Vect>             VectVect ;
	typedef std::vector<Integer>                    IntVect ;
	typedef IntVect::iterator                IntVectIterator;
	typedef typename Vect::iterator                 Iterator;
	typedef typename VectVect::iterator         VectIterator;

	typedef Givaro::Modular<double >           ModularField ;
	typedef ModularField::Element                    Element;
	typedef typename std::vector<Element>             pVect ;
	typedef typename pVect::iterator          pVectIterator ;

	Vect primes(Size) ;
	/*  probably not all coprime... */
	PrimeIterator<IteratorCategories::HeuristicTag> RP((unsigned )PrimeSize);
	for (size_t i = 0 ; i < Size ; ++i) {
		primes[i] = *RP;
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

	double LogIntSize = (double)PrimeSize*std::log(2.)+std::log((double)Size) ;

	std::pair<size_t,double> my_pair(Taille,LogIntSize)  ;

	report << "CRABuilderFullMultipFixed (" <<  my_pair.first << ", " << my_pair.second << ')' << std::endl;

	CRABuilderFullMultipFixed<ModularField> cra( my_pair ) ;
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
			report << "CRABuilderFullMultipFixed exiting successfully." << std::endl;
			return EXIT_SUCCESS ; // pas la faute à cra...
		}
		ModularField F(*genprime);
		pVectIterator residue_jt = residue.begin();
		for (size_t i = 0 ; i < Taille; ++i)
			F.init(residue[i],(*residu)[i]);
		cra.progress(F,residue_jt);
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
				report << " *** CRABuilderFullMultipFixed failed. ***" << std::endl;
				return EXIT_FAILURE ;
			}
		}
	}

	report << "CRABuilderFullMultipFixed exiting successfully." << std::endl;

	return EXIT_SUCCESS ;
}
#endif

bool test_CRA_algos(size_t PrimeSize, size_t Size, size_t Taille, size_t iters)
{
	bool pass = true ;
	std::ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT,
							   INTERNAL_DESCRIPTION);



	typedef std::pair<size_t,size_t> Pair ;

	/* EARLY SINGLE */
	_LB_REPEAT( if (test_early_single<double>(report,22,Size))                       pass = false ;  ) ;
	_LB_REPEAT( if (test_early_single<integer>(report,PrimeSize,Size))               pass = false ;  ) ;

	/* PROB SINGLE */
	_LB_REPEAT( if (test_prob_single<double>(report,22,Size))                       pass = false ;  ) ;
	_LB_REPEAT( if (test_prob_single<integer>(report,PrimeSize,Size))               pass = false ;  ) ;

	/* FULL SINGLE */
	_LB_REPEAT( if (test_full_single<double>(report,22,Size))                       pass = false ;  ) ;
	_LB_REPEAT( if (test_full_single<integer>(report,PrimeSize,Size))               pass = false ;  ) ;

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

    /* FULL MULTIPLE RATIONAL */
	_LB_REPEAT( if (test_full_multip_rat<double>(report,22,Size,Taille))                 pass = false ;  ) ;
	_LB_REPEAT( if (test_full_multip_rat<double>(report,22,Size,Taille/4))                 pass = false ;  ) ;

	return pass ;

}

#include "test-common.h"
#include "linbox/util/timer.h"

// launching tests
int main(int ac, char ** av)
{

	/*  Argument parsing/setting */

	static size_t       n = 30;    /*  Taille */
	static size_t       p = 22;    /*  PrimeSize */
	// static size_t    seed =  0;    /*  ! unused */
	static size_t   iters = 2;    /* _LB_REPEAT */

        static Argument as[] = {
                { 'n', "-n N", "Set number of primes.", TYPE_INT , &n },
                { 'p', "-p P", "Set size of test primes.", TYPE_INT , &p },
                { 'i', "-i I", "Perform each test for I iterations.",     TYPE_INT, &iters },
		END_OF_ARGUMENTS
        };

	parseArguments (ac, av, as);

	bool pass ;

	srand((unsigned)time(NULL));             // seeding
	size_t PrimeSize   =  p;       // size of the residues/primes
	size_t Size        =  n ;      // nb of residues/primes
	size_t Taille      =  2*Size ; // nb of vectors of residues


	commentator().start("CRA-Algos test suite", "CRA-Algos");

	pass = test_CRA_algos(PrimeSize,Size,Taille,iters) ;

	commentator().stop(MSG_STATUS (pass), "CRA-Algos test suite");
	return !pass ;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
