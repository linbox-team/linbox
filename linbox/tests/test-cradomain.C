/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2010 LinBox
 *
 * Time-stamp: <05 Apr 11 11:01:44 Jean-Guillaume.Dumas@imag.fr>
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file tests/test-cradomain.C
 * @ingroup tests
 * @brief tests LinBox::ChineseRemainer
 * @test tests LinBox::ChineseRemainer (see \ref CRA)
 */

#include "linbox/algorithms/cra-domain.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/cra-full-multip-fixed.h"
#include "linbox/algorithms/cra-givrnsfixed.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/integer.h"

using namespace LinBox;

struct Interator {
	std::vector<integer> _v;
	double maxsize;

	Interator(const std::vector<integer>& v) :
		_v(v), maxsize(0.0)
	{
		for(std::vector<integer>::const_iterator it=_v.begin();
		    it != _v.end(); ++it) {
			double ds = ::Givaro::naturallog(*it);
			maxsize = (maxsize<ds?ds:maxsize);
		}
	}

	Interator(int n, int s) :
		_v(n), maxsize(0.0)
	{
		for(std::vector<integer>::iterator it=_v.begin();
		    it != _v.end(); ++it) {
			Integer::random<false>(*it, s);
			double ds = ::Givaro::naturallog(*it);
			maxsize = (maxsize<ds?ds:maxsize);
		}
	}

	const std::vector<integer>& getVector()
	{
	       	return _v;
       	}
	double getLogSize() const
	{
	       	return maxsize;
	}

	template<typename Field>
	std::vector<typename Field::Element>& operator()(std::vector<typename Field::Element>& v,
							 const Field& F) const
	{
		v.resize(_v.size());
		std::vector<integer>::const_iterator vit=_v.begin();
		typename std::vector<typename Field::Element>::iterator eit=v.begin();
		for( ; vit != _v.end(); ++vit, ++eit){
			F.init(*eit, *vit);
		}

		return v;
	}
};

struct InteratorIt;
namespace LinBox
{
	template<class Element> struct CRATemporaryVectorTrait<InteratorIt , Element> {
		typedef typename std::vector<double>::iterator Type_t;
	};
}

struct InteratorIt : public Interator {

	mutable std::vector<double> _C;

	InteratorIt(const std::vector<integer>& v) :
		Interator(v), _C(v.size())
	{}
	InteratorIt(int n, int s) :
		Interator(n,s), _C(n)
	{}

	template<typename Iterator, typename Field>
	Iterator& operator()(Iterator& res, const Field& F) const
	{
		std::vector<integer>::const_iterator vit=this->_v.begin();
		std::vector<double>::iterator eit=_C.begin();
		for( ; vit != _v.end(); ++vit, ++eit) {
			F.init(*eit, *vit);
		}

		return res=_C.begin();
	}

};

template<typename Field> struct InteratorBlas;
namespace LinBox
{
	template<class Element,class Field> struct CRATemporaryVectorTrait<InteratorBlas<Field> , Element> {
		typedef typename LinBox::BlasMatrix<Element>::pointer Type_t;
	};
}

template<typename Field>
struct InteratorBlas : public Interator {
	typedef typename Field::Element Element;
	typedef LinBox::BlasMatrix<Element> Matrix;
	typedef typename Matrix::pointer Pointer;
	mutable Matrix _C;

	InteratorBlas(const std::vector<integer>& v) : Interator(v), _C((int)v.size(), (int)1) {}
	InteratorBlas(int n, int s) : Interator(n,s), _C(n,1) {}

	Pointer& operator()(Pointer& res, const Field& F) const
	{
		std::vector<integer>::const_iterator vit=this->_v.begin();
		res = _C.getWritePointer();
		for( ; vit != _v.end(); ++vit, ++res)
			F.init(*res, *vit);

		return res=_C.getWritePointer();
	}

};

#include <typeinfo>


template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRA(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
	std::vector<integer> Res(N);
	cra( Res, iteration, genprime);
	bool locpass = std::equal( Res.begin(), Res.end(), iteration.getVector().begin() );
	if (locpass) report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	else {
		report << "***ERROR***: ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
		std::vector<integer>::const_iterator Rit=Res.begin();
		std::vector<integer>::const_iterator Oit=iteration.getVector().begin();
		for( ; Rit!=Res.end(); ++Rit, ++Oit)
			if (*Rit != *Oit)
				report << *Rit <<  " != " << * Oit << std::endl;

	}
	return locpass;
}

template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRAbegin(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
	std::vector<integer> Res(N);
	std::vector<integer>::iterator ResIT= Res.begin();
	cra( ResIT, iteration, genprime);
	bool locpass = std::equal( Res.begin(), Res.end(), iteration.getVector().begin() );
	if (locpass) report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	else {
		report << "***ERROR***: ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
		std::vector<integer>::const_iterator Rit=Res.begin();
		std::vector<integer>::const_iterator Oit=iteration.getVector().begin();
		for( ; Rit!=Res.end(); ++Rit, ++Oit)
			if (*Rit != *Oit)
				report << *Rit <<  " != " << * Oit << std::endl;

	}
	return locpass;
}

template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRAWritePointer(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
	LinBox::BlasMatrix<integer> Res( (int)N, (int)N);
	cra( Res.getWritePointer(), iteration, genprime);
	bool locpass = std::equal( iteration.getVector().begin(), iteration.getVector().end(), Res.getWritePointer() );

	if (locpass) report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	else {
		report << "***ERROR***: ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
	}
	return locpass;
}


bool TestCra(int N, int S, size_t seed)
{

	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT,
							   INTERNAL_DESCRIPTION);
	// std::ostream &report = std::cout;

	size_t new_seed = (seed?(seed):(BaseTimer::seed())) ;
	report << "TestCra(" << N << ',' << S << ',' << new_seed << ')' << std::endl;
	Integer::seeding(new_seed);

	Interator iteration(N, S);
	InteratorIt iterationIt(iteration.getVector());
	InteratorBlas<LinBox::Modular<double> > iterationBlas(iteration.getVector());
	LinBox::RandomPrimeIterator genprime( 24, new_seed );

	bool pass = true;

	pass &= TestOneCRA< LinBox::EarlyMultipCRA< LinBox::Modular<double> >,
	     Interator, LinBox::RandomPrimeIterator>(
						     report, iteration, genprime, N, 5);

	pass &= TestOneCRA< LinBox::EarlyMultipCRA< LinBox::Modular<double> >,
	     Interator, LinBox::RandomPrimeIterator>(
						     report, iteration, genprime, N, 15);

	pass &= TestOneCRA< LinBox::FullMultipCRA< LinBox::Modular<double> >,
	     Interator, LinBox::RandomPrimeIterator>(
						     report, iteration, genprime, N, iteration.getLogSize()+1);

	pass &= TestOneCRA< LinBox::FullMultipCRA< LinBox::Modular<double> >,
	     Interator, LinBox::RandomPrimeIterator>(
						     report, iteration, genprime, N, 3*iteration.getLogSize()+15);

	pass &= TestOneCRAbegin<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorIt, LinBox::RandomPrimeIterator>(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,iteration.getLogSize()+1));

	pass &= TestOneCRAbegin<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorIt, LinBox::RandomPrimeIterator>(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,3*iteration.getLogSize()+15));


	pass &= TestOneCRAWritePointer<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorIt, LinBox::RandomPrimeIterator>(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,iterationIt.getLogSize()+1) );

	pass &= TestOneCRAWritePointer<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorIt, LinBox::RandomPrimeIterator>(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,3*iterationIt.getLogSize()+15) );

	pass &= TestOneCRAWritePointer<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorBlas< LinBox::Modular<double> >,
	     LinBox::RandomPrimeIterator>(
					  report, iterationBlas, genprime, N, std::pair<size_t,double>(N,iterationIt.getLogSize()+1) );

	pass &= TestOneCRAWritePointer<LinBox::FullMultipFixedCRA< LinBox::Modular<double> >,
	     InteratorBlas< LinBox::Modular<double> >,
	     LinBox::RandomPrimeIterator>(
					  report, iterationBlas, genprime, N, std::pair<size_t,double>(N,3*iterationIt.getLogSize()+15) );


        std::vector<integer> PrimeSet;
        double PrimeSize = 0.0;
        for( ; PrimeSize < (iterationIt.getLogSize()+1); ++genprime ) {
            if (find(PrimeSet.begin(), PrimeSet.end(), *genprime) == PrimeSet.end()) {
                PrimeSet.push_back( *genprime );
                PrimeSize += ::Givaro::naturallog(*genprime);
            }
        }

        std::vector<integer>::iterator psit = PrimeSet.begin();

	pass &= TestOneCRA<
            LinBox::GivaroRnsFixedCRA< LinBox::Modular<double> >,
            Interator,
            std::vector<integer>::iterator,
            std::vector<integer> >(
                 report, iteration, psit, N, PrimeSet);


	if (pass) report << "TestCra(" << N << ',' << S << ')' << ", passed." << std::endl;
	else
		report << "***ERROR***: TestCra(" << N << ',' << S << ')' << " ***ERROR***" << std::endl;

	return pass;
}

#include "test-common.h"
#include "linbox/util/timer.h"

int main (int argc, char **argv)
{
	static size_t n = 10;
	static size_t s = 30;
	static size_t seed = 0;
	static int iterations = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT , &n },
		{ 's', "-s S", "Set size of test integers.", TYPE_INT , &s },
		{ 'z', "-z Z", "Set seed.", TYPE_INT , &seed },
		{ 'i', "-i I", "Perform each test for I iterations.",     TYPE_INT, &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	LinBox::commentator.start("CRA-Domain test suite", "CRADom");
	bool pass = true;

	for(int i=0; pass && i<iterations; ++i)
		pass &= TestCra((int)n,(int)s,seed);

	LinBox::commentator.stop(MSG_STATUS (pass), (const char *) 0,"CRA-Domain test suite");
	return pass ? 0 : -1;
}
