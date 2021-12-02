/* Copyright (C) 2010 LinBox
 *
 * Time-stamp: <27 Aug 20 17:09:01 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file tests/test-cradomain.C
 * @ingroup tests
 * @brief tests LinBox::ChineseRemainer
 * @test tests LinBox::ChineseRemainer (see \ref CRA)
 */

#include "linbox/ring/modular.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-builder-early-multip.h"
#include "linbox/algorithms/cra-builder-full-multip.h"
#include "linbox/algorithms/cra-builder-full-multip-fixed.h"
#include "linbox/algorithms/cra-givrnsfixed.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/integer.h"

using namespace LinBox;

template <class IntVect_t = BlasVector<Givaro::ZRing<Integer>>>
struct Interator {
    using IntVect = IntVect_t;
    using Field = typename IntVect_t::Field;

    Field _f;
	IntVect _v;
	double maxsize;

	Interator(const IntVect& v) :
            _f(v.field()), _v(v), maxsize(0.0)
	{
		for(auto it=_v.begin(); it != _v.end(); ++it) {
			//!@bug bb: *it < 0 ?
			double ds = Givaro::naturallog(*it);
			maxsize = (maxsize<ds?ds:maxsize);
		}
	}

	Interator(int n, int s) :
            _f(), _v(_f,n), maxsize(0.0)
	{
		for(auto it=_v.begin(); it != _v.end(); ++it) {
			Integer::random<false>(*it, s);
			double ds = Givaro::naturallog(*it);
			maxsize = (maxsize<ds?ds:maxsize);
		}
	}

	const IntVect& getVector() const
	{
		return _v;
    }
	double getLogSize() const
	{
		return maxsize;
	}

	template<typename Vect, typename Field>
	IterationResult operator()(Vect& v, const Field& F) const
	{
		v.resize(_v.size());
		auto vit=_v.begin();
		auto eit=v.begin();
		for( ; vit != _v.end(); ++vit, ++eit){
			F.init(*eit, *vit);
		}

		return IterationResult::CONTINUE;
	}
};

template <class IntVect = BlasVector<Givaro::ZRing<Integer>>>
struct InteratorIt : public Interator<IntVect> {

	// could use BlasVector and changeField
	mutable std::vector<double> _vectC;

	InteratorIt(const IntVect& v) :
		Interator<IntVect>(v), _vectC(v.size())
	{}
	InteratorIt(int n, int s) :
		Interator<IntVect>(n,s), _vectC((size_t)n)
	{}

	template<typename Iterator, typename Field>
	IterationResult operator()(Iterator& res, const Field& F) const
	{
		auto vit=this->_v.begin();
		auto eit=_vectC.begin();
		for( ; vit != this->_v.end(); ++vit, ++eit) {
			F.init(*eit, *vit);
		}

		res=_vectC.begin();
		return IterationResult::CONTINUE;
	}

};

template<typename Field, class IntVect = BlasVector<Givaro::ZRing<Integer>>>
struct InteratorBlas : public Interator<IntVect> {
	typedef typename Field::Element Element;
	typedef LinBox::BlasMatrix<Givaro::ZRing<Element> > Matrix;
	typedef typename Matrix::Element_ptr Pointer;
	typename Givaro::ZRing<Element> _field;
	mutable Matrix _vectC;

	InteratorBlas(const IntVect& v) :
		Interator<IntVect>(v),
		_field(),
		_vectC(_field,(int)v.size(), (int)1)
	{}

	InteratorBlas(int n, int s) :
		Interator<IntVect>(n,s),
		_field(),
		_vectC(_field,n,1) {}

	IterationResult operator()(Pointer& res, const Field& F) const
	{
		auto vit=this->_v.begin();
		res = _vectC.getPointer();
		for( ; vit != this->_v.end(); ++vit, ++res)
			F.init(*res, *vit);

		res=_vectC.getPointer();
		return IterationResult::CONTINUE;
	}

};

#include <typeinfo>


template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRA(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
    typename Iter::IntVect Res( typename Iter::Field(), N);
	cra( Res, iteration, genprime);


    Integer base; cra.getModulus(base);
    auto Riter(Res.begin());
    auto Iiter(iteration.getVector().begin());
    bool locpass=true;
    for( ; Riter != Res.end(); ++Riter, ++Iiter) {
        locpass &= isZero( ( *Riter - *Iiter ) % base );
    }

	if (locpass) report << "ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	else {
		report << "***ERROR***: ChineseRemainder<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
		auto Rit=Res.begin();
		auto Oit=iteration.getVector().begin();
		for( ; Rit!=Res.end(); ++Rit, ++Oit)
			if (*Rit != *Oit) {
				report << *Rit <<  " != " << *Oit << std::endl;
            }
	}
	return locpass;
}

template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRAbegin(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	Givaro::ZRing<Integer> Z;
	report << "ChineseRemainderBeg<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
	BlasVector<Givaro::ZRing<Integer> > Res(Z,N);
	BlasVector<Givaro::ZRing<Integer> >::iterator ResIT= Res.begin();
	cra( ResIT, iteration, genprime);

    Integer base; cra.getModulus(base);
    auto Riter(Res.begin()) ;
    auto Iiter(iteration.getVector().begin());
    bool locpass=true;
    for( ; Riter != Res.end(); ++Riter, ++Iiter) {
        locpass &= isZero( ( *Riter - *Iiter ) % base );
    }


//	bool locpass = std::equal( Res.begin(), Res.end(), iteration.getVector().begin() );
	if (locpass) report << "ChineseRemainderBeg<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	else {
		report << "***ERROR***: ChineseRemainderBeg<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
		BlasVector<Givaro::ZRing<Integer> >::const_iterator Rit=Res.begin();
		BlasVector<Givaro::ZRing<Integer> >::const_iterator Oit=iteration.getVector().begin();
		for( ; Rit!=Res.end(); ++Rit, ++Oit)
			if (*Rit != *Oit)
				report << *Rit <<  " != " << * Oit << std::endl;

	}
	return locpass;
}

template<typename Builder, typename Iter, typename RandGen, typename BoundType>
bool TestOneCRAWritePointer(std::ostream& report, Iter& iteration, RandGen& genprime, size_t N, const BoundType& bound)
{
	report << "ChineseRemainderWP<" << typeid(Builder).name() << ">(" << bound << ')' << std::endl;
	LinBox::ChineseRemainder< Builder > cra( bound );
	Givaro::ZRing<Integer> Z ;
	LinBox::BlasMatrix<Givaro::ZRing<Integer> > Res(Z, (int)N, (int)N);
	cra( Res.getPointer(), iteration, genprime);

    Integer base; cra.getModulus(base);
    auto Riter(Res.getPointer());
    auto Iiter(iteration.getVector().begin());
    bool locpass=true;
    for( ; Iiter != iteration.getVector().end(); ++Riter, ++Iiter) {
        locpass &= isZero( ( *Riter - *Iiter ) % base );
    }

//	bool locpass = std::equal( iteration.getVector().begin(), iteration.getVector().end(), Res.getPointer() );

	if (locpass) {
		report << "ChineseRemainderWP<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << ", passed."  << std::endl;
	}
	else {
		report << "***ERROR***: ChineseRemainderWP<" << typeid(Builder).name() << ">(" << iteration.getLogSize() << ')' << "***ERROR***"  << std::endl;
	}
	return locpass;
}


bool TestCra(size_t N, int S, size_t seed)
{

	std::ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT,
							   INTERNAL_DESCRIPTION);
	// std::ostream &report = std::cout;
	Givaro::ZRing<Integer> Z;

	size_t new_seed = (seed?(seed):((size_t)BaseTimer::seed())) ;
	report << "TestCra(" << N << ',' << S << ',' << new_seed << ')' << std::endl;
	Integer::seeding(new_seed);

    // either of these should work
    using IntVect = BlasVector<Givaro::ZRing<Integer>>;
    // using IntVect = std::vector<Integer>;

	Interator<IntVect> iteration((int)N, S);
	InteratorIt<IntVect> iterationIt(iteration.getVector());
        typedef Givaro::ModularBalanced<double> Field;
	InteratorBlas<Field, IntVect> iterationBlas(iteration.getVector());
	PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(N), new_seed );

	bool pass = true;

	pass &= TestOneCRA< LinBox::CRABuilderEarlyMultip< Field > >(
						     report, iteration, genprime, N, 5);

	pass &= TestOneCRA< LinBox::CRABuilderEarlyMultip< Field > >(
						     report, iteration, genprime, N, 15);

	pass &= TestOneCRA< LinBox::CRABuilderFullMultip< Field > >(
						     report, iteration, genprime, N, iteration.getLogSize()+1);

	pass &= TestOneCRA< LinBox::CRABuilderFullMultip< Field > >(
						     report, iteration, genprime, N, 3*iteration.getLogSize()+15);

#if 0
	pass &= TestOneCRAbegin<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorIt, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,iteration.getLogSize()+1));

	pass &= TestOneCRAbegin<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorIt, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,3*iteration.getLogSize()+15));


	pass &= TestOneCRAWritePointer<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorIt, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,iterationIt.getLogSize()+1) );

	pass &= TestOneCRAWritePointer<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorIt, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
						       report, iterationIt, genprime, N, std::pair<size_t,double>(N,3*iterationIt.getLogSize()+15) );

	pass &= TestOneCRAWritePointer<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorBlas< Field >,
	     LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
					  report, iterationBlas, genprime, N, std::pair<size_t,double>(N,iterationIt.getLogSize()+1) );

	pass &= TestOneCRAWritePointer<LinBox::CRABuilderFullMultipFixed< Field >,
	     InteratorBlas< Field >,
	     LinBox::PrimeIterator<IteratorCategories::HeuristicTag> >(
					  report, iterationBlas, genprime, N, std::pair<size_t,double>(N,3*iterationIt.getLogSize()+15) );
#endif


	// XXX fixed prime set doesn't work with openmp version
#ifndef LINBOX_USES_OPENMP
        IntVect PrimeSet(Z,0);
        double PrimeSize = 0.0;
        for( ; PrimeSize < (iterationIt.getLogSize()+1); ++genprime ) {
            if (std::find(PrimeSet.begin(), PrimeSet.end(), *genprime) == PrimeSet.end()) {
                PrimeSet.push_back( *genprime );
                PrimeSize += Givaro::naturallog(*genprime);
            }
        }

	auto psseq = create_prime_sequence(PrimeSet);

	pass &= TestOneCRA< LinBox::GivaroRnsFixedCRA< Field > >(
                 report, iteration, psseq, N, PrimeSet);
#endif


	if (pass) report << "TestCra(" << N << ',' << S << ')' << ", passed." << std::endl;
	else
		report << "***ERROR***: TestCra(" << N << ',' << S << ')' << " ***ERROR***" << std::endl;

	return pass;
}

#include "test-common.h"
#include "linbox/util/timer.h"

int main (int argc, char **argv)
{

//     commentator().setMaxDetailLevel (1);
//     commentator().setMaxDepth (-1);
//     commentator().setReportStream (std::clog);

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

	LinBox::commentator().start("CRA-Domain test suite", "CRADom");
	bool pass = true;

	for(int i=0; pass && i<iterations; ++i)
		pass &= TestCra((size_t)n,(int)s,seed);

	LinBox::commentator().stop(MSG_STATUS (pass), "CRA-Domain test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
