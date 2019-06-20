/* tests/test-det.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */


/*! @file  tests/test-rational-reconstruction-base.C
 * @ingroup tests
 * @ingroup CRA
 * @brief  tests rational reconstruction using rational-cra-var-prec.h .
 * @test no doc.
 */



#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/linbox-config.h"

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
//#include "linbox/field/gmp-rational.h"

#include "linbox/algorithms/rational-cra-var-prec.h"
#include "linbox/algorithms/cra-builder-var-prec-early-single.h"
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/algorithms/classic-rational-reconstruction.h"
#include "linbox/algorithms/fast-rational-reconstruction.h"
#include "linbox/randiter/random-prime.h"

#include "test-common.h"

using namespace LinBox;

/* Test: Rational reconstruction of random fraction using rational-cra-var-prec.h
 *
 * Constructs a random pair of numerator/denominator
 * Reconstructs it based on rational reconstruction
 *
 * n - size of numerator
 * d - size of denominator
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

struct ModularFraction {
	integer a_,b_;
	ModularFraction(const integer& a, const integer& b) :
		a_(a), b_(b)
	{}
	template<typename Field>
	typename Field::Element& operator()(typename Field::Element& d, const Field& F) const
	{
		F.init(d,a_);
		F.divin(d,b_);
		return d;
	}
};

static bool testRandomFraction (size_t n, size_t d, int iterations)
{
	commentator().start ("Testing rational reconstruction on random fractions", "testRandFrac", (unsigned int)iterations);

	bool ret = true;
	// bool done;
	int i;
	// size_t j;
	// size_t k;
	integer num,den;

	for (i = 0; i < iterations; i++) {
		commentator().startIteration ((unsigned int)i);

		integer::nonzerorandom(num, n);
		integer::nonzerorandom(den, d);
		integer g; gcd(g,num, den);
		num /= g; den /= g;
		if (i %2 ) integer::negin(num);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "True fraction: " << num << " / " << den;
		report << endl;

		ModularFraction iteration(num,den);
                typedef Givaro::ModularBalanced<double> Field;
		PrimeIterator<IteratorCategories::HeuristicTag> genprime(FieldTraits<Field>::bestBitSize(n));

		Givaro::ZRing<Integer> Z;
		ClassicRationalReconstruction<Givaro::ZRing<Integer> > RRB1(Z,false,false);
		ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > RRB2(Z,false,false);

		integer a1_1, b1_1, a2_1, b2_1, a3_1, b3_1, a4_1, b4_1;
		integer a1_2, b1_2, a2_2, b2_2, a3_2, b3_2, a4_2, b4_2;
		integer a1_3, b1_3, a2_3, b2_3, a3_3, b3_3, a4_3, b4_3;
		integer a1_4, b1_4, a2_4, b2_4, a3_4, b3_4, a4_4, b4_4;

		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > RR1_1(Z,INCREMENTAL,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > > cra1_1(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR1_1);
		cra1_1(a1_1,b1_1,iteration,genprime);
		if ((a1_1 != num)  || (b1_1 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, incremental, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR2_1(Z,INCREMENTAL,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra2_1(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR2_1);
		cra2_1(a2_1,b2_1,iteration,genprime);
		if ((a2_1 != num)  || (b2_1 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, incremental, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > RR3_1(RRB1,INCREMENTAL,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > > cra3_1(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR3_1);
		cra3_1(a3_1,b3_1,iteration,genprime);
		if ((a3_1 != num)  || (b3_1 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, incremental, classic) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR4_1(RRB2,INCREMENTAL,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra4_1(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR4_1);
		cra4_1(a4_1,b4_1,iteration,genprime);
		if ((a4_1 != num)  || (b4_1 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, incremental, classic) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > RR1_2(Z,QUADRATIC,0,10);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > > cra1_2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR1_2);
		cra1_2(a1_2,b1_2,iteration,genprime);
		if ((a1_2 != num)  || (b1_2 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, quadratic, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR2_2(Z,QUADRATIC,0,10);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra2_2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR2_2);
		cra2_2(a2_2,b2_2,iteration,genprime);
		if ((a2_2 != num)  || (b2_2 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, quadratic, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > RR3_2(RRB1,QUADRATIC,0,10);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > > cra3_2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR3_2);
		cra3_2(a3_2,b3_2,iteration,genprime);
		if ((a3_2 != num)  || (b3_2 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, quadratic, classic) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR4_2(RRB2,QUADRATIC,0,10);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra4_2(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR4_2);
		cra4_2(a4_2,b4_2,iteration,genprime);
		if ((a4_2 != num)  || (b4_2 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, quadratic, classic) failed" << endl;
		}


		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > RR1_3(Z,GEOMETRIC,0,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > > cra1_3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR1_3);
		cra1_3(a1_3,b1_3,iteration,genprime);
		if ((a1_3 != num)  || (b1_3 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, geometric, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR2_3(Z,GEOMETRIC,0,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra2_3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR2_3);
		cra2_3(a2_3,b2_3,iteration,genprime);
		if ((a2_3 != num)  || (b2_3 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, geometric, fast) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > RR3_3(RRB1,GEOMETRIC,0,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > > cra3_3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR3_3);
		cra3_3(a3_3,b3_3,iteration,genprime);
		if ((a3_3 != num)  || (b3_3 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, geometric, classic) failed" << endl;
		}

		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR4_3(RRB2,GEOMETRIC,0,5);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > > cra4_3(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR4_3);
		cra4_3(a4_3,b4_3,iteration,genprime);
		if ((a4_3 != num)  || (b4_3 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (MaxQ, geometric, classic) failed" << endl;
		}

		size_t H = (n > d) ? n : d;
		H /= (size_t)d;
		++H;
		H *=2;
		++H;
		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > RR1_4(Z,CERTIFIED,0,H);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, FastRationalReconstruction<Givaro::ZRing<Integer> > > > cra1_4(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR1_4);
		cra1_4(a1_4,b1_4,iteration,genprime);
		if ((a1_4 != num)  || (b1_4 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, certified, fast) failed" << endl;
		}

		//RReconstruction<Givaro::ZRing<Integer>, FastMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR2_4(Z,CERTIFIED);
		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > RR3_4(RRB1,CERTIFIED,0,H);
		RationalChineseRemainderVarPrec< CRABuilderVarPrecEarlySingle< Field >,
		RReconstruction<Givaro::ZRing<Integer>, ClassicRationalReconstruction<Givaro::ZRing<Integer> > > > cra3_4(LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD, RR3_4);
		cra3_4(a3_4,b3_4,iteration,genprime);
		if ((a3_4 != num)  || (b3_4 != den) ) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: rational reconstruction (Wang, certified, classic) failed" << endl;
		}

		//RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > RR4_4(RRB2,CERTIFIED);

		commentator().stop ("done");
		commentator().progress ();
		//commentator().progress (i, iterations);
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRationalDeterminantGeneric");

	return ret;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 5;
	//static integer q = 4093U;
	static int iterations = 2;

	static Argument args[] = {
		{ 'n', "-n N", "Set size of test numerator/denominator to N", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	//Givaro::Modular<int> F (q);

	commentator().start("Rational reconstruction test suite", "rr");

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!testRandomFraction          (n, n,iterations)) pass = false;

	commentator().stop("Rational reconstruction test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
