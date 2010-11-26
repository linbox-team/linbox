/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* tests/test-det.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/linbox-config.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/gmp-rational.h"

#include "linbox/algorithms/rational-cra2.h"
#include "linbox/algorithms/varprec-cra-early-single.h"
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/algorithms/classic-rational-reconstruction.h"
#include "linbox/algorithms/fast-rational-reconstruction.h"
#include "linbox/randiter/random-prime.h"

#include "test-common.h"

using namespace LinBox;

/* Test: Rational reconstruction of random fraction using rational-cra2.h
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
	ModularFraction(const integer& a, const integer& b): a_(a), b_(b) {}
	template<typename Field>
	typename Field::Element& operator()(typename Field::Element& d, const Field& F) const {
		F.init(d,a_);
		F.divin(d,b_);
		return d;
	}
};

static bool testRandomFraction (size_t n, size_t d, int iterations) 
{
	commentator.start ("Testing rational reconstruction on random fractions", "testRandFrac", iterations);

	bool ret = true;
	bool done;
	int i;
	size_t j, k;
	integer num,den;

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		integer::nonzerorandom(num, n);
		integer::nonzerorandom(den, d);
		integer g; gcd(g,num, den);
		num /= g; den /= g;
		if (i %2 ) integer::negin(num);
		
		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "True fraction: " << num << " / " << den; 
		report << endl;

		ModularFraction iteration(num,den);
		int bits = 26;
		RandomPrimeIterator genprime( bits);

		PID_integer Z;
		ClassicRationalReconstruction<PID_integer> RRB1(Z,false,false);
		ClassicMaxQRationalReconstruction<PID_integer> RRB2(Z,false,false);

		integer a1_1, b1_1, a2_1, b2_1, a3_1, b3_1, a4_1, b4_1;
		integer a1_2, b1_2, a2_2, b2_2, a3_2, b3_2, a4_2, b4_2;
		integer a1_3, b1_3, a2_3, b2_3, a3_3, b3_3, a4_3, b4_3;
		integer a1_4, b1_4, a2_4, b2_4, a3_4, b3_4, a4_4, b4_4;

		RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > RR1_1(Z,INCREMENTAL,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >, 
			RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > > cra1_1(4UL, RR1_1);
		cra1_1(a1_1,b1_1,iteration,genprime);
		if ((a1_1 != num)  || (b1_1 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, incremental, fast) failed" << endl;
		}

		RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > RR2_1(Z,INCREMENTAL,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > > cra2_1(4UL, RR2_1);
		cra2_1(a2_1,b2_1,iteration,genprime);
		if ((a2_1 != num)  || (b2_1 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, incremental, fast) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > RR3_1(RRB1,INCREMENTAL,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > > cra3_1(4UL, RR3_1);
		cra3_1(a3_1,b3_1,iteration,genprime);
		if ((a3_1 != num)  || (b3_1 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, incremental, classic) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR4_1(RRB2,INCREMENTAL,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > > cra4_1(4UL, RR4_1);
		cra4_1(a4_1,b4_1,iteration,genprime);
		if ((a4_1 != num)  || (b4_1 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, incremental, classic) failed" << endl;
		}

		RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > RR1_2(Z,QUADRATIC,0,10);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > > cra1_2(4UL, RR1_2);
		cra1_2(a1_2,b1_2,iteration,genprime);
		if ((a1_2 != num)  || (b1_2 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, quadratic, fast) failed" << endl;
		}

		RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > RR2_2(Z,QUADRATIC,0,10);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > > cra2_2(4UL, RR2_2);
		cra2_2(a2_2,b2_2,iteration,genprime);
		if ((a2_2 != num)  || (b2_2 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, quadratic, fast) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > RR3_2(RRB1,QUADRATIC,0,10);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > > cra3_2(4UL, RR3_2);
		cra3_2(a3_2,b3_2,iteration,genprime);
		if ((a3_2 != num)  || (b3_2 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, quadratic, classic) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR4_2(RRB2,QUADRATIC,0,10);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > > cra4_2(4UL, RR4_2);
		cra4_2(a4_2,b4_2,iteration,genprime);
		if ((a4_2 != num)  || (b4_2 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, quadratic, classic) failed" << endl;
		}


		RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > RR1_3(Z,GEOMETRIC,0,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > > cra1_3(4UL, RR1_3);
		cra1_3(a1_3,b1_3,iteration,genprime);
		if ((a1_3 != num)  || (b1_3 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, geometric, fast) failed" << endl;	
		}

		RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > RR2_3(Z,GEOMETRIC,0,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > > cra2_3(4UL, RR2_3);
		cra2_3(a2_3,b2_3,iteration,genprime);
		if ((a2_3 != num)  || (b2_3 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, geometric, fast) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > RR3_3(RRB1,GEOMETRIC,0,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > > cra3_3(4UL, RR3_3);
		cra3_3(a3_3,b3_3,iteration,genprime);
		if ((a3_3 != num)  || (b3_3 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, geometric, classic) failed" << endl;
		}

		RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR4_3(RRB2,GEOMETRIC,0,5);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > > cra4_3(4UL, RR4_3);
		cra4_3(a4_3,b4_3,iteration,genprime);
		if ((a4_3 != num)  || (b4_3 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (MaxQ, geometric, classic) failed" << endl;
		}

		size_t H = (n > d) ? n : d;
		H /= bits;
		++H; 
		H *=2;
		++H;
		RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > RR1_4(Z,CERTIFIED,0,H);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, FastRationalReconstruction<PID_integer> > > cra1_4(4UL, RR1_4);
		cra1_4(a1_4,b1_4,iteration,genprime);
		if ((a1_4 != num)  || (b1_4 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, certified, fast) failed" << endl;
		}

//RReconstruction<PID_integer, FastMaxQRationalReconstruction<PID_integer> > RR2_4(Z,CERTIFIED);
		RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > RR3_4(RRB1,CERTIFIED,0,H);
		RationalRemainder2< VarPrecEarlySingleCRA< Modular<double> >,
			RReconstruction<PID_integer, ClassicRationalReconstruction<PID_integer> > > cra3_4(4UL, RR3_4);
		cra3_4(a3_4,b3_4,iteration,genprime);
		if ((a3_4 != num)  || (b3_4 != den) ) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rational reconstruction (Wang, certified, classic) failed" << endl;
}

//RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > RR4_4(RRB2,CERTIFIED);

		commentator.stop ("done");
	 	commentator.progress ();
	 	//commentator.progress (i, iterations);
 	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRationalDeterminantGeneric");

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
		{ '\0' }
	};

	parseArguments (argc, argv, args);
	//Modular<int> F (q);

	commentator.start("Rational reconstruction test suite", "rr"); 

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!testRandomFraction          (n, n,iterations)) pass = false;

	commentator.stop("Rational reconstruction test suite");
	return pass ? 0 : -1;
}
