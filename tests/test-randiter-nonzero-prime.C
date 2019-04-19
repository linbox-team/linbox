/* tests/test-randiter-nonzero-prime.C
 * Copyright (C) 2001-2019 the LinBox group
 *
 * ------------------------------------
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
 *.
 */


/*! @file  tests/test-randiter-nonzero-prime.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */



#include "linbox/linbox-config.h"
#include <givaro/givranditer.h>

#include <iostream>
#include <fstream>


#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

/* Test 1: Nonzero random elements
 *
 * Construct a sequence of random elements from the nonzero random iterator and
 * check that they are all nonzero
 *
 * F - Field over which to perform computations
 * iterations - Number of iterations to perform
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNonzeroRandom (Field &F, unsigned int iterations)
{
	int i;

	commentator().start ("Testing nonzero random elements", "testNonzeroRandom", iterations);

	bool ret = true;

	typename Field::Element x;
	typename Field::RandIter r (F);

	//NonzeroRandIter <Field, typename Field::RandIter> rp (F, r);
	Givaro::GeneralRingNonZeroRandIter <Field> rp (r);

	for (i = 0; i <(int) iterations; i++) {
		commentator().startIteration ((unsigned int)i);

		rp.random (x);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Element: ";
		F.write (report, x);
		report << endl;

		if (F.isZero (x)) {
			ret = false;
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Element is zero" << endl;
		}

		commentator().progress ();
		commentator().stop ("done");
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testNonzeroRandom");

	return ret;
}


#include <linbox/randiter/random-prime.h>

template <class PGenerator>
bool testPrimeIterator(PGenerator& genprime, unsigned int iterations)
{
	commentator().start ("Testing prime random elements", "testPrimeRandom", iterations);

    Givaro::IntPrimeDom IPD;
    bool pass(true);

    for(size_t i=0; i<iterations; ++i, ++genprime) {
		commentator().startIteration ((unsigned int)i);
        if (!IPD.isprime(*genprime, _GIVARO_ISPRIMETESTS_<<3)) {
            std::cerr << "***** ERROR ***** Iteration: " << i << ", supposed prime: " << *genprime << std::endl;
            pass = false;
        }
		commentator().progress ();
		commentator().stop ("done");
    }
	commentator().stop (MSG_STATUS (pass), (const char *) 0, "testPrimeRandom");
    return pass;
}

template <class PGenerator>
bool testMask(PGenerator& genprime, unsigned int iterations, size_t numproc, size_t maxprocs) 
{
	commentator().start ("Testing masked prime random elements", "testMaskedPrimeRandom", iterations);

    bool pass(true);
    const uint64_t MAXMASK = ( (UINT64_C(1)<<genprime.getShift())-1);
    const uint64_t EXPECTS = (numproc<<1) | 0x1;
    if (EXPECTS != genprime.getMask()) pass = false;

    for(size_t i=0; i<iterations; ++i, ++genprime) {
		commentator().startIteration ((unsigned int)i);
        if ( (*genprime & MAXMASK) != EXPECTS ) {
            std::cerr << "***** ERROR ***** Iteration: " << i << ", supposed masked prime: " << *genprime << ", with low bits: " << (*genprime & MAXMASK) << ", when expected: " << EXPECTS << std::endl;
            pass = false;
        }
		commentator().progress ();
		commentator().stop ("done");
    }
	commentator().stop (MSG_STATUS (pass), (const char *) 0, "testMaskedPrimeRandom");
    return pass;
}

template <class MPGenerator>
bool testMaskedPrimeIterators(size_t maxprocs, size_t s, unsigned int iterations)
{
    bool pass(true);
    for(size_t numproc=0; numproc<maxprocs; ++numproc) {
        MPGenerator genprime(numproc,maxprocs,s);
        pass &= testPrimeIterator(genprime,iterations);
        pass &= testMask(genprime,iterations,numproc,maxprocs);
        MPGenerator genprime2(numproc,maxprocs,s);
        pass &= testMask(genprime2,iterations,numproc,maxprocs);
    }

    return pass;
}


template <class PGenerator>
bool testPrimeIterators(size_t s, unsigned int iterations)
{
    bool pass(true);
    PGenerator genprime(s);
    pass &= testPrimeIterator(genprime,iterations);

    return pass;
}




int main (int argc, char **argv)
{
	bool pass = true;

	static integer q = 101;
	static unsigned int iterations = 1000;
	static size_t size= 23;
	static size_t maxprocs= 10;
    static int rseed = (int)time(NULL);

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'r', "-r R", "bitsize is R.", TYPE_INT,     &size },
        { 's', "-s S", "Random generator seed.", TYPE_INT,     &rseed }	,
		{ 'm', "-m M", "Maximum number of masked generators.", TYPE_INT,     &maxprocs },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	Givaro::Modular<uint32_t> F (q);

	std::srand(rseed);

	commentator().start("Nonzero&Prime random iterator test suite", "Givaro::GeneralRingNonZeroRandIter");

	if (!testNonzeroRandom (F, iterations)) pass = false;

	pass &= testPrimeIterators< PrimeIterator<IteratorCategories::HeuristicTag> > (size, iterations);
	pass &= testPrimeIterators< PrimeIterator<IteratorCategories::DeterministicTag> > (size, iterations);
	pass &= testPrimeIterators< PrimeIterator<IteratorCategories::UniformTag> > (size, iterations);
	pass &= testMaskedPrimeIterators< MaskedPrimeIterator<IteratorCategories::HeuristicTag> > (maxprocs, size, iterations);
	pass &= testMaskedPrimeIterators< MaskedPrimeIterator<IteratorCategories::DeterministicTag> > (maxprocs, size, iterations);
	pass &= testMaskedPrimeIterators< MaskedPrimeIterator<IteratorCategories::UniformTag> > (maxprocs, size, iterations);

	commentator().stop("Nonzero&Prime random iterator test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
