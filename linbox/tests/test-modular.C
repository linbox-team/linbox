/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-modular.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Modular interace.
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/modular.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Random number test
 *
 * Test that the random iterator over the given field works
 */

template <class Field>
bool testRandomIterator (const Field &F,
			 const char *text,
			 unsigned int num_trials,
			 unsigned int num_categories,
			 unsigned int hist_len) 
{
	std::ostringstream str;

	str << "Testing " << text << "::RandIter" << ends;

	commentator.start (str.str ().c_str (), "testRandomIterator");
	ostream &report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);

	bool ret = true;

	integer card;
	unsigned int i;
	std::vector<int> categories1 (num_categories, 0);
	std::vector<int> categories2 (num_categories, 0);
	std::list<std::vector<int> > diff_categories;
	std::list<typename Field::Element> x_queue;

	F.cardinality (card);

	typename Field::RandIter iter (F);
	typename Field::Element x, x_prev, x_prev2, d;

	std::list<std::vector<int> >::iterator diff_cat_iter;

	for (i = 0; i < hist_len; ++i)
		diff_categories.push_back (std::vector<int> (num_categories, 0));

	// I make the simplifying assumption that field elements are actually
	// C++ ints. Otherwise, I don't know how to place the numbers into
	// categories in any well-defined manner.
	for (i = 0; i < num_trials; ++i) {
		F.assign (x_prev2, x_prev);
		F.assign (x_prev, x);

		iter.random (x);

		categories1[x % num_categories]++;
		categories2[(unsigned int) (double (x) / double (card) * num_categories)]++;

		typename std::list<typename Field::Element>::iterator x_queue_iter = x_queue.begin ();
		diff_cat_iter = diff_categories.begin ();

		for (; x_queue_iter != x_queue.end (); ++x_queue_iter, ++diff_cat_iter) {
			F.sub (d, *x_queue_iter, x);
			(*diff_cat_iter)[d % num_categories]++;
		}

		x_queue.push_front (x);

		while (x_queue.size () > hist_len)
			x_queue.pop_back ();
	}

	double p, chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories1[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

	report << "Test of distribution uniformity (low-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (low-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Random iterator's values do not appear to be uniformly distributed"
			<< std::endl;
		ret = false;
	}

	chi_squared = 0.0;

	for (i = 0; i < num_categories; ++i)
		chi_squared += pow (double (categories2[i]) -
				    double (num_trials) / double (num_categories), 2);

	p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

	report << "Test of distribution uniformity (high-order): chi^2 = "
	       << chi_squared * num_categories / num_trials << std::endl;
	report << "Test of distribution uniformity (high-order):     p = " << p << std::endl;

	if (p < 0.05 || p > 0.95) {
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Random iterator's values do not appear to be uniformly distributed"
			<< std::endl;
		ret = false;
	}

	diff_cat_iter = diff_categories.begin ();

	int idx = 0;

	for (; diff_cat_iter != diff_categories.end (); ++diff_cat_iter, ++idx) {
		chi_squared = 0.0;

		for (i = 0; i < num_categories; ++i)
			chi_squared += pow (double ((*diff_cat_iter)[i]) -
					    double (num_trials) / double (num_categories), 2);

		p = chiSquaredCDF (chi_squared * num_categories / num_trials, num_categories - 1);

		report << "Test of " << idx + 1 << " spacing: chi^2 = "
		       << chi_squared * num_categories / num_trials << std::endl;
		report << "Test of " << idx + 1 << " spacing:     p = " << p << std::endl;

		if (p < 0.05 || p > 0.95) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Difference values do not appear to be uniformly distributed"
				<< std::endl;
			ret = false;
		}
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomIterator");

	return ret;
}

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 1000000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus (default 18446744073709551557)", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus (default 2147483647)", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus (default 65521)", TYPE_INTEGER, &q3 },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test (default 1000000)", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test (default 100)", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test (default 1)", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "Modular field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	Modular<integer> F_integer (q1);
	Modular<uint32> F_uint32 ((uint32) q2);
	Modular<uint16> F_uint16 ((uint16) q3);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F_integer, "Modular<integer>", iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint32,  "Modular<uint32>",  iterations, n, false)) pass = false;
	if (!runFieldTests (F_uint16,  "Modular<uint16>",  iterations, n, false)) pass = false;

	if (!testRandomIterator (F_integer, "Modular<integer>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint32,  "Modular<uint32>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F_uint16,  "Modular<uint16>", trials, categories, hist_level)) pass = false;

#if 0
	FieldArchetype K(new LargeModular(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of Modular field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
