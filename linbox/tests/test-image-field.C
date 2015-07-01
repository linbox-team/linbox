/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-hom.C
 *
 * Written by Dave Saunders <saunders@cis.udel.edu>
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "linbox/field/modular.h"
#include "linbox/field/image-field.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 65521U;
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 1000000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus (default 65521)", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test (default 1000000)", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test (default 100)", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test (default 1)", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << endl << "ImageField test suite" << endl;
	cout.flush ();
	bool pass = true;

	// need to test generic def, test nonsense attempts, etc.

	Modular<uint32> F_uint32 ((uint32) q);
	Modular<uint16> F_uint16 ((uint16) q);
	ImageField<Modular<uint16>, Modular<uint32> > IF32(F_uint16, F_uint32);

	uint16 a=2, b;
	uint32 z=2, w;
	IF32.image(w, a);
	pass = pass && IF32.areEqual(z, w);
	IF32.preimage(b, z);
	pass = pass && F_uint16.areEqual(a, b);

	uint32 x, y;
	IF32.smul(x, 2, 3);
	IF32.mul(y, 2, 3);
	pass = pass && IF32.areEqual(x, y);
	IF32.smulin(x, 5);
	IF32.mulin(y, 5);
	pass = pass && IF32.areEqual(x, y);

	IF32.saxpy(z, 7, x, 11);
	IF32.axpy(w, 7, x, 11);
	pass = pass && IF32.areEqual(z, w);
	IF32.saxpyin(z, 7, x);
	IF32.axpyin(w, 7, x);
	pass = pass && IF32.areEqual(z, w);
	/*
	*/

	cout << endl << "ImageField " << (pass ? "pass" : "FAIL") << endl;
	cout.flush ();
	return pass ? 0 : -1;
}
