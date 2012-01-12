/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/test-hom.C
 * Copyright (C) LinBox
 *
 * Written by Dave Saunders <saunders@cis.udel.edu>
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */


/*! @file  tests/test-image-field.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

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
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16_t modulus.", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	cout << endl << "ImageField test suite" << endl;
	cout.flush ();
	bool pass = true;

	// need to test generic def, test nonsense attempts, etc.

	Modular<uint32_t> F_uint32_t ((uint32_t) q);
	Modular<uint16_t> F_uint16_t ((uint16_t) q);
	ImageField<Modular<uint16_t>, Modular<uint32_t> > IF32(F_uint16_t, F_uint32_t);

	uint16_t a=2, b;
	uint32_t z=2, w;
	IF32.image(w, a);
	pass = pass && IF32.areEqual(z, w);
	IF32.preimage(b, z);
	pass = pass && F_uint16_t.areEqual(a, b);

	uint32_t x, y;
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
