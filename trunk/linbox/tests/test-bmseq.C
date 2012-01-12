/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-dense.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by George Yuhasz <yuhasz@gmail.com>
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <set>
#include <list>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/bm-seq.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;
using namespace std;

/* Test 1: Identity matrix in dense representation
 *
 * Construct a dense representation of an n x n identity matrix and check
 * whether the output of its application to a series of random vectors is equal
 * to the input.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentity (Field &F, long n, int iterations)
{
	typedef typename Vector<Field>::Dense Vector;
	typedef BlasMatrix <Field> Blackbox;

	commentator.start ("Testing identity apply", "testIdentity", iterations);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = true;
	bool iter_passed = true;

	int i, j;

	Blackbox I (F, n, n);
	typename Field::Element one;

	F.init (one, 1);

	for (i = 0; i < n; i++)
		I.setEntry (i, i, one);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		report << "Input vector: ";
		printVector<Field> (F, report, v);

		I.apply (w, v);

		report << "Output vector: ";
		printVector<Field> (F, report, w);

		for (j = 0; j < n; j++)
			if (!F.areEqual (w[j], v[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentity");

	return ret;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 2; // was 100
	//static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.start("Dense matrix black box test suite", "BlasMatrix");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);



	Field::Element one, zero;
	F.init(one,1);
	F.init(zero,0);
	BlasMatrix<Field> D(F,2,2);
	BlasMatrix<Field> zero24(F,2,4);
	for(size_t i=0; i<2; i++)
		D.setEntry(i,i,one);
	D.setEntry(1,0,one);
	BM_Seq<Field> seq(2,D);
	BlasMatrix<Field> S2(F,2,2);
	MatrixDomain<Field> MD(F);
	BM_Seq<Field>::BM_iterator bmit(seq, 0), bmit2(seq.BM_begin());
	bmit.setDelta(4);
	BM_Seq<Field>::BM_iterator::TerminationState check = bmit.state();
	while(!check.IsGeneratorFound() && !check.IsSequenceExceeded()){
					bmit++;
					check = bmit.state();
	}
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	MD.add(S2,D,D);
	seq.push_back(S2);
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	MD.addin(S2,D);
	seq.push_back(S2);
	bmit++;
	check = bmit.state();
	if(check.IsSequenceExceeded())
					report << "Sequence Exceeded" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit->begin(); it != bmit->end(); it++)
					(*it).write(report);
	if(check.IsGeneratorFound())
					report << "Generator Found" << endl;
	report << "mu = " << bmit.get_mu() << endl;
	report << "sigma = " << bmit.get_sigma() << endl;
	report << "beta = " << bmit.get_beta() << endl;
	BM_Seq<Field>::BM_iterator::TerminationState check2 = bmit2.state();
	while(!check2.IsGeneratorFound() && !check2.IsSequenceExceeded()){
					++bmit2;
					check2 = bmit2.state();
	}
	if(bmit==bmit2)
					report << "Iterators are equal" << endl;
	if(bmit2==seq.BM_end())
					report << "bmit2 is equal to end" << endl;
	for(list<BlasMatrix<Field> >::iterator it = bmit2->begin(); it != bmit2->end(); it++)
					(*it).write(report);
	BM_Seq<Field>::BM_iterator bmit3 = seq.BM_begin();
	bmit3 = bmit;
	if(bmit==bmit3)
					report << "Iterators are equal" << endl;
	vector<BlasMatrix<Field> >gen(bmit.GetGenerator());
	int d = bmit.get_mu();
	for(int j = 0; j <= d; j++)
					gen[j].write(report);

	commentator.stop("dense matrix black box test suite");
	return pass ? 0 : -1;
}
