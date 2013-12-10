
/* tests/test-vector-domain.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

/*! @file  tests/test-vector-domain.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
 */


#include "linbox-config.h"

#include <iostream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/gf2.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/vector-domain-gf2.h"
#include "linbox/vector/stream.h"

#include "test-vector-domain.h"

using namespace std;
using namespace LinBox;

template <class Field>
bool testVectorDomain (const Field &F, const char *text, size_t n, unsigned int iterations)
{
	typedef BlasVector<Field> DenseVector;
	typedef BlasVector<Field> SparseSeqVector;
	typedef BlasVector<Field> SparseMapVector;
	typedef BlasVector<Field> SparseParVector;

	ostringstream str;
	str << "Testing VectorDomain <" << text << ">" << ends;
	commentator().start (str.str ().c_str ());

	bool pass = true;

	RandomDenseStream<Field, DenseVector> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomSparseStream<Field, SparseSeqVector> stream3 (F, 0.1, n, iterations), stream4 (F, 0.1, n, iterations);
	RandomSparseStream<Field, SparseMapVector> stream5 (F, 0.1, n, iterations), stream6 (F, 0.1, n, iterations);
	RandomSparseStream<Field, SparseParVector> stream7 (F, 0.1, n, iterations), stream8 (F, 0.1, n, iterations);

	if (!testDotProduct (F, "dense/dense", stream1, stream2)) pass = false;
	if (!testDotProduct (F, "sparse sequence/dense", stream3, stream1)) pass = false;
	if (!testDotProduct (F, "sparse associative/dense", stream5, stream1)) pass = false;
	if (!testDotProduct (F, "sparse parallel/dense", stream7, stream1)) pass = false;
	if (!testDotProduct (F, "sparse sequence/sparse sequence", stream3, stream4)) pass = false;
	if (!testDotProduct (F, "sparse associative/sparse sequence", stream5, stream3)) pass = false;
	if (!testDotProduct (F, "sparse parallel/sparse sequence", stream7, stream3)) pass = false;
	if (!testDotProduct (F, "sparse associative/sparse associative", stream5, stream6)) pass = false;
	if (!testDotProduct (F, "sparse parallel/sparse associative", stream7, stream6)) pass = false;
	if (!testDotProduct (F, "sparse parallel/sparse parallel", stream7, stream8)) pass = false;

	if (!testAddMul (F, "dense", stream1, stream2)) pass = false;
	if (!testAddMul (F, "sparse sequence", stream3, stream4)) pass = false;
	if (!testAddMul (F, "sparse associative", stream5, stream6)) pass = false;
	if (!testAddMul (F, "sparse parallel", stream7, stream8)) pass = false;

	if (!testSubMul (F, "dense", stream1, stream2)) pass = false;
	if (!testSubMul (F, "sparse sequence", stream3, stream4)) pass = false;
	if (!testSubMul (F, "sparse associative", stream5, stream6)) pass = false;
	if (!testSubMul (F, "sparse parallel", stream7, stream8)) pass = false;

	if (!testAXPY (F, "dense", stream1, stream2)) pass = false;
	if (!testAXPY (F, "sparse sequence", stream3, stream4)) pass = false;
	if (!testAXPY (F, "sparse associative", stream5, stream6)) pass = false;
	if (!testAXPY (F, "sparse parallel", stream7, stream8)) pass = false;

	if (!testCopyEqual (F, "dense/dense", stream1, stream1)) pass = false;
	if (!testCopyEqual (F, "dense/sparse sequence", stream1, stream3)) pass = false;
	if (!testCopyEqual (F, "dense/sparse associative", stream1, stream5)) pass = false;
	if (!testCopyEqual (F, "dense/sparse parallel", stream1, stream7)) pass = false;
	if (!testCopyEqual (F, "sparse sequence/dense", stream3, stream1)) pass = false;
	if (!testCopyEqual (F, "sparse sequence/sparse sequence", stream3, stream3)) pass = false;
	if (!testCopyEqual (F, "sparse sequence/sparse associative", stream3, stream5)) pass = false;
	if (!testCopyEqual (F, "sparse sequence/sparse parallel", stream3, stream7)) pass = false;
	if (!testCopyEqual (F, "sparse associative/dense", stream5, stream1)) pass = false;
	if (!testCopyEqual (F, "sparse associative/sparse sequence", stream5, stream3)) pass = false;
	if (!testCopyEqual (F, "sparse associative/sparse associative", stream5, stream5)) pass = false;
	if (!testCopyEqual (F, "sparse associative/sparse parallel", stream5, stream7)) pass = false;
	if (!testCopyEqual (F, "sparse parallel/dense", stream7, stream1)) pass = false;
	if (!testCopyEqual (F, "sparse parallel/sparse sequence", stream7, stream3)) pass = false;
	if (!testCopyEqual (F, "sparse parallel/sparse associative", stream7, stream5)) pass = false;
	if (!testCopyEqual (F, "sparse parallel/sparse parallel", stream7, stream8)) pass = false;

	commentator().stop (MSG_STATUS (pass));

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static unsigned int n = 100;
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static int q4 = 101;
	static unsigned int iterations = 2;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to N.", TYPE_INT,     &n },
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32_t modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16_t modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"field\" GF(P) [1] for uint8_t modulus.", TYPE_INTEGER, &q4 },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	Modular<integer> F_integer (q1);
	Modular<uint32_t> F_uint32_t ((uint32_t) q2);
	Modular<uint16_t> F_uint16_t ((uint16_t) q3);
	Modular<uint8_t> F_uint8_t ((uint8_t) q4);
	GF2 gf2(2);

	commentator().start("Vector domain test suite", "VectorDomain");

	commentator().setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	commentator().getMessageClass (TIMING_MEASURE).setMaxDepth (3);

	if (!testVectorDomain (F_integer, "Modular <integer>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint32_t, "Modular <uint32_t>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint16_t, "Modular <uint16_t>", n, iterations)) pass = false;
	if (!testVectorDomain (F_uint8_t, "Modular <uint8_t>", n, iterations)) pass = false;
//	if (!testVectorDomain (gf2, "GF2", n, iterations)) pass = false;

	commentator().stop("Vector domain test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
