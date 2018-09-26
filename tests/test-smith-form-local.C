/* tests/test-local-smith.C
 * Copyright (C) LinBox
 *
 * Written by David Saunders
 *
 * --------------------------------------------------------
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


/*! @file  tests/test-smith-form-local.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */

#include "linbox/linbox-config.h"


#include <functional>

#include "test-common.h"
#include "test-field.h"

#include "linbox/util/commentator.h"
#include "linbox/ring/local-pir-modular.h"
#include "linbox/ring/pir-modular-int32.h"
#include "linbox/ring/local2_32.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/timer.h"

using namespace LinBox;

/** @brief Test 1: Invariant factors of random dense matrices.
 *
 * Construct a random matrix which is equivalent to a random diagonal matrix,
 * and check its Smith form.
 *
 * R - PIR over which to perform computations
 * stream - Stream that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class LocalPIR>
typename LocalPIR::Element &
normal(typename LocalPIR::Element &d, LocalPIR & R)
{	typename LocalPIR::Element x; R.init(x, 2);  R.mulin(x, d);
	return R.gcdin(d, x);
}

template <class LocalPIR>
class pplt { // prime power less than
public:
	pplt(LocalPIR R) : _R_(R){}
	bool operator() (typename LocalPIR::Element a, typename LocalPIR::Element b)
	{
		if ( b == 0 ) return true;
		else if ( a == 0 ) return false;
		else return a <= b;
	}
	//protected:
	LocalPIR _R_;
};

template <class LocalPIR>
static bool testLocalSmith (const LocalPIR &R, vector<typename LocalPIR::Element>& d, string s)
{
	typedef typename LocalPIR::Element Elt;
	typedef BlasMatrix<LocalPIR> Blackbox;

	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << s << endl;

	MatrixDomain<LocalPIR> MR (R);
	VectorDomain<LocalPIR> VD (R);

	bool ret = true;
	size_t i,j;
	size_t n = d.size();

	report << "Input vector:  ";
	VD.write (report, d);
	report << endl;

	// set up A equiv diag d.
	Blackbox L (R, n, n), D (R, n, n), U (R, n, n), A (R, n, n);
	for( i = 0; i < n; ++i )
		{ R.assign(D[i][i], d[i]); L[i][i]=U[i][i]=R.one; }
	for (i = 0; i < n; ++ i)
		for ( j = 0; j < i; ++ j) {
			R.assign(D[i][j], R.assign(D[j][i], 0));
			R.assign(L[i][j], std::rand() % 10);
			R.assign(L[j][i], 0);
			R.assign(U[j][i], std::rand() % 10);
			R.assign(U[i][j], 0);
		}
		
	MR.mul(A,L,D);
	MR.mulin(A,U);

	list< Elt > Inv;
	SmithFormLocal< LocalPIR > SmithForm;
	//timer.start();
	SmithForm( Inv, A, R );
	//timer.stop();
	//report << "Time " << timer <<"\n"; report.flush();

	report << "Computed invariants: ";
	report << "[";
	typedef typename list<Elt>::iterator listptr;
	for (listptr p = Inv.begin(); p != Inv.end(); ++p)
		report << *p << ", ";
	//report << "\b\b]" << endl;
	report << "normalize done" << endl; report.flush();

	// figure true invariants
	pplt<LocalPIR> lt(R);
	for (size_t i = 0; i < d.size(); ++i) normal(d[i], R);
	stable_sort(d.begin(), d.end(), lt);
	report << "True invariants: ";
	VD.write (report, d) << endl; report.flush();

	typename vector<Elt>::iterator q;
	listptr p;
	for (p = Inv.begin(), q = d.begin(); q != d.end(); ++p, ++q)
	{
		if ( !R.areEqual (*p, *q ) ) {
			report << "ERROR: Computed invariants incorrect" << endl;
			ret = false;
		}
		//commentator().progress();
	}
	return ret;
}


int main (int argc, char **argv) {
    bool pass1(true), pass2(true);
    static int64_t n = 6;
    static int  q = 3;
    static int32_t e = 4;
    static int rseed = (int)time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the ring Z/q^eZ.", TYPE_INT, &q },
        { 'e', "-e e", "Operate over the ring Z/q^eZ.", TYPE_INT, &e },
        { 's', "-s S", "Random generator seed.", TYPE_INT,     &rseed }	,
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	std::srand(rseed);
	//FFLAS::writeCommandString(std::cout << argv[0] << ' ', args) << std::endl;

	commentator().start("Local Smith Form test suite", "LocalSmith");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "q = " << q << std::endl;

#if 1
  { // zero-th local ring type: modulus p^e as integer.
	typedef LocalPIRModular<integer> Ring;
	Ring R (q, e);
	vector<Ring::Element> d(n);

	commentator().start ("Testing local smith over LocalPIRModular<integer>", "local smith");
	for( int32_t i = 0; i < n; ++i ) R.init(d[i],i);
	pass1 = testLocalSmith<Ring> (R, d, "LocalPIRModular<integer>");
	commentator().stop ("local smith");
	if (not pass1) report << "LocalPIRModular smith FAIL" << std::endl;

  }
#endif
  { // first local ring type: modulus p^e as int32_t.
	typedef PIRModular<int32_t> LocalPIR;
	//typedef PIRModular<dense> LocalPIR;
	//LocalPIR R (81); 
	LocalPIR R (q, e);
	vector<LocalPIR::Element> d(n);

	commentator().start ("Testing local smith on singular dense mat over PIRModular", "testSingular");
	for( int32_t i = 0; i < n; ++i ) R.init(d[i],i);
	pass1 = testLocalSmith<LocalPIR> (R, d, "PIRModular<int32_t>");
	commentator().stop ("testSingular");
	if (not pass1) report << "PIRModular sing FAIL" << std::endl;

	commentator().start ("Testing local smith on nonsingular dense mat over PIRModular", "testNonsingular");
	LocalPIR::RandIter r(R,rseed);
	LocalPIR::Element e; R.init(e);
	for( int32_t i = 0; i < n; ++i ) {	
		r.random(e);
		do { r.random(e); } while (R.isZero(e));
		R.assign(d[i], e);
	}
	bool p = testLocalSmith<LocalPIR> (R, d, "PIRModular<int32_t>");
	if (not p) report << "PIRModular nonsing FAIL" << std::endl;
	commentator().stop ("testNonsingular");
	pass1 = pass1 and p;
	if (not pass1) report << "PIRModular FAIL" << std::endl;
  }

  { // second local ring type: m = 2^32
	typedef Local2_32 LocalPIR;
	LocalPIR R;
	vector<LocalPIR::Element> d(n);

	commentator().start ("Testing local smith on singular dense mat over Local2_32", "testSingular");
	for( int64_t i = 0; i < n; ++i )
		d[i] = (LocalPIR::Element) i;
	if (!testLocalSmith<LocalPIR> (R, d, "Local2_32")) pass2 = false;
	commentator().stop ("testSingular");

	commentator().start ("Testing local smith on nonsingular dense mat over Local2_32", "testNonsingular");
	for( int64_t i = 0; i < n; ++i )
		d[i] = (LocalPIR::Element) i+1;
	if (!testLocalSmith<LocalPIR> (R, d, "Local2_32")) pass2 = false;
	commentator().stop ("testNonsingular");
	if (not pass2) report << "Local2_32 FAIL" << std::endl;
  }

	commentator().stop("Local Smith Form test suite");
	return pass1 and pass2 ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
