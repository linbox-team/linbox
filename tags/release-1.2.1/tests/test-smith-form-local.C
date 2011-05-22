/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* tests/test-local-smith.C
 * Copyright (C) LinBox
 *
 * Written by David Saunders
 *
 * --------------------------------------------------------
 * See COPYING for license information
 */


/*! @file  tests/test-smith-form-local.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */



#include "linbox/linbox-config.h"

#include <vector>
#include <functional>

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/field/PIR-modular-int32.h"
//#include "linbox/field/PIR-modular-double.h"
#include "linbox/field/local2_32.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/smith-form-local.h"
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/timer.h>

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
class foobar {
public:
	typedef typename LocalPIR::Element first_argument_type;
	typedef LocalPIR second_argument_type;
	typedef void result_type;
	void operator()(typename LocalPIR::Element& d, const LocalPIR& R) const
	{
		typename LocalPIR::Element x; R.init(x, 2);  R.mulin(x, d);
		if (R.isUnit(d)) R.divin(d, d);
		else R.gcd(d, d, x);
	}
};

template<>
class foobar<LinBox::Local2_32> {
public:
	typedef LinBox::Local2_32 LocalPIR;

	typedef LocalPIR::Element first_argument_type;
	typedef LocalPIR second_argument_type;
	typedef void result_type;
	void operator()(LocalPIR::Element& d, const LocalPIR& R) const
	{
		if(d != 0)    {
			int r = 1;
			while ( !(d & 1) ) {
				d >>= 1;
				r <<= 1;
			}
			d = r;
		}
	}
};

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

#if 0
template<>
class pplt<LinBox::NTL_PID_zz_p> {
public:
	typedef LinBox::NTL_PID_zz_p LocalPIR;

	pplt(LocalPIR R) : _R_(R){}
	bool operator() (LocalPIR::Element a, LocalPIR::Element b)
	{
		if ( b == 0 ) return true;
		else if ( a == 0 ) return false;
		else return NTL::rep(a) <= NTL::rep(b);
	}
	//protected:
	LocalPIR _R_;
};
#endif

template <class LocalPIR>
static bool testLocalSmith (const LocalPIR &R, vector<typename LocalPIR::Element>& d, string s)
{
	typedef typename LocalPIR::Element Elt;
	typedef DenseMatrix<LocalPIR> Blackbox;

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
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
		{ D[i][i] = d[i]; L[i][i]=U[i][i]=1; }
	for (i = 0; i < n; ++ i)
		for ( j = 0; j < i; ++ j) {
			D[i][j] = D[j][i] = 0;
			L[i][j] = rand() % 10;
			L[j][i] = 0;
			U[j][i] = rand() % 10;
			U[i][j] = 0;
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
	for_each(d.begin(), d.end(), bind2nd(foobar<LocalPIR>(), R));
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
		commentator.progress();
	}
	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true, pass1 = true;

	static size_t n = 6;
	static integer q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator.start("Local Smith Form test suite", "LocalSmith");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

#if 0
//PIRModular does not currently support the local ring interface -bds Mar2011
  { // first local ring type
	typedef PIRModular<int32_t> LocalPID;
	LocalPID R (536870912); // 2^32
	//typedef PIRModular<dense> LocalPID;
	//LocalPID R (32768);
	vector<LocalPID::Element> d(n);

	commentator.start ("Testing local smith on singular dense mat over PIRModular", "testSingular");
	for( size_t i = 0; i < n; ++i ) d[i] = i;
	if (!testLocalSmith<LocalPID> (R, d, "PIRModular<int32_t>")) pass1 = false;
	commentator.stop ("testSingular");

	commentator.start ("Testing local smith on nonsingular dense mat over PIRModular", "testNonsingular");
	for( size_t i = 0; i < n; ++i ) d[i] = i+1;
	if (!testLocalSmith<LocalPID> (R, d, "PIRModular<int32_t>")) pass1 = false;
	commentator.stop ("testNonsingular");
  }
  if (not pass1) report << "PIRModular FAIL" << std::endl;
#endif

  { // second local ring type
	typedef Local2_32 LocalPID;
	LocalPID R;
	vector<LocalPID::Element> d(n);

	commentator.start ("Testing local smith on singular dense mat over Local2_32", "testSingular");
	for( size_t i = 0; i < n; ++i )
		d[i] = (LocalPID::Element) i;
	if (!testLocalSmith<LocalPID> (R, d, "Local2_32")) pass = false;
	commentator.stop ("testSingular");

	commentator.start ("Testing local smith on nonsingular dense mat over Local2_32", "testNonsingular");
	for( size_t i = 0; i < n; ++i )
		d[i] = (LocalPID::Element) i+1;
	if (!testLocalSmith<LocalPID> (R, d, "Local2_32")) pass = false;
	commentator.stop ("testNonsingular");
  }
  if (not pass) report << "PIRModular FAIL" << std::endl;

	commentator.stop("Local Smith Form test suite");
	return pass and pass1 ? 0 : -1;
}

