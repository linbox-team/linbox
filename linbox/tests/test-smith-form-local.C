/* tests/test-local-smith.C
 * Copyright (C) LinBox
 *
 * Written by David Saunders 
 *
 * --------------------------------------------------------
 * See COPYING for license information
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/field/PIR-modular-int32.h"
//#include "linbox/field/PIR-modular-double.h"
#include "linbox/field/local2_32.h"
#include "linbox/blackbox/dense.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/vector/stream.h"
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
class foobar 
{
	public:
	typedef typename LocalPIR::Element first_argument_type;
	typedef LocalPIR second_argument_type;
	typedef void result_type;
	void operator()(typename LocalPIR::Element& d, const LocalPIR& R) const
	{ 
		typename LocalPIR::Element x = d;
		if (R.isUnit(d)) R.divin(d, d);
		else R.gcd(d, x, x);
	}
};

template<>
class foobar<LinBox::Local2_32> 
{
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
class pplt
{   
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
class pplt<LinBox::NTL_PID_zz_p> 
{
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
static bool testLocalSmith (const LocalPIR &R, VectorStream<vector<typename LocalPIR::Element> > &stream, string s) 
{
	typedef vector <typename LocalPIR::Element> Vector;
	typedef typename LocalPIR::Element Elt;
	typedef DenseMatrix<LocalPIR> Blackbox;

	commentator.start ("Testing local smith on random dense matrices", "testLocalSmith", stream.m ());
	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << s << endl;

	VectorDomain<LocalPIR> VD (R);

	bool ret = true;
	size_t i;
	size_t n = stream.n();

	Vector d;

	VectorWrapper::ensureDim (d, stream.dim ());

	while (stream) {
		commentator.startIteration (stream.j ());

		stream.next (d);

 		report << "Input vector:  ";
 		VD.write (report, d);
 		report << endl;

		Blackbox Lm (R, n, n), D (R, n, n), U (R, n, n), A (R, n, n);
		for( i = 0; i < n; ++i ) {D[i][i] = d[i];Lm[i][i]=U[i][i]=1;}
		
		size_t j;
		
		for (i = 0; i < n; ++ i) 
		       for ( j = 0; j < i; ++ j) {
			       
			       D[i][j] = D[j][i] = 0;
			       
			       Lm[i][j] = rand() % 10;
			       Lm[j][i] = 0;
			       
			       U[j][i] = rand() % 10;
			       U[i][j] = 0;
		       }

		MatrixDomain<LocalPIR> MR(R);
		
		Timer timer;
		
		timer.start();
		MR.mul(A,Lm,D);

		MR.mulin(A,U);
		timer.stop();
		report << "Two matrix multiplication: " << timer << "\n";
		
		//for( i = 0; i < n; ++i ) D[i][i] = rand() % 10 + 1;

		list< typename LocalPIR::Element > L;
		SmithFormLocal< LocalPIR > SmithForm;
		timer.start();
		SmithForm( L, A, R );
		timer.stop();
		report << "Time " << timer <<"\n";
			
		report.flush();
		report << "Computed invariants: ";
		
		report << "[";
		typedef typename list<Elt>::iterator listptr;
		for (listptr p = L.begin(); p != L.end(); ++p)
		    report << *p << ", ";
		report << "\b\b]" << endl;

		pplt<LocalPIR> lt(R);
		report << "normalize done" << endl;
		report.flush();

		for_each(d.begin(), d.end(), bind2nd(foobar<LocalPIR>(), R));
		timer.start();
		stable_sort(d.begin(), d.end(), lt);
		timer.stop();
		report << "Sorting " << timer <<"\n";

		report << "sort done" << endl;
		report.flush();

		report << "True invariants: ";
		VD.write (report, d);
		report << endl;
		report << flush;

		if ( L.size() != D.rowdim() ) {ret = false; break;}
		typedef typename Vector::iterator vectptr;
		listptr p; vectptr q;
		for (p = L.begin(), q = d.begin(); 
		     q != d.end(); 
		     ++p, ++q)
		    if ( !R.areEqual (*p, *q ) )
		    {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed invariants incorrect" << endl;
			 ret = false;
		    }
		commentator.stop("done");
		commentator.progress();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalTrace");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 2;
	static integer q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	typedef PIRModular<int32> Ring;
	//typedef PIRModular<dense> Ring;
	typedef vector<Ring::Element> Vector;

	Ring R (32768);
// 	Ring R (536870912);

	commentator.start("Local Smith Form test suite", "LocalSmith");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	RandomDenseStream<Ring, Vector> stream (R, n, iterations);

	if (!testLocalSmith<Ring> (R, stream, "PIRModular<int32>")) pass = false;

	// power of 2 test
	Local2_32 R2;
	RandomDenseStream<Local2_32, vector<Local2_32::Element> > 
		stream2 (R2, n, iterations);
	if (!testLocalSmith<Local2_32> (R2, stream2, "Local2_32")) pass = false;

	commentator.stop("Local Smith Form test suite");
	return pass ? 0 : -1;
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
