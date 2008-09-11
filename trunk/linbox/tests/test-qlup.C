
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-rank.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * -----------------------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General
 * Public License. See COPYING for more information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/vector/sparse.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

// tests 1 and 2 were certain diagonals - now deemed unnecessary.  -bds 2005Mar15
/* Test 3: Rank of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its rank using Gaussian
 * elimination (direct and blas) and Wiedemann's algorithm. Checks that the results match.
 */

template <class Field>
bool testQLUP(const Field &F, size_t n, unsigned int iterations, double sparsity = 0.05) 
{
	bool res = true;
	typedef SparseMatrix<Field, Sparse_Vector<typename Field::Element> > Blackbox;

	commentator.start ("Testing Sparse elimination qlup", "testQLUP", iterations);

	unsigned long Ni = n;
	unsigned long Nj = n;
	typename Field::RandIter generator (F);

	for (size_t i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		RandomSparseStream<Field, Sparse_Vector<typename Field::Element> > stream (F, generator, sparsity, n, n);
		

		Blackbox A (F, stream);

		std::ostream & report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);	

		F.write( report ) << endl; 
		A.write( report,FORMAT_MAPLE ) << endl;

		std::vector<typename Field::Element> u(Nj), v(Ni), w1(Nj), w2(Ni), w3(Ni), w(Ni);
		for(typename std::vector<typename Field::Element>::iterator it=u.begin();it!=u.end();++it) 
			generator.random (*it);

		A.apply(v,u);


		unsigned long rank;
		
		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );
		typename Field::Element determinant;
		Blackbox L(F, A.rowdim(), A.coldim());
		Permutation<Field> Q(A.rowdim(),F);
		Permutation<Field> P(A.coldim(),F);
		
		GD.QLUPin(rank, determinant, 
			  Q, L, A, P, 
			  A.rowdim(), A.coldim() );

		Q.apply(w, L.apply(w3, A.apply(w2, P.apply(w1,u) ) ) );
        
		bool error = false;
		typename std::vector<typename Field::Element>::const_iterator itv=v.begin();
		typename std::vector<typename Field::Element>::const_iterator itw=w.begin();
		for( ; itw!=w.end();++itw,++itv) {
			if (! F.areEqual(*itw,*itv) ) {
				error = true;
			}
		}
        
		if (error) {
			res = false;

			report << "ERROR : matrix(" << u.size() << ",1,[";
			for(typename std::vector<typename Field::Element>::const_iterator itu=u.begin(); itu!=u.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(typename std::vector<typename Field::Element>::const_iterator itv=v.begin(); itv!=v.end();++itv)
				report << *itv << ' ';
			report << "]  !=  [";
			for(typename std::vector<typename Field::Element>::const_iterator itw=w.begin(); itw!=w.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;
			
			
			report << "w1: [";
			for(typename std::vector<typename Field::Element>::const_iterator itw=w1.begin(); itw!=w1.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;
			report << "w2: [";
			for(typename std::vector<typename Field::Element>::const_iterator itw=w2.begin(); itw!=w2.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;
			report << "w3: [";
			for(typename std::vector<typename Field::Element>::const_iterator itw=w3.begin(); itw!=w3.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;
		}
            		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (res), (const char *) 0, "testQLUP");

	return res;
}

/* Test 4: Rank of zero and identity matrices by Wiedemann variants
 *
 */



int main (int argc, char **argv)
{

//     commentator.setMaxDetailLevel( 100000 );
//     commentator.setMaxDepth( 100000 );
   
	bool pass = true;

	static size_t n = 80;
	static integer q = 65519U;
	//static integer q = 1000003U;
	static int iterations = 2;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
        { 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("QLUP  test suite", "qlup");
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over Modular<uint32>" << endl; 
	Modular<LinBox::uint32> F (q);
	if (!testQLUP (F, n, iterations, sparsity)) pass = false;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
	<< "over Modular<double>" << endl; 
	Modular<double> G (q);
	if (!testQLUP (G, n, iterations, sparsity)) pass = false;


// 	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
// 	<< "over PID_integer" << endl; 
//         PID_integer R;
// 	if (!testRankMethods (R, n, iterations, sparsity)) pass = false;

	commentator.stop("QLUP test suite");
	return pass ? 0 : -1;
}
