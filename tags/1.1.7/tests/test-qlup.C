/* tests/test-qlup.C
 * Copyright (C) The LinBox group
 *
 * Time-stamp: <22 Jun 10 15:59:56 Jean-Guillaume.Dumas@imag.fr>
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
#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/givaro-zpz.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: LQUP decomposition of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination. Checks that the results match.
 */
template <class Field, class Blackbox, class RandStream >
bool testQLUP(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05) 
{
	bool res = true;

	commentator.start ("Testing Sparse elimination qlup", "testQLUP", iterations);

	unsigned long Ni = n;
	unsigned long Nj = n;
        integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n);

	for (size_t i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		
		stream.reset();

		Blackbox A (F, stream);

		std::ostream & report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);	

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

/* Test 2: LQUP solve of a random sparse matrix and a random dense vector
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination. 
 * Then solve using the decomposition and checks that the results match.
 */
template <class Field, class Blackbox, class RandStream>
bool testQLUPsolve(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05) 
{
	bool res = true;

	commentator.start ("Testing Sparse elimination qlup solve", "testQLUPsolve", iterations);

	unsigned long Ni = n;
	unsigned long Nj = n;
        integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n);
		
        GF2 F2; GF2::RandIter bitgenerator(F2,2,rseed); GF2::Element randomsolve;

	for (size_t i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		stream.reset();
		Blackbox A (F, stream);

		std::ostream & report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);	

		F.write( report ) << endl; 
		A.write( report, FORMAT_MAPLE ) << endl;

		std::vector<typename Field::Element> u(Nj), v(Ni), x(Nj), y(Ni);
		for(typename std::vector<typename Field::Element>::iterator it=u.begin();it!=u.end();++it) 
			generator.random (*it);

		A.apply(v,u);


		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );
		
		Blackbox CopyA ( A );

		GD.solvein(x, A, v, bitgenerator.random(randomsolve) );
                report << "Random solving: " << randomsolve << std::endl;
		
		CopyA.apply(y, x);
        
		VectorDomain<Field> VD(F);
        
		
		if (! VD.areEqual(v,y)) {
			res=false;
                        A.write( report, FORMAT_MAPLE ) << endl;
			
			report << "ERROR v: matrix(" << v.size() << ",1,[";
			for(typename std::vector<typename Field::Element>::const_iterator itu=v.begin(); itu!=v.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			report << "ERROR y: matrix(" << y.size() << ",1,[";
			for(typename std::vector<typename Field::Element>::const_iterator itu=y.begin(); itu!=y.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(typename std::vector<typename Field::Element>::const_iterator itv=x.begin(); itv!=x.end();++itv)
				report << *itv << ' ';
			report << "]  !=  [";
			for(typename std::vector<typename Field::Element>::const_iterator itw=y.begin(); itw!=y.end();++itw)
				report << *itw << ' ';
			report << "]" << std::endl;
			
		}
            		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (res), (const char *) 0, "testQLUPsolve");

	return res;
}



/* Test 2: LQUP nullspacebasis of a random sparse matrix 
 *
 * Constructs a random sparse matrix and computes its QLUP decomposition
 * using Sparse Gaussian elimination (stores only U and P). 
 * Then solve using the decomposition and checks that the results match.
 */
template <class Field, class Blackbox, class RandStream>
bool testQLUPnullspace(const Field &F, size_t n, unsigned int iterations, int rseed, double sparsity = 0.05) 
{
	bool res = true;

	commentator.start ("Testing Sparse elimination qlup nullspacebasis", "testQLUPnullspace", iterations);

	unsigned long Ni = n;
	unsigned long Nj = n;
        integer card; F.cardinality(card);
	typename Field::RandIter generator (F,card,rseed);
	RandStream stream (F, generator, sparsity, n, n, rseed);
		
	for (size_t i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		stream.reset();
		Blackbox A (F, stream);

		std::ostream & report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);	

		F.write( report ) << endl; 
		A.write( report, FORMAT_MAPLE ) << endl;


		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		GaussDomain<Field> GD ( F );
		
		Blackbox CopyA ( A );
                Blackbox X(F, A.coldim(), A.coldim() );
                
		GD.nullspacebasisin(X, CopyA );

                unsigned long nullity = X.coldim();
                
                std::vector<typename Field::Element> u(nullity);
                for(typename std::vector<typename Field::Element>::iterator it=u.begin();it!=u.end();++it) 
			generator.random (*it);
                std::vector<typename Field::Element> v(Nj);
                X.apply(v,u);
                report << "Random combination of the rows of the NullSpace basis" << std::endl;
                
		std::vector<typename Field::Element> w(Ni);
                A.apply(w, v);
		
		VectorDomain<Field> VD(F);

		if (! VD.isZero(w)) {
			res=false;
                        A.write( report, FORMAT_MAPLE ) << endl;
			
			report << "ERROR u: matrix(" << u.size() << ",1,[";
			for(typename std::vector<typename Field::Element>::const_iterator itu=u.begin(); itu!=u.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			report << "ERROR v: matrix(" << v.size() << ",1,[";
			for(typename std::vector<typename Field::Element>::const_iterator itu=v.begin(); itu!=v.end();++itu)
				report << *itu << ',';
			report << "]);\n[";
			for(typename std::vector<typename Field::Element>::const_iterator itv=w.begin(); itv!=w.end();++itv)
				report << *itv << ' ';
			report << "]  !=  0" << std::endl;
			
		}
            		
		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (res), (const char *) 0, "testQLUPnullspace");

	return res;
}



int main (int argc, char **argv)
{

	commentator.setMaxDepth (-1);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
// 	commentator.setMaxDetailLevel( 100000 );
// 	commentator.setMaxDepth( 100000 );
   
	bool pass = true;

	static size_t n = 80;
	static integer q = 65519U;
        static integer bigQ("1234567890123456789012345678901234568123");
	//static integer q = 1000003U;
	static int iterations = 2;
        static double sparsity = 0.05;
        static int rseed = time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
                { 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
                { 'r', "-r R", "Random generator seed.", TYPE_INT,     &rseed },
				{ '\0' }
	};

	parseArguments (argc, argv, args);
	srand (rseed);

	commentator.start("QLUP  test suite", "qlup");
        commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
            << "Seed: " << rseed << endl; 
        
        { 
            commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
                << "over Modular<uint32>" << endl; 
            typedef Modular<LinBox::uint32> Field;
            Field F (q);
            typedef SparseMatrix<Field, Sparse_Vector<Field::Element> > Blackbox;
            typedef RandomSparseStream<Field, Sparse_Vector<Field::Element> > RandStream;
            if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
        }
        
        {
            commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
                << "over Modular<double>" << endl; 
            typedef Modular<double> Field;
            Field F (q);
            typedef SparseMatrix<Field, Sparse_Vector<Field::Element> > Blackbox;
            typedef RandomSparseStream<Field, Sparse_Vector<Field::Element> > RandStream;
            if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
        }
        
        {
            
            commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
                << "over GivaroZpz<Integer>" << endl; 
            typedef GivaroZpz<Integer> Field;
            Field F (bigQ);
            typedef SparseMatrix<Field, Sparse_Vector<Field::Element> > Blackbox;
            typedef RandomSparseStream<Field, Sparse_Vector<Field::Element> > RandStream;
            if (!testQLUP<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPsolve<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPnullspace<Field, Blackbox, RandStream> (F, n, iterations, rseed, sparsity)) pass = false;
        }
 
        {
            commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
                << "specialized over GF2>" << endl; 
            typedef GF2 Field;
            Field F2;
            typedef LinBox::GaussDomain<LinBox::GF2>::Matrix Blackbox;
            typedef RandomSparseStreamGF2<Blackbox::Row_t> RandStream;
            if (!testQLUP<Field, Blackbox, RandStream> (F2, n, iterations, rseed, sparsity)) pass = false;
            if (!testQLUPsolve<Field, Blackbox, RandStream> (F2, n, iterations, rseed, sparsity)) pass = false;
        }

	commentator.stop("QLUP test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
