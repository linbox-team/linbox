/* linbox/algorithms/diophantine-solver.inl
 * Copyright (C) 2004 David Pritchard
 *
 * Written by David Pritchard <daveagp@mit.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_diophantine_solver_INL
#define __LINBOX_diophantine_solver_INL

#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h>

#include <linbox/linbox-config.h>

//#define DEBUG_DIO
//#define INFO_DIO

#define MONTE_CARLO_BOREDOM 21

namespace LinBox 
{


	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus DiophantineSolver<QSolver>::solve 
	(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes, const SolverLevel level) {
		
		SolverReturnStatus result = _rationalSolver.solve(x, den, A, b, false, maxPrimes, level);
		if (result == SS_INCONSISTENT && level >= SL_CERTIFIED) 
			lastCertificate.copy(_rationalSolver.lastCertificate);
		return result;
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::randomSolve 
	(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes, const SolverLevel level) {
		
		SolverReturnStatus result = _rationalSolver.findRandomSolution(x, den, A, b, maxPrimes, level);
		if (result == SS_INCONSISTENT && level >= SL_CERTIFIED) 
			lastCertificate.copy(_rationalSolver.lastCertificate);
		return result;
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::diophantineSolve 
	(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes, const SolverLevel level) {

		//here maxPrimes is only used to bound trials of initial solution
		SolverReturnStatus status;
		
		//this should eliminate all inconsistent systems; when level == SL_MONTECARLO maybe not.
		status = _rationalSolver.monolithicSolve(x, den, A, b, (level >= SL_LASVEGAS), true, maxPrimes, level);
		if (status != SS_OK) {
			if (status == SS_FAILED && maxPrimes > 2) 
				std::cout << "ERROR, failed to find original solution and maxPrimes is not too small!" << std::endl;
			if (status == SS_INCONSISTENT && level >= SL_CERTIFIED) 
				lastCertificate.copy(_rationalSolver.lastCertificate);
			return status;
		}

		VectorFraction<Ring> y(_R,x.size());
		y. numer = x;
		y. denom = den;
		VectorFraction<Ring> y0(y);

		Integer ODB = y0.denom, n1; //ODB -- original denominator bound. equal to g(y0) from Muld+Storj. 
		if (level >= SL_CERTIFIED) {
			lastCertificate.copy(_rationalSolver.lastCertificate);
			_R.assign(n1, _rationalSolver.lastZBNumer);
		}

		Integer upperDenBound = ODB;
		Integer lowerDenBound;
		if (level >= SL_LASVEGAS) 
			lowerDenBound = _rationalSolver.lastCertifiedDenFactor;
		else
			_R.init(lowerDenBound, 1);
#ifdef DEBUG_DIO	       
		std::cout << "lower bound on denominator: " << lowerDenBound << std::endl;
		std::cout << "upper bound on denominator: " << upperDenBound << std::endl;
#endif
		numSolutionsNeeded     = 1;
		numFailedCallsToSolver = 0;
		numRevelantSolutions=1;
		int boredom = 0; //used in monte carlo, when we assume there's a diophantine solution
		while (! _R.areEqual(upperDenBound, lowerDenBound)) {
			_rationalSolver.chooseNewPrime();
			status = _rationalSolver.monolithicSolve(x, den, A, b, (level >= SL_LASVEGAS), true, 1, level);
			numSolutionsNeeded++;
#ifdef DEBUG_DIO	       
			std::cout << '.' ;
#endif
			if (status != SS_OK) {
				numFailedCallsToSolver++;
				continue;
			}
			VectorFraction<Ring> yhat(_R, x.size());
			yhat. numer = x;
			yhat. denom = den;
			// goodCombination first represents whether a decrease in upperDenBound is achieved
			bool goodCombination = y.boundedCombineSolution(yhat, ODB, upperDenBound); 

			if (goodCombination) {
				numRevelantSolutions++;
#ifdef DEBUG_DIO
				std::cout << "new gcd(denom, y0.denom): " << upperDenBound << std::endl;
#endif
			}
			// now, goodCombination will be updated as to whether there is an increase in lowerDenBound
			if (level == SL_MONTECARLO) { 
				if (goodCombination)
					boredom = 0;
				else 
					boredom++;
				if (boredom > MONTE_CARLO_BOREDOM) 
					break; //exit while loop
				goodCombination = false;          //since we dont update lowerDenBound, no increase happens
			}
			else if (level == SL_LASVEGAS) {
#ifdef DEBUG_DIO
				goodCombination =
					!_R.isDivisor(lowerDenBound, _rationalSolver.lastCertifiedDenFactor);
#endif
				_R.lcmin(lowerDenBound, _rationalSolver.lastCertifiedDenFactor);
			}
			else { //level == SL_CERTIFIED

// 				paranoid check
// 				if (_R.isZero(_rationalSolver.lastCertifiedDenFactor)) {
// 					std::cout << "ERROR: got a 0 den factor" << std::endl;
// 					return SS_FAILED;
// 				}

				goodCombination = lastCertificate.combineCertificate
					(_rationalSolver.lastCertificate, n1, lowerDenBound,
					 _rationalSolver.lastZBNumer, 
					 _rationalSolver.lastCertifiedDenFactor);
			}
#ifdef DEBUG_DIO
			if (goodCombination) 
				std::cout << "new certified denom factor: " << lowerDenBound << std::endl;
#endif
		}
#ifdef INFO_DIO
		std::cout << "number of solutions needed in total: " << numSolutionsNeeded << std::endl;
		std::cout << "number of failed calls to solver: " << numFailedCallsToSolver << std::endl;
#endif		
		y.combineSolution(y0);
		//y.toFVector(x);
		x   = y.numer;
		den = y.denom;
		return SS_OK;
	}
	
} //end of namespace LinBox

#undef MONTE_CARLO_BOREDOM

#endif //__LINBOX_diophantine_solver_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
