/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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

#ifndef __LINBOX_DIOPHANTINE_SOLVER_INL
#define __LINBOX_DIOPHANTINE_SOLVER_INL

#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h>

#include <linbox-config.h>

//#define DEBUG_DIO
//#define INFO_DIO

#define MONTE_CARLO_BOREDOM 21

namespace LinBox {


	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus DiophantineSolver<QSolver>::solve 
	(Vector1& x, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) {
		
		SolverReturnStatus result = _rationalSolver.solve(x, A, b, false, maxPrimes, level);
		if (result == SS_INCONSISTENT && level >= SL_CERTIFIED) 
			lastCertificate.copy(_rationalSolver.lastCertificate);
		return result;
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::randomSolve 
	(Vector1& x, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) {
		
		SolverReturnStatus result = _rationalSolver.findRandomSolution(x, A, b, maxPrimes, level);
		if (result == SS_INCONSISTENT && level >= SL_CERTIFIED) 
			lastCertificate.copy(_rationalSolver.lastCertificate);
		return result;
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::diophantineSolve 
	(Vector1& x, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) {

		//here maxPrimes is only used to bound trials of initial solution
		SolverReturnStatus status;
		
		//this should eliminate all inconsistent systems; when level == SL_MONTECARLO maybe not.
		status = _rationalSolver.monolithicSolve(x, A, b, (level >= SL_LASVEGAS), true, maxPrimes, level);
		if (status != SS_OK) {
			if (status == SS_FAILED) 
				cout << "WARNING, failed to find original solution; is maxPrimes > 1?" << endl;
			if (status == SS_INCONSISTENT && level >= SL_CERTIFIED) 
				lastCertificate.copy(_rationalSolver.lastCertificate);
			return status;
		}

		VectorFraction<Ring> y(_R, x), y0(y);

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
		cout << "lower bound on denominator: " << lowerDenBound << endl;
		cout << "upper bound on denominator: " << upperDenBound << endl;
#endif
		int numSolutionsNeeded = 1;
		int numFailedCallsToSolver = 0;
		int boredom = 0; //used in monte carlo, when we assume there's a diophantine solution
		while (! _R.areEqual(upperDenBound, lowerDenBound)) {
			_rationalSolver.chooseNewPrime();
			status = _rationalSolver.monolithicSolve(x, A, b, (level >= SL_LASVEGAS), true, 1, level);
			numSolutionsNeeded++;
			if (status != SS_OK) {
				numFailedCallsToSolver++;
				continue;
			}
			VectorFraction<Ring> yhat(_R, x);
#ifdef DEBUG_DIO	       
			cout << "random solution has denom: " << yhat.denom << endl;
#endif
			// goodCombination first represents whether a decrease in upperDenBound is achieved
			bool goodCombination = y.boundedCombineSolution(yhat, ODB, upperDenBound); 
#ifdef DEBUG_DIO
			if (goodCombination) 
				cout << "new gcd(denom, y0.denom): " << upperDenBound << endl;
#endif
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
// 					cout << "ERROR: got a 0 den factor" << endl;
// 					return SS_FAILED;
// 				}

				goodCombination = lastCertificate.combineCertificate
					(_rationalSolver.lastCertificate, n1, lowerDenBound,
					 _rationalSolver.lastZBNumer, 
					 _rationalSolver.lastCertifiedDenFactor);
			}
#ifdef DEBUG_DIO
			cout << "jonx";
			if (goodCombination) cout << "new certified denom factor: " << lowerDenBound << endl;
#endif
		}
#ifdef INFO_DIO
		cout << "number of solutions needed in total: " << numSolutionsNeeded << endl;
		cout << "number of failed calls to solver: " << numFailedCallsToSolver << endl;
#endif		
		y.combineSolution(y0);
		y.toFVector(x);
		return SS_OK;
	}
	
} //end of namespace LinBox
#endif
