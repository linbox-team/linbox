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

#define DEBUG_DIO

namespace LinBox {


	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>
	SolverReturnStatus DiophantineSolver<QSolver>::solve (Vector1& x,
							       const IMatrix& A,
							       const Vector2& b, 
							       int maxPrimes = DEFAULT_MAXPRIMES) {
		
		return _rationalSolver.solve(x, A, b, false, maxPrimes);
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::randomSolve (Vector1& x,
								     const IMatrix& A,
								     const Vector2& b, 
								     int maxPrimes = DEFAULT_MAXPRIMES) {
		
		return _rationalSolver.findRandomSolution(x, A, b, maxPrimes);
	}
	
	template<class QSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<QSolver>::diophantineSolve (Vector1& x,
										const IMatrix& A,
										const Vector2& b,
										bool makeCertificate = false,
										int maxPrimes = DEFAULT_MAXPRIMES) {
		//here maxPrimes is only used to bound trials of initial solution
		SolverReturnStatus status, status2;
		
		status = _rationalSolver.findRandomSolutionAndCertificate(x, A, b, makeCertificate, true, maxPrimes);
		if (status != SS_OK) {
			if (status == SS_FAILED) {
				cout << "WARNING, failed to find original solution; is maxPrimes > 1?" << endl;
			}
			return status;
		}

		bool error;		
		VectorFraction<Ring> y(_R, x, error), y0(y);
		if (error) {cout << "dammit\n"; return SS_FAILED;}

		Integer ODB = y0.denom, n1; //ODB -- original denominator bound. equal to g(y0) from M+S. 
		if (makeCertificate) {
			lastCertificate.copy(_rationalSolver.lastCertificate);
			_R.assign(n1, _rationalSolver.lastZBNumer);
		}

		Integer upperDenBound = ODB;
		Integer lowerDenBound = _rationalSolver.lastCertifiedDenFactor;
#ifdef DEBUG_DIO	       
		cout << "lower bound on denominator: " << lowerDenBound << endl;
		cout << "upper bound on denominator: " << upperDenBound << endl;
#endif
		int numSolutionsNeeded = 1;
		while (! _R.areEqual(upperDenBound, lowerDenBound)) {
			_rationalSolver.chooseNewPrime();
			status2 = _rationalSolver.findRandomSolutionAndCertificate(x, A, b, makeCertificate, true, 1);
			numSolutionsNeeded++;
			if (status2 == SS_OK) {
				VectorFraction<Ring> yhat(_R, x, error);
				if (!error) {
#ifdef DEBUG_DIO	       
					// cout << "random solution has denom: " << yhat.denom << endl;
#endif
					// reduce denominator of solution
					bool goodCombination = y.boundedCombineSolution(yhat, ODB, upperDenBound); 
#ifdef DEBUG_DIO
					if (goodCombination) 
						cout << "new gcd(denom, y0.denom): " << upperDenBound << endl;
#endif
					// increase size of lower bound on denominator
					if (!makeCertificate) {
#ifdef DEBUG_DIO
						goodCombination =
							!_R.isDivisor(lowerDenBound, _rationalSolver.lastCertifiedDenFactor);
#endif
						_R.lcmin(lowerDenBound, _rationalSolver.lastCertifiedDenFactor);
					}
					else {
						if (_R.isZero(_rationalSolver.lastCertifiedDenFactor)) {
							cout << "ERROR: got a 0 den factor" << endl;
							return SS_FAILED;
						}
// 						cout << "new z: ";
// 						_rationalSolver.lastCertificate.write(cout);
// 						cout << endl << "zb = " << _rationalSolver.lastZBNumer <<'/' << _rationalSolver.lastCertifiedDenFactor << endl;
						goodCombination = lastCertificate.combineCertificate(_rationalSolver.lastCertificate, 
												     n1, lowerDenBound,
												     _rationalSolver.lastZBNumer, 
												     _rationalSolver.lastCertifiedDenFactor);
						}
// 					cout << "new certified denom factor: " << lowerDenBound << endl;
// 					{
// 						VectorFraction<Ring> tmp(y);
// 						tmp.combine(y0);
// 						cout << "solution right now: ";
// 						tmp.write(cout) << endl;
						
// 					}
#ifdef DEBUG_DIO
					if (goodCombination) cout << "new certified denom factor: " << lowerDenBound << endl;
#endif
				}
				else {
					cout << "dammit2\n";
				}
			}
		}
		y.combineSolution(y0);
		cout << "number of solutions needed in total: " << numSolutionsNeeded << endl;
		
		y.toFVector(x);
		return status;
	}
	
} //end of namespace LinBox
#endif
