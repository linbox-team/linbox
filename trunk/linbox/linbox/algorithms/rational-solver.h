/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lifting-container.h
 * Copyright (C) 2004 Zhendong Wan, Pascal Giorgi
 *
 * Written by Zhendong Wan  <wan@mail.eecis.udel.edu> 
 *         and Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * Modified by David Pritchard  <daveagp@mit.edu>
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

#ifndef __LINBOX_RATIONAL_SOLVER_H
#define __LINBOX_RATIONAL_SOLVER_H


#include <linbox/solutions/methods.h>
#include <linbox/blackbox/archetype.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/blackbox/compose.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/util/timer.h>
#define RSTIMING
#define DEFAULT_PRIMESIZE 14

namespace LinBox {
	
#define SINGULARITY_THRESHOLD 5
#define BAD_PRECONTITIONER_THRESHOLD 5
#define DEFAULT_MAXPRIMES 5
#define SL_DEFAULT SL_LASVEGAS
	
	/** _Ring integer ring
	 *  _Field, finite field for lifting
	 */

 	template<class Ring, class Field,class RandomPrime, class MethodTraits = DixonTraits>		
 	class RationalSolver {};

	// used as return type for solving routines
	enum SolverReturnStatus {
		SS_OK, SS_FAILED, SS_SINGULAR, SS_INCONSISTENT, SS_BAD_PRECONDITIONER
	};
    
	// used to determine what level of solving should be done
	// Monte Carlo: Try to solve if possible, but result is not guaranteed. 
	//              In any case a 0 denominator should not be returned.
	// Las Vegas  : Result should be guaranteed correct.
	// Certified  : Additionally, provide certificates that the result returned is correct.
	//              - if the return value is SS_INCONSISTENT, this means
	//                   lastCertificate satisfies lC.A = 0, lC.b != 0
	//              - if diophantine solving was called and the return value is SS_OK, this means
	//                   lastCertificate satisfies den(lC.A) = 1, den(lC.b) = den(answer)
	enum SolverLevel {
		SL_MONTECARLO, SL_LASVEGAS, SL_CERTIFIED
	};    // note: code may assume that each level is 'stronger' than the previous one

	/* RationalSolver for linears systems over a Ring
	 * using p-adic lifting and Wiedemann algorithm.
	 */
	template<class Ring, class Field,class RandomPrime>		
	class RationalSolver<Ring, Field, RandomPrime, WiedemannTraits> {

	public: 
		typedef Ring                                 RingType;
		typedef typename Ring::Element                Integer;
		typedef typename Field::Element               Element;
		typedef typename RandomPrime::Prime_Type        Prime;
		typedef std::vector<Element>              FPolynomial;

	protected:
		
		RandomPrime      _genprime;
		mutable Prime       _prime;
		WiedemannTraits    _traits;
		Ring                    _R; 
    	
	public:

		/* Constructor
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE), const WiedemannTraits& traits=WiedemannTraits()) : 
			_R(r), _genprime(rp), _traits(traits){_prime=_genprime.randomPrime();}
    
		/* Constructor with a prime
		 * @param p   , a Prime
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Prime& p, const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE), const WiedemannTraits& traits=WiedemannTraits()) : 
			_R(r), _genprime(rp), _prime(p), _traits(traits){}
    
		/** Solve a linear system Ax=b over quotient field of a ring		 
		 * giving a random solution if the system is singular and consistent.
		 * giving the unique solution if the system is non-singular.
		 *
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A    , Matrix of linear system
		 * @param b    , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b,const bool, int maxPrimes = DEFAULT_MAXPRIMES) const;
    
		/** Solve a nonsingular linear system Ax=b over quotient field of a ring.
		 * giving the unique solution of the system.
		 *
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A   , Matrix of linear system
		 * @param b   , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveNonsingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES) const;         

		/** Solve a singular linear system Ax=b over quotient field of a ring.
		 * giving a random solution if the system is singular and consistent.
		 *
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A   , Matrix of linear system
		 * @param b   , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveSingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES) const;	


		template <class IMatrix, class FMatrix, class IVector>
		void sparseprecondition (const Field&,
					 const IMatrix* ,
					 Compose< LambdaSparseMatrix<Ring>,Compose<IMatrix, LambdaSparseMatrix<Ring> > > *&,
					 const FMatrix*,
					 Compose<LambdaSparseMatrix<Field>,Compose<FMatrix,LambdaSparseMatrix<Field> > > *&,
					 const IVector&,
					 IVector&,
					 LambdaSparseMatrix<Ring> *&,
					 LambdaSparseMatrix<Ring> *&,
					 LambdaSparseMatrix<Field> *&,
					 LambdaSparseMatrix<Field> *&) const;


 
		/*
		template <class IMatrix, class FMatrix, class IVector, class FVector>
		void precondition (const Field&,
				   const IMatrix&,
				   BlackboxArchetype<IVector>*&,
				   const FMatrix*,
				   BlackboxArchetype<FVector>*&,
				   const IVector&,				   
				   IVector&,
				   BlackboxArchetype<IVector>*&,
				   BlackboxArchetype<IVector>*&) const; 
	*/
			

	}; // end of specialization for the class RationalSover with Wiedemann traits

#ifdef RSTIMING
	class DixonTimer {
	public: 
		mutable Timer ttSetup, ttRecon, ttGetDigit, ttGetDigitConvert, ttRingApply, ttRingOther;
		void clear() const {
			ttSetup.clear();
			ttRecon.clear();
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
			ttRingOther.clear();
			ttRingApply.clear();
		}

		template<class RR, class LC>
		void update(RR& rr, LC& lc) const {
			ttSetup += lc.ttSetup;
			ttRecon += rr.ttRecon;
			ttGetDigit += lc.ttGetDigit;
			ttGetDigitConvert += lc.ttGetDigitConvert;
			ttRingOther += lc.ttRingOther;
			ttRingApply += lc.ttRingApply;
		}
	};
#endif

	/* RationalSolver for linears systems over a Ring
	 * using p-adic lifting and Dixon algorithm.
	 */
	template<class Ring, class Field,class RandomPrime>		
	class RationalSolver<Ring, Field, RandomPrime, DixonTraits> {
	
	public:          
		
		typedef Ring                                 RingType;
		typedef typename Ring::Element               Integer;
		typedef typename Field::Element              Element;
		typedef typename RandomPrime::Prime_Type     Prime;

		// polymorphic 'certificate' generated when level >= SL_CERTIFIED
		// certificate of inconsistency when any solving routine returns SS_INCONSISTENT
		// certificate of minimal denominator when findRandomSolutionAndCertificate is called & return is SS_OK
		mutable VectorFraction<Ring>                 lastCertificate;     

		//next 2 fields generated only by findRandomSolutionAndCertificate, when return is SS_OK
		mutable Integer                              lastZBNumer;               //filled in if level >= SL_CERTIFIED
		mutable Integer                              lastCertifiedDenFactor;    //filled in if level >= SL_LASVEGAS
		//note: lastCertificate * b = lastZBNumer / lastCertifiedDenFactor, in lowest form
		
	protected:
		
		RandomPrime                     _genprime;
		mutable Prime                   _prime;
		Ring                            _R; 
#ifdef RSTIMING
		mutable Timer       
   		        tSetup,           ttSetup,
			tLQUP,            ttLQUP,
			tFastInvert,      ttFastInvert,        //only done in deterministic or inconsistent
			tCheckConsistency,ttCheckConsistency,        //includes lifting the certificate
			tMakeConditioner, ttMakeConditioner,
			tInvertBP,        ttInvertBP,              //only done in random
			tCheckAnswer,     ttCheckAnswer,
			tCertSetup,       ttCertSetup,        //remaining 3 only done when makeMinDenomCert = true
			tCertMaking,      ttCertMaking,

			tNonsingularSetup,ttNonsingularSetup,
			tNonsingularInv,  ttNonsingularInv,

			totalTimer;
		
		mutable DixonTimer
   		        ttConsistencySolve, ttSystemSolve, ttCertSolve, ttNonsingularSolve;
#endif
		
	public:
		
		/** Constructor
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE)) : 
			lastCertificate(r, 0), _genprime(rp), _R(r) 
		{
			_prime=_genprime.randomPrime(); 
#ifdef RSTIMING
			clearTimers();
#endif
		}
    
		
		/** Constructor, trying the prime p first
		 * @param p   , a Prime
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Prime& p, const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE)) : 
			lastCertificate(r, 0), _genprime(rp), _prime(p), _R(r) 
		{
#ifdef RSTIMING
			clearTimers();
#endif
		}
    
		
		/** Solve a linear system Ax=b over quotient field of a ring
		 * 
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A        , Matrix of linear system
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED), and (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. 
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate is in lastCertificate if (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, const bool = false, 
					 const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;
		
		/** overload so that the bool 'oldMatrix' argument is not accidentally set to true */
		template <class IMatrix, class Vector1, class Vector2>	
		SolverReturnStatus solve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes, 
					 const SolverLevel level = SL_DEFAULT) const {
			return solve (num, den, A, b, false, maxPrimes, level);
		}

		/** Solve a nonsingular, square linear system Ax=b over quotient field of a ring
		 * 
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A        , Matrix of linear system (it must be square)
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 *
		 * @return status of solution. 
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found, guaranteed correct. 
		 *   SS_SINGULAR - system appreared singular mod all primes. 
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveNonsingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, bool = false, 
						    int maxPrimes = DEFAULT_MAXPRIMES) const;

		/** Solve a general rectangular linear system Ax=b over quotient field of a ring. 
		 *  If A is known to be square and nonsingular, calling solveNonsingular is more efficient.
		 * 
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A        , Matrix of linear system
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED), and (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. 
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate is in lastCertificate if (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveSingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, 
						 int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;

		/** Find a random solution of the general linear system Ax=b over quotient field of a ring.
		 * 
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param A        , Matrix of linear system
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED), and (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. 
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate is in lastCertificate if (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus findRandomSolution(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, 
						      int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;
		
		/** Big solving routine to perform random solving and certificate generation.
		 * Same arguments and return as findRandomSolution, except
		 *
		 * @param num  , Vector of numerators of the solution
		 * @param den  , The common denominator. 1/den * num is the rational solution of Ax = b.
		 * @param randomSolution,  parameter to determine whether to randomize or not (since solveSingular calls this function as well)
		 * @param makeMinDenomCert,  determines whether a partial certificate for the minimal denominator of a rational solution is made
		 *
		 * When (randomSolution == true && makeMinDenomCert == true), 
		 *   If (level == SL_MONTECARLO) this function has the same effect as calling findRandomSolution.
		 *   If (level >= SL_LASVEGAS && return == SS_OK), lastCertifiedDenFactor contains a certified factor of the min-solution's denominator.
		 *   If (level >= SL_CERTIFIED && return == SS_OK), lastZBNumer and lastCertificate are updated as well.
		 *
		 */
		template <class IMatrix, class Vector1, class Vector2>	
		SolverReturnStatus monolithicSolve (Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, 
						    bool makeMinDenomCert, bool randomSolution,
						    int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;

		Ring getRing() const {return _R;}

		void chooseNewPrime() const {_prime = _genprime.randomPrime();}

#ifdef RSTIMING
		void clearTimers() const
		{
			ttSetup.clear();
			ttLQUP.clear();
			ttFastInvert.clear();
			ttCheckConsistency.clear();
			ttMakeConditioner.clear();
			ttInvertBP.clear();
			ttCheckAnswer.clear();
			ttCertSetup.clear();
			ttCertMaking.clear();
			ttNonsingularSetup.clear();
			ttNonsingularInv.clear();

   		        ttConsistencySolve.clear();
			ttSystemSolve.clear();
			ttCertSolve.clear();
			ttNonsingularSolve.clear();
		}

	public:

		inline std::ostream& printTime(const Timer& timer, const char* title, std::ostream& os, const char* pref = "") const {
			if (&timer != &totalTimer)
				totalTimer += timer;
			if (timer.count() > 0) {
				os << pref << title;
				for (int i=strlen(title)+strlen(pref); i<28; i++) 
					os << ' ';
				return os << timer << endl;
			}
			else
				return os;
		}

		inline std::ostream& printDixonTime(const DixonTimer& timer, const char* title, std::ostream& os) const{
			if (timer.ttSetup.count() > 0) {
				printTime(timer.ttSetup, "Setup", os, title);
				printTime(timer.ttGetDigit, "Field Apply", os, title);
				printTime(timer.ttGetDigitConvert, "Ring-Field-Ring Convert", os, title);
				printTime(timer.ttRingApply, "Ring Apply", os, title);
				printTime(timer.ttRingOther, "Ring Other", os, title);
				printTime(timer.ttRecon, "Reconstruction", os, title);
			}
			return os;
		}

		std::ostream& reportTimes(std::ostream& os) const {
			totalTimer.clear();
			printTime(ttNonsingularSetup, "NonsingularSetup", os);
			printTime(ttNonsingularInv, "NonsingularInv", os);
			printDixonTime(ttNonsingularSolve, "NS ", os);
			printTime(ttSetup , "Setup", os);
			printTime(ttLQUP , "LQUP", os);
			printTime(ttFastInvert , "FastInvert", os);
			printTime(ttCheckConsistency , "CheckConsistency", os);
			printDixonTime(ttConsistencySolve, "INC ", os);
			printTime(ttMakeConditioner , "MakeConditioner", os);
			printTime(ttInvertBP , "InvertBP", os);
			printDixonTime(ttSystemSolve, "SYS ", os);
			printTime(ttCheckAnswer , "CheckAnswer", os);
			printTime(ttCertSetup , "CertSetup", os);
			printDixonTime(ttCertSolve, "CER ", os);
			printTime(ttCertMaking , "CertMaking", os);
			printTime(totalTimer , "TOTAL", os);
			return os;
		}
#endif
		
	}; // end of specialization for the class RationalSover with Dixon traits

}
#include <linbox/algorithms/rational-solver.inl>

#endif
