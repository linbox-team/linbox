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
#define DEFAULT_PRIMESIZE 14

namespace LinBox {
	
#define SINGULARITY_THRESHOLD 5
#define BAD_PRECONTITIONER_THRESHOLD 5
#define DEFAULT_MAXPRIMES 5
	
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
#define  SL_DEFAULT SL_LASVEGAS

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
		 * @param A    , Matrix of linear system
		 * @param x    , Vector in which to store solution
		 * @param b    , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& answer, const IMatrix& A, const Vector2& b,const bool, int maxPrimes = DEFAULT_MAXPRIMES) const;
    
		/** Solve a nonsingular linear system Ax=b over quotient field of a ring.
		 * giving the unique solution of the system.
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveNonsingular(Vector1& answer, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES) const;         

		/** Solve a singular linear system Ax=b over quotient field of a ring.
		 * giving a random solution if the system is singular and consistent.
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 * @param maxPrimes , maximum number of moduli to try
		 *
		 * @return status of solution
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveSingular(Vector1& answer, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES) const;	


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
			

	}; // end of specialization for the class RationalSover with Wiedemann traits


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
		
	public:

		/** Constructor
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE)) : 
			lastCertificate(r, 0), _genprime(rp), _R(r) {_prime=_genprime.randomPrime(); }
    
		
		/** Constructor, trying the prime p first
		 * @param p   , a Prime
		 * @param r   , a Ring, set by default
		 * @param rp  , a RandomPrime generator, set by default		 
		 */
		RationalSolver (const Prime& p, const Ring& r = Ring(), const RandomPrime& rp = RandomPrime(DEFAULT_PRIMESIZE)) : 
			lastCertificate(r, 0), _genprime(rp), _prime(p), _R(r) {}
    
		
		/** Solve a linear system Ax=b over quotient field of a ring
		 * 
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
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
		SolverReturnStatus solve(Vector1& x, const IMatrix& A, const Vector2& b, const bool = false, 
					 const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;
		
		/** overload so that the bool 'oldMatrix' argument is not accidentally set to true */
		template <class IMatrix, class Vector1, class Vector2>	
		SolverReturnStatus solve(Vector1& answer, const IMatrix& A, const Vector2& b, const int maxPrimes, 
					 const SolverLevel level = SL_DEFAULT) const {
			return solve (answer, A, b, false, maxPrimes, level);
		}

		/** Solve a nonsingular, square linear system Ax=b over quotient field of a ring
		 * 
		 * @param A        , Matrix of linear system (it must be square)
		 * @param x        , Vector in which to store solution
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 *
		 * @return status of solution. 
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found, guaranteed correct. 
		 *   SS_SINGULAR - system appreared singular mod all primes. 
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solveNonsingular(Vector1& x, const IMatrix& A, const Vector2& b, bool, 
						    int maxPrimes = DEFAULT_MAXPRIMES) const;

		/** Solve a general rectangular linear system Ax=b over quotient field of a ring. 
		 *  If A is known to be square and nonsingular, calling solveNonsingular is more efficient.
		 * 
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
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
		SolverReturnStatus solveSingular(Vector1& x, const IMatrix& A, const Vector2& b, 
						 int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;

		/** Find a random solution of the general linear system Ax=b over quotient field of a ring.
		 * 
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
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
		SolverReturnStatus findRandomSolution(Vector1& x, const IMatrix& A, const Vector2& b, 
						      int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;
		
		/** Big solving routine to perform random solving and certificate generation.
		 * Same arguments and return as findRandomSolution, except
		 *
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
		SolverReturnStatus monolithicSolve (Vector1& answer, const IMatrix& A, const Vector2& b, 
						    bool makeMinDenomCert, bool randomSolution,
						    int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT) const;

		Ring getRing() const {return _R;}

		void chooseNewPrime() const {_prime = _genprime.randomPrime();}
		
	}; // end of specialization for the class RationalSover with Dixon traits

}
#include <linbox/algorithms/rational-solver.inl>

#endif
