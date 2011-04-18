/* linbox/solutions/solve.h
 * Copyright(C) LinBox
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_solve_H
#define __LINBOX_solve_H

#include <vector>
#include <algorithm>

// must fix this list...
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/algorithms/wiedemann.h"
#include "linbox/algorithms/rational-solver.h"
#include "linbox/algorithms/diophantine-solver.h"
#include "linbox/blackbox/dense.h"
#include "linbox/matrix/factorized-matrix.h"
#include "linbox/util/debug.h"
#include "linbox/util/error.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"
#include "linbox/algorithms/bbsolve.h"

#include "linbox/algorithms/rational-cra2.h"
#include "linbox/algorithms/varprec-cra-early-multip.h"

namespace LinBox 
{

	// for specialization with respect to the DomainCategory
	template< class Vector, class Blackbox, class SolveMethod, class DomainCategory>
	Vector & solve (Vector & 			x,
			const Blackbox &                A,
			const Vector &			b,
			const DomainCategory &        tag,
			const SolveMethod &            M);
	//		SolveStatus * 			s = 0);

	/** \brief Solve Ax = b, for x.
	 *
	 * Vector x such that Ax = b is returned.  
	 In the case of a singular matrix A, if the system is consistent, a random
	 solution is returned by default.  The method parameter may contain
	 an indication that an arbitrary element of the solution space is 
	 acceptable, which can be faster to compute.  
	 If the system is inconsistent the zero vector is returned. 
	 
         \ingroup solutions
        */
	//and the SolveStatus, if non-null, is set to indicate inconsistency.
	template< class Vector, class Blackbox, class SolveMethod>
	Vector & solve (Vector &        		x,
			const Blackbox &                A,
			const Vector &			b,
			const SolveMethod &             M)
	//		SolveStatus * 			s = 0)
	{ 
		return solve(x, A, b, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// the solve with default method
	template< class Vector, class Blackbox>
	Vector& solve(Vector& x, const Blackbox& A, const Vector& b)
	{ return solve(x, A, b, Method::Hybrid()); }

	// in methods.h FoobarMethod and Method::Foobar are the same class.
	// in methods.h template<BB> bool useBB(const BB& A) is defined.

	// specialize this on blackboxes which have local methods
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const Method::Hybrid& m)
	{	
		if (useBB(A)) return solve(x, A, b, Method::Blackbox(m)); 
		else return solve(x, A, b, Method::Elimination(m));
	}

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const Method::Blackbox& m)
	{ 
		// what is chosen here should be best and/or most reliable currently available choice
		// 		integer c; A.field().cardinality(c);
		// 		if (c < 100) return solve(x, A, b, Method::BlockLanczos(m));
		return solve(x, A, b, Method::Wiedemann(m));
	}

	// temporary - fix this
#define inBlasRange(p) true

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const Method::Elimination& m)
	{ 
		integer c, p;
		A.field().cardinality(c);
		A.field().characteristic(p);
		//if ( p == 0 || (c == p && inBlasRange(p)) )
		return solve(x, A, b, 
			     typename FieldTraits<typename BB::Field>::categoryTag(), 
			     Method::BlasElimination(m)); 
  		//else 
		//	return solve(x, A, b, 
		//			typename FieldTraits<typename BB::Field>::categoryTag(), 
		//			Method::NonBlasElimination(m)); 
	}

	template <class Vector, class Field> 
	Vector& solvein(Vector& x, SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq>& A, const Vector& b, const Method::SparseElimination& m)
	{
            commentator.start ("Sparse Elimination Solve In Place", "sesolvein");
            GaussDomain<Field> GD ( A.field() );
            GD.solvein(x, A, b);
            commentator.stop ("done", NULL, "sesolvein");
            return x;
	}


	// Change of representation to be able to call the sparse elimination
	template <class Vector, class Blackbox> 
	Vector& solve(Vector& x, const Blackbox& A, const Vector& b, 
                      const Method::SparseElimination& m)
	{
            typedef typename Blackbox::Field Field;
            typedef SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq> SparseBB;
            SparseBB SpA(A.field(), A.rowdim(), A.coldim());
            MatrixHom::map(SpA, A, A.field());
            return solvein(x, SpA, b, m);
	}

	template <class Vector> 
	Vector& solvein(Vector& x, 
                        GaussDomain<GF2>::Matrix    &A,
                        const Vector& b, 
                        const Method::SparseElimination& m)
	{
            commentator.start ("Sparse Elimination Solve In Place over GF2", "GF2sesolvein");
            GaussDomain<GF2> GD ( A.field() );
            GD.solvein(x, A, b);
            commentator.stop ("done", NULL, "GF2sesolvein");
            return x;
	}
	template <class Vector> 
	Vector& solve(Vector& x,
                      GaussDomain<GF2>::Matrix    &A,
                      const Vector& b, 
                      const Method::SparseElimination& m)
	{
                // We make a copy
            GaussDomain<GF2>::Matrix SpA(A.field(), A.rowdim(), A.coldim());
            MatrixHom::map(SpA, A, A.field());
            return solvein(x, SpA, b, m);
	}

	template <class Vector, class Field> 
	Vector& solve(Vector& x, const SparseMatrix<Field>& A, const Vector& b, 
		      const Method::Elimination& m)
	{	
		//             bool consistent = false;
		// sparse elimination based solver can be called here ?
        	// For now we call the dense one
        	
		return solve(x, A, b, 
			     typename FieldTraits<typename SparseMatrix<Field>::Field>::categoryTag(), 
			     Method::BlasElimination(m)); 

#if 0
		 		if ( ! consistent ) {  // we will return the zero vector
		 			typename Field::Element zero; A.field().init(zero, 0);
		 			for (typename Vector::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		 		}
		 		return x;
#endif
	}
	// BlasElimination section ///////////////////

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		BlasBlackbox<typename BB::Field> B(A); // copy A into a BlasBlackbox
		return solve(x, B, b, tag, m);
	} 

	template <class Vector, class Field> 
	Vector& solve(Vector& x, const BlasBlackbox<Field>& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Solving linear system (FFLAS LQUP)", "LQUP::left_solve");
		//		bool consistent = false;
		LQUPMatrix<Field> LQUP(A);
		//FactorizedMatrix<Field> LQUP(A);

		LQUP.left_solve(x, b);

		// this should be implemented directly in left_solve 
		//if ( ! consistent ) {  // we will return the zero vector
		//	typename Field::Element zero; A.field().init(zero, 0);
		//	for (typename Vector::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		//}
		commentator.stop ("done", NULL, "LQUP::left_solve");
		return x;
	} 

	
	/* Integer tag Specialization for Dixon method:
	 * 2 interfaces: 
	 *   - the output is a common denominator and a vector of numerator (no need of rational number)
	 *   - the output is a vector of rational
	 */  


	// error handler for bad use of the integer solver API
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::IntegerTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		std::cout<<"try to solve system over the integer\n"
			 <<"the API need either \n"
			 <<" - a vector of rational as the solution \n"
			 <<" - or an integer for the common denominator and a vector of integer for the numerators\n\n";
		throw LinboxError("bad use of integer API solver\n");
		
	} 
#if 0
	template <class RatVector, class Vector, class BB, class MethodTraits> 
	Vector& solve(RatVector& x, const BB& A, const Vector& b, 
		      const RingCategories::RationalTag & tag, 
		      const MethodTraits& m)
	{
                if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
	                throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Rational CRA Solve", "Rsolve");
		size_t bits = 26 -(int)ceil(log((double)A.rowdim())*0.7213475205);
		RandomPrimeIterator genprime( bits);

		RationalRemainder2< VarPrecEarlyMultipCRA< Modular<double> > > rra(3UL);//using default RR method
		IntegerModularSolve<BB,Vector,MethodTraits > iteration(A, b, m);
                integer den;
		std::vector< integer > num(A.coldim());
		
		rra(num, den, iteration, genprime);

		typename RatVector::iterator it_x= x.begin();
		typename std::vector<integer>::const_iterator it_num= num.begin();

		for (; it_x != x.end(); ++it_x, ++it_num){
			integer g = gcd( *it_num, den);
			*it_x = typename RatVector::value_type(*it_num/g, den/g);
		}

		commentator.stop ("done", NULL, "Rsolve");
		return x;
	}
#endif
	// error handler for non defined solver over rational domain
#if 0
	template <class Vector, class BB, class MethodTraits> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::RationalTag & tag, 
		      const MethodTraits& m)
	{ 
		throw LinboxError("LinBox ERROR: solver not yet defined over rational domain");
	}
#endif
	
	/*
	 * 1st integer solver API :
	 * solution is a vector of rational numbers
	 * RatVector is assumed to be the type of a vector of rational number
	 */	

	// default API (method is BlasElimination)
	template<class RatVector, class Vector, class BB>	
	RatVector& solve(RatVector& x, const BB &A, const Vector &b){
		return solve(x, A, b, Method::BlasElimination());
	}

	// API with Hybrid method
	template<class RatVector, class Vector, class BB>	
	RatVector& solve(RatVector& x, const BB &A, const Vector &b, const Method::Hybrid &m){
		if (useBB(A)) 
			return solve(x, A, b, Method::Blackbox(m)); 
		else 
			return solve(x, A, b, Method::Elimination(m));
	}

	// API with Blackbox method
	template<class RatVector, class Vector, class BB>	
	RatVector& solve(RatVector& x, const BB &A, const Vector &b, const Method::Blackbox &m){
		return solve(x, A, b, Method::Wiedemann(m));
	}

	// API with Elimination method
	template<class RatVector, class Vector, class BB>	
	RatVector& solve(RatVector& x, const BB &A, const Vector &b, const Method::Elimination &m){
		return solve(x, A, b,  Method::BlasElimination(m));
	}


	// launcher of specialized solver depending on the MethodTrait
	template<class RatVector, class Vector, class BB, class MethodTraits>	
	RatVector& solve(RatVector& x, const BB &A, const Vector &b, const MethodTraits &m){
		return solve(x, A, b, typename FieldTraits<typename BB::Field>::categoryTag(),  m);
	}


	/* Specializations for BlasElimination over the integers
	 */

	// input matrix is generic (copying it into a BlasBlackbox)
	template <class RatVector, class Vector, class BB> 
	RatVector& solve(RatVector& x, const BB& A, const Vector& b, 
			 const RingCategories::IntegerTag & tag, 
			 const Method::BlasElimination& m)
	{ 
		BlasBlackbox<typename BB::Field> B(A); // copy A into a BlasBlackbox
		return solve(x, B, b, tag, m);
	} 
	
	// input matrix is a BlasBlackbox (no copy)
	template <class RatVector, class Vector, class Ring> 
	RatVector& solve(RatVector& x, const BlasBlackbox<Ring>& A, const Vector& b, 
			 const RingCategories::IntegerTag & tag, 
			 const Method::BlasElimination& m)
	{ 
	
		Method::Dixon mDixon(m);
		typename Ring::Element d;
		std::vector< typename Ring::Element> num(A.coldim());
		solve (num, d, A, b, tag, mDixon);
		
		typename RatVector::iterator it_x= x.begin();
		typename std::vector< typename Ring::Element>::const_iterator it_num= num.begin();
		integer n,den;
		A.field().convert(den,d);
		for (; it_x != x.end(); ++it_x, ++it_num){			
			A.field().convert(n, *it_num);
			*it_x = typename RatVector::value_type(n, den);
		}
			
		return x;
	} 

	// input matrix is a DenseMatrix (no copy)
	template <class RatVector, class Vector, class Ring> 
	RatVector& solve(RatVector& x, const DenseMatrix<Ring>& A, const Vector& b, 
			 const RingCategories::IntegerTag & tag, 
			 const Method::BlasElimination& m)
	{ 
		Method::Dixon mDixon(m);
		typename Ring::Element d;
		std::vector< typename Ring::Element> num(A.coldim());
		solve (num, d, A, b, tag, mDixon);
		typename RatVector::iterator it_x= x.begin();
		typename std::vector< typename Ring::Element>::const_iterator it_num= num.begin();
		integer n,den;
		A.field().convert(den,d); 
		for (; it_x != x.end(); ++it_x, ++it_num){			
			A.field().convert(n, *it_num);
			*it_x = typename RatVector::value_type(n, den);
		}
		
		return x;
	} 

	/*
	 * 2nd integer solver API :
	 * solution is a formed by a common denominator and a vector of integer numerator
	 * solution is num/d
	 */
	
	
	// default API (method is BlasElimination)
	template< class Vector, class BB>
	Vector& solve(Vector &x, typename BB::Field::Element &d, const BB &A, const Vector &b){
		return solve(x, d, A, b, typename FieldTraits<typename BB::Field>::categoryTag(),  Method::BlasElimination());
	}
		
	// launcher of specialized solver depending on the MethodTraits
	template< class Vector, class BB, class MethodTraits>
	Vector& solve(Vector &x, typename BB::Field::Element &d, const BB &A, const Vector &b, const MethodTraits &m){
		return solve(x, d, A, b, typename FieldTraits<typename BB::Field>::categoryTag(), m);
	}
	
	/* Specialization for BlasElimination over the integers
	 */
	
	// input matrix is generic (copying it into a BlasBlackbox)
	template <class Vector, class BB> 
	Vector& solve(Vector& x, typename BB::Field::Element &d, const BB& A, const Vector& b, 
		      const RingCategories::IntegerTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		BlasBlackbox<typename BB::Field> B(A); // copy A into a BlasBlackbox
		return solve(x, d, B, b, tag, m);
	} 
	
	// input matrix is a BlasBlackbox (no copy)
	template <class Vector, class Ring> 
	Vector& solve(Vector& x, typename Ring::Element &d, 
                      const BlasBlackbox<Ring>& A, const Vector& b, 
		      const RingCategories::IntegerTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		Method::Dixon mDixon(m);
		return solve(x, d, A, b, tag, mDixon);
	} 

	// input matrix is a DenseMatrix (no copy)
	template <class Vector, class Ring> 
	Vector& solve(Vector& x, typename Ring::Element &d, 
                      const DenseMatrix<Ring>& A, const Vector& b, 
		      const RingCategories::IntegerTag & tag, 
		      const Method::BlasElimination& m)
	{ 
		Method::Dixon mDixon(m);
		return solve(x, d, A, b, tag, mDixon);
	}

	// input matrix is a SparseMatrix (no copy)
	template <class Vect, class Ring> 
	Vect& solve(Vect& x, typename Ring::Element &d, 
                    const SparseMatrix<Ring, typename Vector<Ring>::SparseSeq>& A, 
                    const Vect& b, 
                    const RingCategories::IntegerTag & tag, 
                    const Method::SparseElimination& m)
	{ 
            Method::Dixon mDixon(m);
            return solve(x, d, A, b, tag, mDixon);
	}



	/** \brief solver specialization with the 2nd API and DixonTraits over integer (no copying)
	 */
	template <class Vector, class Ring> 
	Vector& solve(Vector& x, typename Ring::Element &d, const BlasBlackbox<Ring>& A, const Vector& b, 
		      const RingCategories::IntegerTag tag, Method::Dixon& m)
	{ 
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Padic Integer Blas-based Solving", "solving");
		
		typedef Modular<double> Field;
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		RationalSolver<Ring, Field, RandomPrimeIterator, DixonTraits> rsolve(A.field(), genprime); 			
		SolverReturnStatus status = SS_OK;

		// if singularity unknown and matrix is square, we try nonsingular solver
		switch ( m.singular() ) {
		case Specifier::SINGULARITY_UNKNOWN:
			switch (A.rowdim() == A.coldim() ? 
				status=rsolve.solveNonsingular(x, d, A, b, false ,m.maxTries()) : SS_SINGULAR) {				
			case SS_OK:
				m.singular(Specifier::NONSINGULAR);				
				break;					
			case SS_SINGULAR:
				switch (m.solution()){
				case DixonTraits::DETERMINIST:
					status= rsolve.monolithicSolve(x, d, A, b, false, false, m.maxTries(), 
								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					break;					
				case DixonTraits::RANDOM:
					status= rsolve.monolithicSolve(x, d, A, b, false, true, m.maxTries(), 
								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					break;					
				case DixonTraits::DIOPHANTINE:
					{ 
						DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
						status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
										(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
                                        }
                                        break;					
				default:
					break;
				}			
				break;
			default:
                                break;
			}
			break;
			
		case Specifier::NONSINGULAR:
			rsolve.solveNonsingular(x, d, A, b, false ,m.maxTries());
			break;
			    
		case Specifier::SINGULAR:
			switch (m.solution()){
			case DixonTraits::DETERMINIST:
				status= rsolve.monolithicSolve(x, d, A, b, 
							       false, false, m.maxTries(), 
							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				break;
				
			case DixonTraits::RANDOM:
				status= rsolve.monolithicSolve(x, d, A, b, 
							       false, true, m.maxTries(), 
							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				break;
				
			case DixonTraits::DIOPHANTINE:
				{
					DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
					status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
									(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
                                }
				break;
				
			//default:
			//	break;
			}		
		default:			    
			break;
		}

		commentator.stop("done", NULL, "solving");

		if ( status == SS_INCONSISTENT ) {  
			throw LinboxMathInconsistentSystem("Linear system is inconsistent");
//			typename Ring::Element zero; A.field().init(zero, 0);
// 			for (typename Vector::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		}
		return x;
	}	

	/** \brief solver specialization with the 2nd API and DixonTraits over integer (no copying)
	 */
	template <class Vector, class Ring> 
	Vector& solve(Vector& x, typename Ring::Element &d, const DenseMatrix<Ring>& A, const Vector& b, 
		      const RingCategories::IntegerTag tag, Method::Dixon& m)
	{  
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");
		
		commentator.start ("Padic Integer Blas-based Solving", "solving");
		
		typedef Modular<double> Field;
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		RationalSolver<Ring, Field, RandomPrimeIterator, DixonTraits> rsolve(A.field(), genprime); 			
		SolverReturnStatus status = SS_OK;
		// if singularity unknown and matrix is square, we try nonsingular solver
		switch ( m.singular() ) {
		case Specifier::SINGULARITY_UNKNOWN:
			switch (A.rowdim() == A.coldim() ? 
				status=rsolve.solveNonsingular(x, d, A, b, false ,m.maxTries()) : SS_SINGULAR) {				
			case SS_OK:
				m.singular(Specifier::NONSINGULAR);				
				break;					
			case SS_SINGULAR:
				switch (m.solution()){
				case DixonTraits::DETERMINIST:
					status= rsolve.monolithicSolve(x, d, A, b, false, false, m.maxTries(), 
								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					break;					
				case DixonTraits::RANDOM:
					status= rsolve.monolithicSolve(x, d, A, b, false, true, m.maxTries(), 
								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					break;					
				case DixonTraits::DIOPHANTINE:
					DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
					status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
									(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					break;					
				//default:
				//	break;
				}			
				break;
			}
			
		case Specifier::NONSINGULAR:
			rsolve.solveNonsingular(x, d, A, b, false ,m.maxTries());
			break;
			    
		case Specifier::SINGULAR:
			switch (m.solution()){
			case DixonTraits::DETERMINIST:
				status= rsolve.monolithicSolve(x, d, A, b, 
							       false, false, m.maxTries(), 
							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				break;
				
			case DixonTraits::RANDOM:
				status= rsolve.monolithicSolve(x, d, A, b, 
							       false, true, m.maxTries(), 
							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				break;
				
			case DixonTraits::DIOPHANTINE:
				DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
				status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
								(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				break;
				
			//default:
			//	break;
			}		
		default:			    
			break;
		}

		commentator.stop("done", NULL, "solving");
		if ( status == SS_INCONSISTENT ) {  
			throw LinboxMathInconsistentSystem("Linear system is inconsistent");
// 			typename Ring::Element zero; A.field().init(zero, 0);
// 			for (typename Vector::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		}
		return x;	
	}	



	/** \brief solver specialization with the 2nd API and DixonTraits over integer (no copying)
	 */
	template <class Vect, class Ring> 
	Vect& solve(Vect& x, typename Ring::Element &d, 
                    const SparseMatrix<Ring, typename Vector<Ring>::SparseSeq> & A,
                    const Vect& b, 
                    const RingCategories::IntegerTag tag, 
                    Method::Dixon& m)
	{ 
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Padic Integer Sparse Elimination Solving", "solving");
		
		typedef Modular<double> Field;
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		RationalSolver<Ring, Field, RandomPrimeIterator, SparseEliminationTraits> rsolve(A.field(), genprime); 			
		SolverReturnStatus status = SS_OK;

		// if singularity unknown and matrix is square, we try nonsingular solver
		switch ( m.singular() ) {
		case Specifier::SINGULARITY_UNKNOWN:
			switch (A.rowdim() == A.coldim() ? 
				status=rsolve.solveNonsingular(x, d, A, b,m.maxTries()) : SS_SINGULAR) {				
			case SS_OK:
				m.singular(Specifier::NONSINGULAR);				
				break;					
#if 0
 			case SS_SINGULAR:
 				switch (m.solution()){
 				case DixonTraits::DETERMINIST:
 					status= rsolve.monolithicSolve(x, d, A, b, false, false, m.maxTries(), 
 								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
 					break;					
 				case DixonTraits::RANDOM:
 					status= rsolve.monolithicSolve(x, d, A, b, false, true, m.maxTries(), 
 								       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
 					break;					
 				case DixonTraits::DIOPHANTINE:
 					{ 
 						DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
 						status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
 										(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
					 }
					 break;					
 				default:
 					break;
 				}			
 				break;
#endif
			default:
                                break;
			}
			break;
			
		case Specifier::NONSINGULAR:
			rsolve.solveNonsingular(x, d, A, b, m.maxTries());
			break;
			    
		case Specifier::SINGULAR:
#if 0
 			switch (m.solution()){
 			case DixonTraits::DETERMINIST:
 				status= rsolve.monolithicSolve(x, d, A, b, 
 							       false, false, m.maxTries(), 
 							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
 				break;
				
 			case DixonTraits::RANDOM:
 				status= rsolve.monolithicSolve(x, d, A, b, 
 							       false, true, m.maxTries(), 
 							       (m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
 				break;
				
 			case DixonTraits::DIOPHANTINE:
 				{
 					DiophantineSolver<RationalSolver<Ring,Field,RandomPrimeIterator, DixonTraits> > dsolve(rsolve);
 					status= dsolve.diophantineSolve(x, d, A, b, m.maxTries(),
 									(m.certificate()? SL_LASVEGAS: SL_MONTECARLO));
				 }
 				break;
				
 			//default:
 			//	break;
 			}		
#endif
		default:			    
			break;
		}

		commentator.stop("done", NULL, "solving");

		if ( status == SS_INCONSISTENT ) {  
			throw LinboxMathInconsistentSystem("Linear system is inconsistent");
//			typename Ring::Element zero; A.field().init(zero, 0);
// 			for (typename Vect::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		}
		return x;
	}	



#if 0
	  struct BlasEliminationCRASpecifier;
	  // Extra case put in (1) for timing comparison or (2) for parallelism or 
	  // (3) as an example of how we might leave an abandoned choice around in a 
	  // callable state for future reference 
	  template <class Vector, class Field> 
	  Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
			const RingCategories::IntegerTag & tag, 
			const BlasEliminationCRASpecifier & m)
	  { // (low priority) J-G puts in code using CRA object CRA and solve(x, A, b, ModularTag, Method::BlasElimination) 
		  typename Field::Element zero; A.field().init(zero, 0);
		  for (typename Vector::iterator i = x.begin(); i != x.end(); ++i) *i = zero;
		  return x;
	  } 
#endif

	// NonBlasElimination section ////////////////

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::NonBlasElimination& m)
	{	DenseMatrix<typename BB::Field> B(A); // copy
		return solve(x, B, b, tag, m);
	}

	// specialization when no need to copy
	template <class Vector, class Field> 
	Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::NonBlasElimination& m)
	{ //Do we have a non blas elimination?  There was not one in the original solve.h (now algorithms/bbsolve.h).
		return x;
	}

	// note: no need for NonBlasElimination when RingCategory is integer

	// Lanczos ////////////////
	// may throw SolverFailed or InconsistentSystem
#if 0
	
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::Lanczos& m)
	{
		solve(A, x, b, A.field(), m);
		return x;
	}



	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::BlockLanczos& m)
	{
		try { 
			solve(A, x, b, A.field(), m); 
		} catch (SolveFailed) {
			typename BB::Field::Element zero; A.field().init(zero, 0);
			for (typename Vector::iterator i = x.begin(); 
			     i != x.end(); ++i) 
				*i = zero;
		}
		return x;
	}
#endif

	// Wiedemann section ////////////////

	// may throw SolverFailed or InconsistentSystem
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
		      const RingCategories::ModularTag & tag, 
		      const Method::Wiedemann& m)
	{
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");
		
		// adapt to earlier signature of wiedemann solver
		solve(A, x, b, A.field(), m);
		return x;
	}


	/* remark 1.  I used copy constructors when switching method types.
	   But if the method types are (empty) child classes of a common  parent class containing
	   all the information, then casts can be used in place of copies.
	*/ 

} // LinBox



#include "linbox/field/modular.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/algorithms/rational-cra-early-multip.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/vector/vector-traits.h"

namespace LinBox 
{
   

	template <class Blackbox, class Vector, class MyMethod>
	struct IntegerModularSolve {       
		const Blackbox &A;
		const Vector &B;
		const MyMethod &M;

		IntegerModularSolve(const Blackbox& b, const Vector& v, const MyMethod& n) 
			: A(b), B(v), M(n) {}
        
        
		template<typename Field>
		typename Rebind<Vector, Field>::other& operator()(typename Rebind<Vector, Field>::other& x, const Field& F) const {
			typedef typename Blackbox::template rebind<Field>::other FBlackbox;
			FBlackbox Ap(F, A.rowdim(), A.coldim());
			MatrixHom::map(Ap, A, F);

			typedef typename Rebind<Vector, Field>::other FVector;
			Hom<typename Blackbox::Field, Field> hom(A.field(), F);
			FVector Bp(B.size());
			typename Vector::const_iterator Bit = B.begin();
			typename FVector::iterator      Bpit = Bp.begin();
			for( ; Bit != B.end(); ++Bit, ++Bpit)
				hom.image (*Bpit, *Bit);

			VectorWrapper::ensureDim (x, A.coldim());
			return solve( x, Ap, Bp, M);
		}            
	};

	// may throw SolverFailed or InconsistentSystem
	template <class Vector, class BB, class MyMethod> 
	Vector& solve(Vector& x, typename BB::Field::Element& d, const BB& A, const Vector& b, 
		      const RingCategories::IntegerTag & tag, 
		      const MyMethod& M)
	{
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Integer CRA Solve", "Isolve");

		RandomPrimeIterator genprime( 26 -(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		//         RationalRemainder< Modular<double> > rra((double)
		//                                                  ( A.coldim()/2.0*log((double) A.coldim()) ) );
	
		RationalRemainder< EarlyMultipRatCRA< Modular<double> > > rra(3UL);
		IntegerModularSolve<BB,Vector,MyMethod> iteration(A, b, M);

		// use of integer due to non genericity of rra (PG 2005-09-01)
		Integer den;
		std::vector< Integer > num(A.coldim());
		rra(num, den, iteration, genprime);
		//rra(x, d, iteration, genprime);
		
		typename Vector::iterator it_x= x.begin();
		typename std::vector<Integer>::const_iterator it_num= num.begin();
		
		// convert the result
		for (; it_x != x.end(); ++it_x, ++it_num)
			A.field().init(*it_x, *it_num);
		A.field().init(d, den);

		commentator.stop ("done", NULL, "Isolve");
		return x;
	}
	
	template <class RatVector, class Vector, class BB, class MyMethod> 
	RatVector& solve(RatVector& x, const BB& A, const Vector& b, 
			 const RingCategories::IntegerTag & tag, 
			 const MyMethod& M)
	{
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");

		commentator.start ("Rational CRA Solve", "Rsolve");
		typename BB::Field::Element den;
		std::vector<typename BB::Field::Element > num(A.coldim());
		solve (num, den, A, b, tag, M);
		typename RatVector::iterator it_x= x.begin();
		typename std::vector<typename BB::Field::Element>::const_iterator it_num= num.begin();
		integer n,d;
		A.field().convert(d,den); 
		for (; it_x != x.end(); ++it_x, ++it_num){			
			A.field().convert(n, *it_num);
			*it_x = typename RatVector::value_type(n, d);
		}
		commentator.stop ("done", NULL, "Rsolve");
		return x;
	}

        template <class RatVector, class Vector, class BB, class MethodTraits>
        RatVector& solve(RatVector& x, const BB& A, const Vector& b,
		              const RingCategories::RationalTag & tag,
		              const MethodTraits& m)
	{
		if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
			throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");
		commentator.start ("Rational CRA Solve", "Rsolve");
		size_t bits = 26 -(int)ceil(log((double)A.rowdim())*0.7213475205);
		RandomPrimeIterator genprime( bits);
	        RationalRemainder2< VarPrecEarlyMultipCRA< Modular<double> > > rra(3UL);//using default RR method
		IntegerModularSolve<BB,Vector,MethodTraits > iteration(A, b, m);
		integer den;
		std::vector< integer > num(A.coldim());
                rra(num, den, iteration, genprime);
                typename RatVector::iterator it_x= x.begin();
                typename std::vector<integer>::const_iterator it_num= num.begin();
                for (; it_x != x.end(); ++it_x, ++it_num){
			integer g = gcd( *it_num, den);
			*it_x = typename RatVector::value_type(*it_num/g, den/g);
		}
                commentator.stop ("done", NULL, "Rsolve");
                return x;
	}

	template <class RatVector, class BB, class MethodTraits>
	        RatVector& solve(RatVector& x, const BB& A, const RatVector& b,
		              const RingCategories::RationalTag & tag,
		              const MethodTraits& m)
		{
			if ((A.coldim() != x.size()) || (A.rowdim() != b.size()))
				throw LinboxError("LinBox ERROR: dimension of data are not compatible in system solving (solving impossible)");
			commentator.start ("Rational CRA Solve", "Rsolve");
			size_t bits = 26 -(int)ceil(log((double)A.rowdim())*0.7213475205);
			RandomPrimeIterator genprime( bits);
			RationalRemainder2< VarPrecEarlyMultipCRA< Modular<double> > > rra(3UL);//using default RR method
			IntegerModularSolve<BB,RatVector,MethodTraits > iteration(A, b, m);
			integer den;
			std::vector< integer > num(A.coldim());
			rra(num, den, iteration, genprime);
			typename RatVector::iterator it_x= x.begin();
			typename std::vector<integer>::const_iterator it_num= num.begin();
			for (; it_x != x.end(); ++it_x, ++it_num){
				integer g = gcd( *it_num, den);
				*it_x = typename RatVector::value_type(*it_num/g, den/g);
			}
			commentator.stop ("done", NULL, "Rsolve");
			return x;
		}
    
} // LinBox




#endif // __LINBOX_solve_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
