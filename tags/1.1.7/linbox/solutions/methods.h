/* linbox/solutions/methods.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2003-02-03  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Reorganization: moved all the method-specific traits to the
 * corresponding structures, out of SolverTraits. Added a class
 * BlockLanczosTraits.
 * ------------------------------------
 * 2002-07-08  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Changed the name _DEFAULT_EarlyTerm_THRESHOLD_ to the more
 * standard-consistent DEFAULT_EARLY_TERM_THRESHOLD; changed the name
 * Early_Term_Threshold to earlyTermThreshold, also in keeping with the
 * standard.
 *
 * Added method traits for elimination and lanczos
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __LINBOX_method_H
#define __LINBOX_method_H

#ifndef DEFAULT_EARLY_TERM_THRESHOLD
#  define DEFAULT_EARLY_TERM_THRESHOLD 20
#endif

#include "linbox/blackbox/dense.h"
#ifdef __LINBOX_HAVE_MPI
#include "linbox/util/mpicpp.h"
#endif


namespace LinBox
{
// JGD 22.01.2007 : adapted from Lidzhade Fhulu's
    template <typename EnumT, typename BaseEnumT>
    struct InheritEnum {
        InheritEnum() {}
        InheritEnum(EnumT e) : enum_(e) {}
        InheritEnum(BaseEnumT e) : baseEnum_(e) {}
        explicit InheritEnum( int val ) : enum_(static_cast<EnumT>(val)) {}
        operator EnumT() const { return enum_; }
    private:
        union { 
            EnumT enum_;
            BaseEnumT baseEnum_;
        };
};



	struct Specifier {
		/** Whether the system is known to be singular or nonsingular */
		enum SingularState {
			SINGULARITY_UNKNOWN, SINGULAR, NONSINGULAR
		};

		/** Which preconditioner to use to ensure generic rank profile
		 *
		 * NO_PRECONDITIONER - Do not use any preconditioner
		 * BUTTERFLY - Use a butterfly network, see @ref{Butterfly}
		 * SPARSE - Use a sparse preconditioner, c.f. (Mulders 2000)
		 * TOEPLITZ - Use a Toeplitz preconditioner, c.f. (Kaltofen and Saunders
		 * 1991)
		 * SYMMETRIZE - Use A^T A (Lanczos only)
		 * PARTIAL_DIAGONAL - Use AD, where D is a random nonsingular diagonal
		 * matrix (Lanczos only)
		 * PARTIAL_DIAGONAL_SYMMETRIZE - Use A^T D A, where D is a random
		 * nonsingular diagonal matrix (Lanczos only)
		 * FULL_DIAGONAL - Use D_1 A^T D_2 A D_1, where D_1 and D_2 are random
		 * nonsingular diagonal matrices (Lanczos only)
		 * DENSE (Dixon use)
		 */
        
		enum Preconditioner {
			NO_PRECONDITIONER, BUTTERFLY, SPARSE, TOEPLITZ, SYMMETRIZE, PARTIAL_DIAGONAL, PARTIAL_DIAGONAL_SYMMETRIZE, FULL_DIAGONAL, DENSE
		};

                /** Other shapes : 
                 *  UNIMOD_UT -- unimodular upper triang. Toeplitz
                 *  UNIMOD_LT -- unimodular lower triang. Toeplitz
                 *  UNIMOD_UH -- unimodular upper triang. Hankel
                 *  UNIMOD_LH -- unimodular lower triang. Hankel
                 **/
            	enum BlackboxShape {
                    DIAGONAL = 15, HANKEL, UNIMOD_UT, UNIMOD_LT,  UNIMOD_UH,  UNIMOD_LH,  BLKVECTOR, TRI_SUP, TRI_INF
            	};

            	/** Shape of a Blackbox
                 *  Precontioner shapes and 
                 *  other blackbox shape are 
                 *  combined
                 **/ 
            	typedef InheritEnum<BlackboxShape, Preconditioner> Shape;

		/** Whether the rank of the system is known (otherwise its value) */
		enum {
			RANK_UNKNOWN = 0
		};

		/** Whether the system is known to be symmetric */
		enum {
			SYMMETRIC = true, NON_SYMMETRIC = false
		};
    
		/** Whether the probabilistic computation has to be certified Las-Vegas */
		enum {
			CERTIFY = true, DONT_CERTIFY = false
		};

		/** Linear-time pivoting or not for eliminations */
                enum PivotStrategy {
			PIVOT_LINEAR, PIVOT_NONE
		};
            
		Specifier ( ) 
			: _preconditioner(NO_PRECONDITIONER),
			  _rank(RANK_UNKNOWN),
			  _singular(SINGULARITY_UNKNOWN),
			  _symmetric(NON_SYMMETRIC),
			  _certificate(CERTIFY),
			  _maxTries(1),
			  _ett(DEFAULT_EARLY_TERM_THRESHOLD),
			  _blockingFactor(16),
			  _strategy(PIVOT_LINEAR),
                          _shape(SPARSE),
			  _provensuccessprobability( 0.0 )
#ifdef __LINBOX_HAVE_MPI
			  , _communicatorp( 0 )
#endif
		{}
  
		Specifier (const Specifier& s): 
			_preconditioner( s._preconditioner),
			_rank( s._rank),
			_singular( s._singular),
			_symmetric( s._symmetric),
			_certificate( s._certificate),
			_maxTries( s._maxTries),
			_ett( s._ett),
			_blockingFactor( s._blockingFactor),
			_strategy( s._strategy),
                        _shape( s._shape),
			_provensuccessprobability( s._provensuccessprobability)
#ifdef __LINBOX_HAVE_MPI
			, _communicatorp(s._communicatorp)
#endif
		{}

		/** Accessors
		 * 
		 * These functions just return the corresponding parameters from the
		 * structure
		 */

		Preconditioner	preconditioner ()	const { return _preconditioner; }
		size_t		rank ()			const { return _rank; }
		SingularState	singular ()		const { return _singular; }
		bool		symmetric ()		const { return _symmetric; }
		bool		certificate ()		const { return _certificate; }
		unsigned long	maxTries ()		const { return _maxTries; }
		unsigned long	earlyTermThreshold ()	const { return _ett; }
		unsigned long	blockingFactor ()	const { return _blockingFactor; }
		PivotStrategy	strategy ()		const { return _strategy; }
            	Shape		shape ()		const { return _shape; }
		double		trustability ()		const { return _provensuccessprobability; }
		bool		checkResult ()		const { return _checkResult; }
#ifdef __LINBOX_HAVE_MPI
		Communicator* communicatorp ()		const { return _communicatorp; }
#endif


		/** Manipulators
		 *
		 * These functions allow on-the-fly modification of a SolverTraits
		 * structure. Note that it is guaranteed that your SolverTraits
		 * structure will not be modified during @ref{solve}.
		 */

		void preconditioner (Preconditioner p) { _preconditioner = p; }
		void rank           (size_t r)         { _rank = r; }
		void singular       (SingularState s)  { _singular = s; }
		void symmetric      (bool s)           { _symmetric = s; }
		void certificate    (bool s)           { _certificate = s; }
		void maxTries       (unsigned long n)  { _maxTries = n; }
		void earlyTermThreshold (unsigned long e) { _ett = e; }
		void blockingFactor (unsigned long b)  { _blockingFactor = b; }
		void strategy (PivotStrategy strategy) { _strategy = strategy; }
		void shape          (Shape s)          { _shape = s; }
		void trustability   (double p)         { _provensuccessprobability = p; }
		void checkResult    (bool s)           { _checkResult = s; }
#ifdef __LINBOX_HAVE_MPI
		void communicatorp  (Communicator* cp) { _communicatorp = cp; }
#endif


	protected:
		Preconditioner _preconditioner;
		size_t         _rank;
		SingularState  _singular;
		bool           _symmetric;
		bool           _certificate;
		unsigned long  _maxTries;
		unsigned long  _ett;
		unsigned long  _blockingFactor;
		PivotStrategy  _strategy;
            	Shape          _shape;
		double         _provensuccessprobability;
		bool           _checkResult;
#ifdef __LINBOX_HAVE_MPI
		Communicator*   _communicatorp;
#endif
	};
    
	struct HybridSpecifier :public Specifier {
		HybridSpecifier(){};
		HybridSpecifier (const Specifier& m): Specifier(m){};
#ifdef __LINBOX_HAVE_MPI
		HybridSpecifier (Communicator& C): Specifier()
		{ _communicatorp = &C; };
#endif
	};
	struct BlackboxSpecifier :public Specifier {
		BlackboxSpecifier(){};
		BlackboxSpecifier (const Specifier& m): Specifier(m){};
	};
	struct EliminationSpecifier :public Specifier {
		EliminationSpecifier(){};
		EliminationSpecifier (const Specifier& m): Specifier(m){};
	};

	struct WiedemannTraits : public Specifier {
		/** Constructor
		 *
		 * @param precond Preconditioner to use, default is sparse
		 * @param rank Rank, if known; otherwise use RANK_UNKNOWN
		 * @param singular Whether the system is known to be singular or
		 * nonsingular; default is UNKNOWN
		 * @param symmetric True only if the system is symmetric. This improves
		 * performance somewhat, but will yield incorrect results if the system
		 * is not actually symmetric. Default is false.
		 * @param certificate True if the solver should attempt to find a
		 * certificate of inconsistency if it suspects the system to be
		 * inconsistent; default is true
		 * @param maxTries Maximum number of trials before giving up and
		 * returning a failure; default is 100
		 */
		WiedemannTraits (
				 bool           symmetric      = NON_SYMMETRIC,
				 unsigned long  thres          = DEFAULT_EARLY_TERM_THRESHOLD,
				 size_t         rank           = RANK_UNKNOWN,
				 Preconditioner preconditioner = SPARSE,
				 SingularState  singular       = SINGULARITY_UNKNOWN,
				 bool           certificate    = CERTIFY,
				 unsigned long  maxTries       = 100,
				 bool           checkResult    = true
				 )

		{ Specifier::_preconditioner = preconditioner;
			Specifier::_rank =(rank);
			Specifier::_singular =(singular);
			Specifier::_symmetric =(symmetric);
			Specifier::_certificate =(certificate);
			Specifier::_maxTries =(maxTries);
			Specifier::_ett =(thres);
			Specifier::_checkResult = (checkResult);
		}

		WiedemannTraits( const Specifier& S) :  Specifier(S) {}   
	};
    
	struct WiedemannExtensionTraits : public WiedemannTraits {
		WiedemannExtensionTraits (
				 bool           symmetric      = NON_SYMMETRIC,
				 unsigned long  thres          = DEFAULT_EARLY_TERM_THRESHOLD,
				 size_t         rank           = RANK_UNKNOWN,
				 Preconditioner preconditioner = SPARSE,
				 SingularState  singular       = SINGULARITY_UNKNOWN,
				 bool           certificate    = CERTIFY,
				 unsigned long  maxTries       = 100,
				 bool           checkResult    = true
				 )
                        : WiedemannTraits(symmetric,thres,rank,preconditioner,singular,certificate,maxTries,checkResult) 
                {}
		WiedemannExtensionTraits( const Specifier& S) :  WiedemannTraits(S) {}   
	};
    



	struct LanczosTraits : public Specifier {
		/** Constructor
		 *
		 * @param precond Preconditioner to use, default is sparse
		 * @param maxTries Maximum number of trials before giving up and
		 * returning a failure; default is 100
		 */
		LanczosTraits (Preconditioner preconditioner = FULL_DIAGONAL,
			       unsigned long maxTries       = 100)
		{ Specifier::_preconditioner =(preconditioner);
			Specifier::_maxTries =(maxTries);    
		}
		LanczosTraits( const Specifier& S) :  Specifier(S) {}   
	};

	struct BlockLanczosTraits : public Specifier {
		/** Constructor
		 *
		 * @param precond Preconditioner to use, default is sparse
		 * @param maxTries Maximum number of trials before giving up and
		 * returning a failure; default is 100
		 * @param blockingFactor Blocking factor to use
		 */
		BlockLanczosTraits (Preconditioner preconditioner = FULL_DIAGONAL,
				    unsigned long  maxTries       = 100,
				    int            blockingFactor = 16)
		{ Specifier::_preconditioner =(preconditioner);
            
			Specifier::_maxTries = (maxTries);
            
			Specifier::_blockingFactor = (blockingFactor);
		}
        
		BlockLanczosTraits( const Specifier& S) :  Specifier(S) {}   
	};
    
	struct SparseEliminationTraits  : public Specifier {
		/** Constructor
		 *
		 * @param strategy Pivoting strategy to use
		 */
		SparseEliminationTraits (PivotStrategy strategy = PIVOT_LINEAR) { Specifier::_strategy = (strategy) ;}
		SparseEliminationTraits( const EliminationSpecifier& S) :  Specifier(S) {}   
	};


    	struct DixonTraits : public Specifier {
		
		enum SolutionType {
			DETERMINIST, RANDOM, DIOPHANTINE 
		};

		DixonTraits ( SolutionType   solution       = DETERMINIST,
			      SingularState  singular       = SINGULARITY_UNKNOWN,
			      bool           certificate    = DONT_CERTIFY,
			      int            maxTries       = 10,
			      Preconditioner preconditioner = DENSE,
			      size_t          rank          = RANK_UNKNOWN)
		{ 
			_solution= (solution);
			Specifier::_singular= (singular);
			Specifier::_certificate= (certificate);
			Specifier::_maxTries= (maxTries);
			Specifier::_preconditioner=(preconditioner);		    
			Specifier::_rank=(rank);
		}

		DixonTraits( const Specifier& S) :  Specifier(S) {
			_solution= RANDOM;
		}   
		
          	SolutionType solution () const { return _solution;}
		
		void solution (SolutionType s) { _solution= (s);}

	protected:
		SolutionType _solution;
	};


	struct BlockWiedemannTraits : public Specifier {
		BlockWiedemannTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
				       size_t          rank            = RANK_UNKNOWN)
		{
			Specifier::_preconditioner = preconditioner;
			Specifier::_rank=rank;
		}
		BlockWiedemannTraits( const Specifier& S) :  Specifier(S) {}   
	};

	//Using numerical methods to symbolically solve linear systems. 
	//based on a preprinted article, submitted to JSC 2004
	struct NumericalTraits : public Specifier{
		NumericalTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
				  size_t          rank          = RANK_UNKNOWN)
		{ Specifier::_preconditioner=(preconditioner);
			
			Specifier::_rank=(rank) ;
		}
		NumericalTraits( const Specifier& S) :  Specifier(S) {}   
	};

	struct BlockHankelTraits : public Specifier {
		BlockHankelTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
				    size_t          rank            = RANK_UNKNOWN)
		{
			Specifier::_preconditioner = preconditioner;
			Specifier::_rank=rank;
		}
		BlockHankelTraits( const Specifier& S) :  Specifier(S) {}   
	};


	struct BlasEliminationTraits : public Specifier {
		BlasEliminationTraits() {}
		BlasEliminationTraits( const Specifier& S) :  Specifier(S) {}   

	};
    
	struct BlasExtensionTraits : public BlasEliminationTraits {
            	BlasExtensionTraits (bool           //certificate    = CERTIFY
									 , unsigned long  //maxTries       = 100
									 , bool         //checkResult    = true
                                     )
                        : BlasEliminationTraits() 
                {}
            	BlasExtensionTraits( const Specifier& S) :  BlasEliminationTraits(S) {}   
	};
    

	struct NonBlasEliminationTraits : public Specifier {
		NonBlasEliminationTraits() {}
		NonBlasEliminationTraits( const Specifier& S) :  Specifier(S) {}   
        
	};

	/// Method specifiers for controlling algorithm choice
	struct Method {
		typedef HybridSpecifier		Hybrid;
		typedef BlackboxSpecifier	Blackbox;
		typedef EliminationSpecifier	Elimination;
		typedef WiedemannTraits		Wiedemann;
		typedef WiedemannExtensionTraits ExtensionWiedemann;
		typedef LanczosTraits		Lanczos;
		typedef BlockLanczosTraits	BlockLanczos;
		typedef SparseEliminationTraits	SparseElimination;       
		typedef NumericalTraits		Numerical;
            	typedef BlasEliminationTraits 	BlasElimination;
		typedef BlasExtensionTraits ExtensionBlasElimination;
		typedef NonBlasEliminationTraits NonBlasElimination;
		typedef DixonTraits             Dixon;
		typedef BlockHankelTraits       BlockHankel;
		Method(){}
	};

	template<class BB>
	bool useBB(const BB& A) 
	{  return (A.coldim() > 1000) && (A.rowdim() > 1000);
	}

	template<class Field>
	bool useBB(const DenseMatrix<Field>& A) { return false; }

	/** Solver traits
	 *
         * User-specified parameters for solving a linear system.
         */
	struct SolverTraits : public Specifier
	{
		/** Constructor
		 *
		 * @param checkResult True if and only if the solution should be checked
		 * for correctness after it is computed (very much recommended for the
		 * randomized algorithms Wiedemann and Lanczos); default is true
		 */
		SolverTraits (bool checkResult = true)
		{                Specifier::_checkResult = checkResult;
		}

		/** Constructor from a MethodTraits structure
		 *
		 * @param traits MethodTraits structure from which to get defaults
		 */
		SolverTraits( const Specifier& S) :  Specifier(S) {}
	};

	/** Exception thrown when the computed solution vector is not a true
	 * solution to the system, but none of the problems cited below exist
	 */
	class SolveFailed {};

	/** Exception thrown when the system to be solved is
	 * inconsistent. Contains a certificate of inconsistency.
	 */
	template <class Vector>
	class InconsistentSystem 
	{
	public:
		InconsistentSystem (Vector &u)
			: _u (u)
		{}

		const Vector &certificate () const { return _u; }

	private:

		Vector _u;
	};

}

#endif // __LINBOX_method_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
