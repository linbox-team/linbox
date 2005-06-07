/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __METHODS_H
#define __METHODS_H

#ifndef DEFAULT_EARLY_TERM_THRESHOLD
#  define DEFAULT_EARLY_TERM_THRESHOLD 20
#endif

#include "linbox/blackbox/dense.h"

namespace LinBox
{

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

	enum PivotStrategy {
            PIVOT_LINEAR, PIVOT_NONE
	};


	Specifier ( ) 
                : _preconditioner(NO_PRECONDITIONER),
                  _rank(RANK_UNKNOWN),
                  _singular(SINGULARITY_UNKNOWN),
                  _symmetric(NON_SYMMETRIC),
                  _certificate(DONT_CERTIFY),
                  _maxTries(1),
                  _ett(DEFAULT_EARLY_TERM_THRESHOLD),
                  _blockingFactor(16),
                  _strategy(PIVOT_LINEAR),
                  _provensuccessprobability( 0.0 )
            {}
  
            /** Accessors
             * 
             * These functions just return the corresponding parameters from the
             * structure
             */

	Preconditioner preconditioner ()     const { return _preconditioner; }
	size_t         rank ()               const { return _rank; }
	SingularState  singular ()           const { return _singular; }
	bool           symmetric ()          const { return _symmetric; }
	bool           certificate ()        const { return _certificate; }
	int            maxTries ()           const { return _maxTries; }
	unsigned long  earlyTermThreshold () const { return _ett; }
        unsigned long  blockingFactor () const { return _blockingFactor; }
	PivotStrategy strategy () const { return _strategy; }
        double trustability ()   const  { return _provensuccessprobability; }


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
	void maxTries       (int n)            { _maxTries = n; }
        void blockingFactor (unsigned long b)  { _blockingFactor = b; }
	void strategy (PivotStrategy strategy) { _strategy = strategy; }
        void trustability   (double p)         { _provensuccessprobability = p; }

    protected:
	Preconditioner _preconditioner;
	size_t         _rank;
	SingularState  _singular;
	bool           _symmetric;
	bool           _certificate;
	int            _maxTries;
	unsigned long  _ett;
	unsigned long  _blockingFactor;
	PivotStrategy _strategy;
        double         _provensuccessprobability;
    };
    
    struct HybridSpecifier {
		HybridSpecifier(){};
		HybridSpecifier(const Specifier& m): _m(m){};
		const Specifier& specifier() const { return _m; } 
		const Specifier _m;
    };
    struct BlackboxSpecifier {
		BlackboxSpecifier(){};
		BlackboxSpecifier (const Specifier& m): _m(m){};
		const Specifier& specifier() const { return _m; }
		const Specifier _m;
    };
    struct EliminationSpecifier {
		EliminationSpecifier(){};
		EliminationSpecifier (const Specifier& m): _m(m){};
		const Specifier& specifier() const { return _m; }
		const Specifier _m;
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
            int            maxTries       = 100)

            { Specifier::_preconditioner = preconditioner;
            Specifier::_rank =(rank);
            Specifier::_singular =(singular);
            Specifier::_symmetric =(symmetric);
            Specifier::_certificate =(certificate);
            Specifier::_maxTries =(maxTries);
            Specifier::_ett =(thres);
            }
        
    };
    
    struct LanczosTraits : public Specifier {
            /** Constructor
             *
             * @param precond Preconditioner to use, default is sparse
             * @param maxTries Maximum number of trials before giving up and
             * returning a failure; default is 100
             */
	LanczosTraits (Preconditioner preconditioner = FULL_DIAGONAL,
		       int            maxTries       = 100)
            { Specifier::_preconditioner =(preconditioner);
            Specifier::_maxTries =(maxTries);    
            }
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
			    int            maxTries       = 100,
			    int            blockingFactor = 16)
            { Specifier::_preconditioner =(preconditioner);
            
            Specifier::_maxTries = (maxTries);
            
            Specifier::_blockingFactor = (blockingFactor);
            }
        
    };
    
    struct SparseEliminationTraits  : public Specifier {
            /** Constructor
             *
             * @param strategy Pivoting strategy to use
             */
	SparseEliminationTraits (PivotStrategy strategy = PIVOT_LINEAR) 
            { Specifier::_strategy = (strategy) ;}
    };


    struct DixonTraits : public Specifier {
	DixonTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
		      size_t          rank          = RANK_UNKNOWN)
            { Specifier::_preconditioner=(preconditioner);
            
            Specifier::_rank=(rank);
            }
    };

    struct BlockWiedemannTraits : public Specifier {
	BlockWiedemannTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
			       size_t          rank            = RANK_UNKNOWN)
            {
                Specifier::_preconditioner = preconditioner;
                Specifier::_rank=rank;
            }
    };

	//Using numerical methods to symbolically solve linear systems. 
	//based on a preprinted article, submitted to JSC 2004
    struct NumericalTraits : public Specifier{
	NumericalTraits ( Preconditioner preconditioner = NO_PRECONDITIONER,
                          size_t          rank          = RANK_UNKNOWN)
            { Specifier::_preconditioner=(preconditioner);
            
            Specifier::_rank=(rank) ;
            }
    };

    struct BlasEliminationTraits : public Specifier {};
    struct NonBlasEliminationSpecifier : public Specifier {};


	/// Method specifiers for controlling algorithm choice
    struct Method {
	typedef HybridSpecifier		Hybrid;
	typedef BlackboxSpecifier	Blackbox;
	typedef EliminationSpecifier	Elimination;
        typedef WiedemannTraits		Wiedemann;
        typedef LanczosTraits		Lanczos;
        typedef BlockLanczosTraits	BlockLanczos;
        typedef SparseEliminationTraits	SparseElimination;       
        typedef NumericalTraits		Numerical;
        typedef BlasEliminationTraits 	BlasElimination;
        typedef NonBlasEliminationSpecifier NonBlasElimination;
	Method(){}
    };

	template<class BB>
	bool useBB(const BB& A)
	{  return A.coldim() > 1000 && A.rowdim() > 1000;
	}

	template<class Field>
	bool useBB(const DenseMatrix<Field>& A) { return false; }

	/** Solver traits
	 *
         * User-specified parameters for solving a linear system.
         */
    template <class MethodTraits>
    struct SolverTraits : public MethodTraits
    {
            /** Constructor
             *
             * @param checkResult True if and only if the solution should be checked
             * for correctness after it is computed (very much recommended for the
             * randomized algorithms Wiedemann and Lanczos); default is true
             */

	SolverTraits (bool checkResult = true)
		: _checkResult (checkResult)
            {}

            /** Constructor from a MethodTraits structure
             *
             * @param traits MethodTraits structure from which to get defaults
             * @param checkResult True if and only if the solution should be checked
             * for correctness after it is computed (very much recommended for the
             * randomized algorithms Wiedemann and Lanczos); default is true
             */
	SolverTraits (MethodTraits traits, bool checkResult = true)
		: MethodTraits (traits), _checkResult (checkResult)
            {}

            /** Accessors
             * 
             * These functions just return the corresponding parameters from the
             * structure
             */

	bool           checkResult ()    const { return _checkResult; }

            /** Manipulators
             *
             * These functions allow on-the-fly modification of a SolverTraits
             * structure. Note that it is guaranteed that your SolverTraits
             * structure will not be modified during @ref{solve}.
             */

	void checkResult    (bool s)           { _checkResult = s; }

    private:
	bool           _checkResult;
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

#endif
