/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/methods.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
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

namespace LinBox
{

/** Solver traits
 *
 * User-specified parameters for solving a linear system.
 */

struct SolverTraits
{
	/** Algorithm to use
	 *
	 * ELIMINATION - Gaussian elimination
	 * WIEDEMANN - Wiedemann's iterative algorithm
	 * LANCZOS - Lanczos iteration
	 * BLOCK_LANCZOS - Block Lanczos iteration
	 */

	enum Method {
		ELIMINATION, WIEDEMANN, LANCZOS, BLOCK_LANCZOS
	};

	/** Whether the system is known to be singular or nonsingular */

	enum SingularState {
		UNKNOWN, SINGULAR, NONSINGULAR
	};

	enum {
		RANK_UNKNOWN = 0
	};

	/** Which preconditioner to use to ensure generic rank profile
	 *
	 * NONE - Do not use any preconditioner
	 * DEFAULT - Use the default for the given solution method
	 * BUTTERFLY - Use a butterfly network, see @ref{Butterfly} (Wiedemann only)
	 * SPARSE - Use a sparse preconditioner, c.f. (Mulders 2000) (Wiedemann only)
	 * TOEPLITZ - Use a Toeplitz preconditioner, c.f. (Kaltofen and Saunders
	 * 1991) (Wiedemann only)
	 * SYMMETRIZE - Use A^T A (Lanczos only)
	 * PARTIAL_DIAGONAL - Use AD, where D is a random nonsingular diagonal
	 * matrix (Lanczos only)
	 * PARTIAL_DIAGONAL_SYMMETRIZE - Use A^T D A, where D is a random
	 * nonsingular diagonal matrix (Lanczos only)
	 * FULL_DIAGONAL - Use D_1 A^T D_2 A D_1, where D_1 and D_2 are random
	 * nonsingular diagonal matrices (Lanczos only)
	 */

	enum Preconditioner {
		NONE, DEFAULT, BUTTERFLY, SPARSE, TOEPLITZ, SYMMETRIZE,
		PARTIAL_DIAGONAL, PARTIAL_DIAGONAL_SYMMETRIZE, FULL_DIAGONAL
	};

	/** Constructor
	 *
	 * @param method Method to use, default is Wiedemann
	 * @param precond Preconditioner to use, default is sparse
	 * @param rank Rank, if known; otherwise use RANK_UNKNOWN
	 * @param singular Whether the system is known to be singular or
	 * nonsingular; default is UNKNOWN
	 * @param symmetric True only if the system is symmetric. This improves
	 * performance somewhat, but will yield incorrect results if the system
	 * is not actually symmetric. Default is false.
	 * @param checkResult True if and only if the solution should be checked
	 * for correctness after it is computed (very much recommended for the
	 * randomized algorithms Wiedemann and Lanczos); default is true
	 * @param certificate True if the solver should attempt to find a
	 * certificate of inconsistency if it suspects the system to be
	 * inconsistent; default is true
	 * @param maxTries Maximum number of trials before giving up and
	 * returning a failure; default is 100
	 * @param blockingFactor Blocking factor to use for blocked algorithms
	 */

	SolverTraits (Method method = WIEDEMANN,
		      Preconditioner precond = DEFAULT,
		      size_t rank = RANK_UNKNOWN,
		      SingularState singular = UNKNOWN,
		      bool symmetric = false,
		      bool checkResult = true,
		      bool certificate = true,
		      int maxTries = 100,
		      size_t blockingFactor = 4)
		: _method (method), _preconditioner (precond), _rank (rank), _singular (singular), _symmetric (symmetric),
		  _checkResult (checkResult), _certificate (certificate), _maxTries (maxTries),
		  _blockingFactor (blockingFactor)
	{}

	/** Accessors
	 * 
	 * These functions just return the corresponding parameters from the
	 * structure
	 */

	Method         method ()         const { return _method; }
	Preconditioner preconditioner () const { return _preconditioner; }
	size_t         rank ()           const { return _rank; }
	SingularState  singular ()       const { return _singular; }
	bool           symmetric ()      const { return _symmetric; }
	bool           checkResult ()    const { return _checkResult; }
	bool           certificate ()    const { return _certificate; }
	int            maxTries ()       const { return _maxTries; }
	size_t         blockingFactor () const { return _blockingFactor; }

	/** Manipulators
	 *
	 * These functions allow on-the-fly modification of a SolverTraits
	 * structure. Note that it is guaranteed that your SolverTraits
	 * structure will not be modified during @ref{solve}.
	 */

	void method         (Method m)         { _method = m; }
	void preconditioner (Preconditioner p) { _preconditioner = p; }
	void rank           (size_t r)         { _rank = r; }
	void singular       (SingularState s)  { _singular = s; }
	void symmetric      (bool s)           { _symmetric = s; }
	void checkResult    (bool s)           { _checkResult = s; }
	void certificate    (bool s)           { _certificate = s; }
	void maxTries       (int n)            { _maxTries = n; }
	void blockingFactor (size_t n)         { _blockingFactor = n; }

    private:
	Method         _method;
	Preconditioner _preconditioner;
	size_t         _rank;
	SingularState  _singular;
	bool           _symmetric;
	bool           _checkResult;
	bool           _certificate;
	int            _maxTries;
	size_t         _blockingFactor;
};

struct WiedemannTraits
{
	WiedemannTraits (unsigned long thres = DEFAULT_EARLY_TERM_THRESHOLD) : _ett (thres) {}
	unsigned long earlyTermThreshold () const { return _ett; }

    private:
	unsigned long _ett;
};

struct LanczosTraits
{
};

struct EliminationTraits
{
	enum PivotStrategy {
		PIVOT_PARTIAL, PIVOT_FULL
	};

	EliminationTraits (PivotStrategy strategy) : _strategy (strategy) {}
	PivotStrategy strategy () const { return _strategy; }

    private:
	PivotStrategy _strategy;
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

	const Vector &u () const { return _u; }

    private:

	Vector _u;
};

struct MethodTrait
{
	typedef WiedemannTraits   Wiedemann;
	typedef LanczosTraits     Lanczos;
	typedef EliminationTraits Elimination;
};
 
}

#endif
