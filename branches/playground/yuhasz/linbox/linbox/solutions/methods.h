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

namespace LinBox
{

struct WiedemannTraits
{
	/** Whether the system is known to be singular or nonsingular */

	enum SingularState {
		UNKNOWN, SINGULAR, NONSINGULAR
	};

	/** Which preconditioner to use to ensure generic rank profile
	 *
	 * NONE - Do not use any preconditioner
	 * BUTTERFLY - Use a butterfly network, see @ref{Butterfly}
	 * SPARSE - Use a sparse preconditioner, c.f. (Mulders 2000)
	 * TOEPLITZ - Use a Toeplitz preconditioner, c.f. (Kaltofen and Saunders
	 * 1991)
	 */

	enum Preconditioner {
		NONE, BUTTERFLY, SPARSE, TOEPLITZ
	};

	enum {
		RANK_UNKNOWN = 0
	};

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
	WiedemannTraits (Preconditioner preconditioner = SPARSE,
			 size_t         rank           = RANK_UNKNOWN,
			 SingularState  singular       = UNKNOWN,
			 bool           symmetric      = false,
			 bool           certificate    = true,
			 int            maxTries       = 100,
			 unsigned long  thres          = DEFAULT_EARLY_TERM_THRESHOLD)
		: _preconditioner (preconditioner),
		  _rank (rank),
		  _singular (singular),
		  _symmetric (symmetric),
		  _certificate (certificate),
		  _maxTries (maxTries),
		  _ett (thres) {}

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

    private:
	Preconditioner _preconditioner;
	size_t         _rank;
	SingularState  _singular;
	bool           _symmetric;
	bool           _certificate;
	int            _maxTries;
	unsigned long  _ett;
};

struct LanczosTraits
{
	/** Which preconditioner to use to ensure generic rank profile
	 *
	 * NONE - Do not use any preconditioner
	 * SYMMETRIZE - Use A^T A (Lanczos only)
	 * PARTIAL_DIAGONAL - Use AD, where D is a random nonsingular diagonal
	 * matrix (Lanczos only)
	 * PARTIAL_DIAGONAL_SYMMETRIZE - Use A^T D A, where D is a random
	 * nonsingular diagonal matrix (Lanczos only)
	 * FULL_DIAGONAL - Use D_1 A^T D_2 A D_1, where D_1 and D_2 are random
	 * nonsingular diagonal matrices (Lanczos only)
	 */

	enum Preconditioner {
		NONE, SYMMETRIZE, PARTIAL_DIAGONAL, PARTIAL_DIAGONAL_SYMMETRIZE, FULL_DIAGONAL
	};

	/** Constructor
	 *
	 * @param precond Preconditioner to use, default is sparse
	 * @param maxTries Maximum number of trials before giving up and
	 * returning a failure; default is 100
	 */
	LanczosTraits (Preconditioner preconditioner = FULL_DIAGONAL,
		       int            maxTries       = 100)
		: _preconditioner (preconditioner),
		  _maxTries (maxTries) {}

	/** Accessors
	 * 
	 * These functions just return the corresponding parameters from the
	 * structure
	 */

	Preconditioner preconditioner () const { return _preconditioner; }
	unsigned int   maxTries ()       const { return _maxTries; }

	/** Manipulators
	 *
	 * These functions allow on-the-fly modification of a SolverTraits
	 * structure. Note that it is guaranteed that your SolverTraits
	 * structure will not be modified during @ref{solve}.
	 */

	void preconditioner (Preconditioner p) { _preconditioner = p; }
	void maxTries       (unsigned int m)   { _maxTries = m; }

    private:
	Preconditioner _preconditioner;
	unsigned int   _maxTries;
};

struct BlockLanczosTraits
{
	/** Which preconditioner to use to ensure generic rank profile
	 *
	 * NONE - Do not use any preconditioner
	 * SYMMETRIZE - Use A^T A (Lanczos only)
	 * PARTIAL_DIAGONAL - Use AD, where D is a random nonsingular diagonal
	 * matrix (Lanczos only)
	 * PARTIAL_DIAGONAL_SYMMETRIZE - Use A^T D A, where D is a random
	 * nonsingular diagonal matrix (Lanczos only)
	 * FULL_DIAGONAL - Use D_1 A^T D_2 A D_1, where D_1 and D_2 are random
	 * nonsingular diagonal matrices (Lanczos only)
	 */

	enum Preconditioner {
		NONE, SYMMETRIZE, PARTIAL_DIAGONAL, PARTIAL_DIAGONAL_SYMMETRIZE, FULL_DIAGONAL
	};

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
		: _preconditioner (preconditioner),
		  _maxTries (maxTries),
		  _blockingFactor (blockingFactor) {}

	/** Accessors
	 * 
	 * These functions just return the corresponding parameters from the
	 * structure
	 */

	Preconditioner preconditioner () const { return _preconditioner; }
	unsigned long  blockingFactor () const { return _blockingFactor; }
	unsigned int   maxTries ()       const { return _maxTries; }

	/** Manipulators
	 *
	 * These functions allow on-the-fly modification of a SolverTraits
	 * structure. Note that it is guaranteed that your SolverTraits
	 * structure will not be modified during @ref{solve}.
	 */

	void preconditioner (Preconditioner p) { _preconditioner = p; }
	void blockingFactor (unsigned long b)  { _blockingFactor = b; }
	void maxTries       (unsigned int m)   { _maxTries = m; }

    private:
	Preconditioner _preconditioner;
	unsigned int   _maxTries;
	unsigned long  _blockingFactor;
};

struct EliminationTraits
{
	enum PivotStrategy {
		PIVOT_LINEAR, PIVOT_NONE
	};

	/** Constructor
	 *
	 * @param strategy Pivoting strategy to use
	 */
	EliminationTraits (PivotStrategy strategy = PIVOT_LINEAR) : _strategy (strategy) {}

	/** Accessors
	 * 
	 * These functions just return the corresponding parameters from the
	 * structure
	 */

	PivotStrategy strategy () const { return _strategy; }

	/** Manipulators
	 *
	 * These functions allow on-the-fly modification of a SolverTraits
	 * structure. Note that it is guaranteed that your SolverTraits
	 * structure will not be modified during @ref{solve}.
	 */

	void strategy (PivotStrategy strategy) { _strategy = strategy; }

    private:
	PivotStrategy _strategy;
};

struct MethodTrait
{
	typedef WiedemannTraits    Wiedemann;
	typedef LanczosTraits      Lanczos;
	typedef BlockLanczosTraits BlockLanczos;
	typedef EliminationTraits  Elimination;
};

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

	const Vector &u () const { return _u; }

    private:

	Vector _u;
};

}

#endif
