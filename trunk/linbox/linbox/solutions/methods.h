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

struct SolverTraits
{
	enum Method {
		ELIMINATION, WIEDEMANN, LANCZOS
	};

	enum SingularState {
		UNKNOWN, SINGULAR, NONSINGULAR
	};

	enum {
		RANK_UNKNOWN = 0
	};

	enum Preconditioner {
		NONE, BUTTERFLY, SPARSE, TOEPLITZ
	};

	SolverTraits (Method method = WIEDEMANN,
		      Preconditioner precond = SPARSE,
		      size_t rank = RANK_UNKNOWN,
		      SingularState singular = UNKNOWN,
		      bool checkResult = true,
		      bool certificate = true,
		      int maxTries = 10)
		: _method (method), _preconditioner (precond), _rank (rank), _singular (singular), _checkResult (checkResult),
		  _certificate (certificate), _maxTries (maxTries)
	{}

	SolverTraits (const char *str)
	{}

	Method         method ()         const { return _method; }
	Preconditioner preconditioner () const { return _preconditioner; }
	size_t         rank ()           const { return _rank; }
	SingularState  singular ()       const { return _singular; }
	bool           checkResult ()    const { return _checkResult; }
	bool           certificate ()    const { return _certificate; }
	int            maxTries ()       const { return _maxTries; }

	void method         (Method m)         { _method = m; }
	void preconditioner (Preconditioner p) { _preconditioner = p; }
	void rank           (size_t r)         { _rank = r; }
	void singular       (SingularState s)  { _singular = s; }
	void checkResult    (bool s)           { _checkResult = s; }
	void certificate    (bool s)           { _certificate = s; }
	void maxTries       (int n)            { _maxTries = n; }

    private:
	Method         _method;
	Preconditioner _preconditioner;
	size_t         _rank;
	SingularState  _singular;
	bool           _checkResult;
	bool           _certificate;
	int            _maxTries;
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

struct MethodTrait
{
	typedef WiedemannTraits   Wiedemann;
	typedef LanczosTraits     Lanczos;
	typedef EliminationTraits Elimination;
};
 
}

#endif
