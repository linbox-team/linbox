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
		METHOD_ELIMINATION, METHOD_WIEDEMANN, METHOD_LANCZOS
	};

	SolverTraits (Method method = METHOD_WIEDEMANN, bool precondition = true, size_t rank = 0)
		: _method (method), _precondition (precondition), _rank (rank) 
	{}

	Method method () const { return _method; }
	bool precondition () const { return _precondition; }
	size_t rank () const { return _rank; }

    private:
	Method _method;
	bool _precondition;
	size_t _rank;
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
