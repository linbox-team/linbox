/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/gauss.h
 * Copyright (C) 1999 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *
 * -----------------------------------------------------------
 * 2003-02-02  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Ported to new matrix archetype; update interface to meet current
 * standards. Rename gauss_foo as foo and gauss_Uin as upperin
 *
 * Move function definitions to gauss.inl
 * -----------------------------------------------------------
 *
 * See COPYING for license information.
 */

// ========================================================================= //
// (C) The Linbox Group 1999
// Calcul de rang par la méthode de Gauss pivot par ligne, sur matrice creuse
// Time-stamp: <03 Nov 00 19:19:06 Jean-Guillaume.Dumas@imag.fr> 
// ========================================================================= //

#ifndef __GAUSS_H
#define __GAUSS_H

#include "linbox/util/debug.h"
#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/archetype.h"

namespace LinBox 
{

/** @memo Repository of functions for rank by elimination 
 on sparse matrices.
 @doc Several versions allow for adjustment of the pivoting strategy
 and for choosing in-place elimination or for not modifying the input matrix.
 Also an LU interface is offered.
 */
template <class _Field>
class GaussDomain {
    public:
	typedef _Field Field;
	typedef typename Field::Element Element;

    private:
	const Field         &_F;

    public:

	/** @memo The field parameter is the domain 
	 * over which to perform computations
	 */
	GaussDomain (const Field &F) : _F (F) {}

	//Copy constructor
	/// 
	GaussDomain (const GaussDomain &M) : _F (M._F) {}

	/** accessor for the field of computation
	 */
	const Field &field () { return _F; }

    protected:
    
	//-----------------------------------------
	// Sparse elimination using a pivot row :
	// lc <-- lc - lc[k]/lp[0] * lp 
	// D is the number of elements per column
	//   it is updated and used for reordering
	// Vector is a vector of Pair (lin_pair.h)
	//-----------------------------------------
	template <class Vector, class D>
	void eliminate (Vector              &lignecourante,
			const Vector        &lignepivot,
			const unsigned long &indcol,
			const unsigned long &indpermut,
			D                   &columns);

	//-----------------------------------------
	// Sparse elimination using a pivot row :
	// lc <-- lc - lc[k]/lp[0] * lp 
	// No density update
	// Vector is a vector of Pair (lin_pair.h)
	//-----------------------------------------
	template <class Vector>
	void eliminate (Vector              &lignecourante,
			const Vector        &lignepivot,
			const unsigned long &indcol,
			const unsigned long &indpermut);

	//-----------------------------------------
	// Dense elimination using a pivot row :
	// lc <-- lc - lc[k]/lp[0] * lp 
	// Computing only for k to n (and not 0 to n in LU)
	//-----------------------------------------
	template<class Vector>
	void Upper (Vector       &lignecur,
		    const Vector &lignepivot,
		    unsigned long indcol,
		    unsigned long indpermut);

	//-----------------------------------------
	// Dense elimination using a pivot row :
	// lc <-- lc - lc[k]/lp[0] * lp 
	//-----------------------------------------
	template <class Vector>
	void LU (Vector       &lignecur,
		 const Vector &lignepivot,
		 unsigned long indcol,
		 unsigned long indpermut);

	//------------------------------------------
	// Looking for a non-zero pivot in a row 
	// Using the column density for reordering
	// Pivot is chosen as to :
	// 1. Row density is minimum
	// 2. Column density is minimum for this row
	//------------------------------------------
	template <class Vector, class D>
	void SparseFindPivot (Vector &lignepivot, unsigned long &indcol, unsigned long &indpermut, D &columns);

	//------------------------------------------
	// Looking for a non-zero pivot in a row 
	// No reordering
	//------------------------------------------
	template <class Vector>
	void SparseFindPivot (Vector &lignepivot, unsigned long &indcol, unsigned long &indpermut);

	//------------------------------------------
	// Looking for a non-zero pivot in a row  
	// Dense search
	//------------------------------------------
	template <class Vector>
	void FindPivot (Vector &lignepivot, unsigned long &k, unsigned long &indpermut);

    public:

	/** @memo
	 Sparse in place Gaussian elimination with reordering to reduce fill-in.
	 @doc pivots are chosen in sparsest column of sparsest row.
	 This runs in linear overhead.
	 It is similar in spirit but different from Markovitz' approach.  

	 \begin{verbatim}
	 Using : SparseFindPivot(..., density) for sparsest column, and 
	         eliminate (..., density)
	 \end{verbatim}

	 The Matrix parameter must meet the LinBox sparse matrix interface.
	 [check details].
	 The storrows indicates whether the algorithm must keep already computed rows.

	 @Ref paper.
	*/
	template <class Matrix>
	void rankinFullPivot (unsigned long &rank,
			      Matrix        &A,
			      bool           storrows = false);

	/** @memo
	 With this signature the Matrix parameter must be a vector of sparse 
	 row vectors.
	 @doc
	 Sparse in place Gaussian elimination with reordering to reduce fill-in.
	 */
	template <class Matrix>
	void rankinFullPivot (unsigned long &rank,
			      Matrix        &LigneA,
			      unsigned long  Ni,
			      unsigned long  Nj,
			      bool           storrows = false);

	/** @memo
	 Sparse Gaussian elimination with reordering to reduce fill-in.  Not in place.
	 @doc
	 Gaussian elimination is done on a copy of the matrix.
	 Using : SparseFindPivot
	         eliminate

	 Requirements : SLA is an array of sparse rows
	 WARNING : NOT IN PLACE, THERE IS A COPY.
	 Without reordering (Pivot is first non-zero in row)

	 */
	template <class Matrix>
	void rank (unsigned long &rank, const Matrix &SLA);

	template <class Matrix>
	void rankin (unsigned long &rank, Matrix &LigneA);

	template <class Matrix>
	void rankin (unsigned long &rank, Matrix &LigneA, unsigned long Ni, unsigned long Nj);

	/** @memo
	 Dense in place
	 Gaussian elimination without reordering
	 @doc
	 Using : FindPivot and LU
	 */
	template <class Matrix>
	long &upperin (unsigned long &rank, Matrix &A);

	/**
	 Dense in place
	 LU factorization without reordering
	 @doc Using : FindPivot and LU
	 */
	template <class Matrix>
	long &LUin (unsigned long &rank, Matrix &A);
};

} // namespace LinBox

#include "linbox/algorithms/gauss.inl"

#endif // __GAUSS_H
