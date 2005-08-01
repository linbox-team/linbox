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
#include "linbox/solutions/methods.h"

namespace LinBox 
{

/** \brief Repository of functions for rank by elimination on sparse matrices.

    Several versions allow for adjustment of the pivoting strategy
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

            /** \brief The field parameter is the domain 
             * over which to perform computations
             */
	GaussDomain (const Field &F) : _F (F) {}

            //Copy constructor
            /// 
	GaussDomain (const GaussDomain &M) : _F (M._F) {}

            /** accessor for the field of computation
             */
        const Field &field () { return _F; }

            /** @name rank
                Callers of the different rank routines\\
                -/ The "in" suffix indicates in place computation\\
                -/ Without Ni, Nj, the Matrix parameter must be a vector of sparse 
                row vectors, NOT storing any zero.\\
                -/ Calls {@link rankinLinearPivoting rankinLinearPivoting} (by default) or {@link rankinNoReordering rankinNoReordering}
            */
            //@{
            ///     
	template <class Matrix> unsigned long& rankin(unsigned long &rank,
                                                      Matrix        &A,
                                                      SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR,
                                                      bool           storrows = false);
            ///
        template <class Matrix> unsigned long& rankin(unsigned long &rank,
                                                      Matrix        &A,
                                                      unsigned long  Ni,
                                                      unsigned long  Nj,
						      SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR,
                                                      bool           storrows = false);
            ///        
	template <class Matrix> unsigned long& rank(unsigned long &rank,
                                                    const Matrix        &A,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR,
                                                    bool           storrows = false);
            ///        
        template <class Matrix> unsigned long& rank(unsigned long &rank,
                                                    const Matrix        &A,
                                                    unsigned long  Ni,
                                                    unsigned long  Nj,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR,
                                                    bool           storrows = false);
            //@}


            /** \brief Sparse in place Gaussian elimination with reordering to reduce fill-in.
                pivots are chosen in sparsest column of sparsest row.
                This runs in linear overhead.
                It is similar in spirit but different from Markovitz' approach.  

                <pre>
                Using : SparseFindPivot(..., density) for sparsest column, and 
                eliminate (..., density)
                </pre>

                The Matrix parameter must meet the LinBox sparse matrix interface.
                [check details].
                The storrows indicates whether the algorithm must keep already computed rows.

                @ref [Jean-Guillaume Dumas and Gilles Villard, 
                Computing the rank of sparse matrices over finite fields.
                In Ganzha et~al. CASC'2002, pages 47--62.]
            */
	template <class Matrix>
	unsigned long& rankinLinearPivoting (unsigned long &rank,
                                             Matrix        &A,
                                             unsigned long Ni, 
                                             unsigned long Nj,
                                             bool           storrows = false);


            /** \brief Sparse Gaussian elimination without reordering. 

                Gaussian elimination is done on a copy of the matrix.
                Using : SparseFindPivot
                eliminate

                Requirements : SLA is an array of sparse rows
                WARNING : NOT IN PLACE, THERE IS A COPY.
                Without reordering (Pivot is first non-zero in row)

            */
	template <class Matrix>
	unsigned long& rankinNoReordering (unsigned long &rank, Matrix &LigneA, unsigned long Ni, unsigned long Nj);

            /** \brief Dense in place LU factorization without reordering

                Using : FindPivot and LU
            */
	template <class Matrix>
	unsigned long &LUin (unsigned long &rank, Matrix &A);


            /** \brief Dense in place Gaussian elimination without reordering

                Using : FindPivot and LU
            */
	template <class Matrix>
	unsigned long &upperin (unsigned long &rank, Matrix &A);
        


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

    };

} // namespace LinBox

#include "linbox/algorithms/gauss.inl"

#endif // __GAUSS_H
