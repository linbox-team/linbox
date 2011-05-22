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

#ifndef __LINBOX_gauss_H
#define __LINBOX_gauss_H

#include "linbox/util/debug.h"
#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/sparse.h"
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
        const Field &field () const { return _F; }

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
                                                      SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///
        template <class Matrix> unsigned long& rankin(unsigned long &rank,
                                                      Matrix        &A,
                                                      unsigned long  Ni,
                                                      unsigned long  Nj,
						      SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
	template <class Matrix> unsigned long& rank(unsigned long &rank,
                                                    const Matrix        &A,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
        template <class Matrix> unsigned long& rank(unsigned long &rank,
                                                    const Matrix        &A,
                                                    unsigned long  Ni,
                                                    unsigned long  Nj,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            //@}

            /** @name det
                Callers of the different determinant routines\\
                -/ The "in" suffix indicates in place computation\\
                -/ Without Ni, Nj, the Matrix parameter must be a vector of sparse 
                row vectors, NOT storing any zero.\\
                -/ Calls {@link LinearPivoting } (by default) or {@link NoReordering}
            */
            //@{
            ///     
	template <class Matrix> Element& detin(Element &determinant,
                                               Matrix        &A,
                                               SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///
        template <class Matrix> Element& detin(Element &determinant,
                                               Matrix        &A,
                                               unsigned long  Ni,
                                               unsigned long  Nj,
                                               SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
	template <class Matrix> Element& det(Element &determinant,
                                             const Matrix        &A,
                                             SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
        template <class Matrix> Element& det(Element &determinant,
                                             const Matrix        &A,
                                             unsigned long  Ni,
                                             unsigned long  Nj,
                                             SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
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
                The computedet indicates whether the algorithm must compute the determionant as it goes

                @ref [Jean-Guillaume Dumas and Gilles Villard, 
                Computing the rank of sparse matrices over finite fields.
                In Ganzha et~al. CASC'2002, pages 47--62.]
            */
	template <class Matrix, class Perm>
	unsigned long& QLUPin(unsigned long &rank,
                              Element& determinant,
                              Perm          &Q,
                              Matrix	    &L,
                              Matrix        &U,
                              Perm	    &P,
                              unsigned long Ni, 
                              unsigned long Nj) const;

        template <class Matrix, class Perm, class Vector1, class Vector2> 
        Vector1& solve(Vector1& x, unsigned long rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b, bool randomsol=false)  const;
        
        template <class Matrix, class Perm, class Vector1, class Vector2> 
        Vector1& solve(Vector1& x, Vector1& w, unsigned long rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b)  const;
        

	template <class Matrix, class Vector1, class Vector2>
	Vector1& solvein(Vector1	&x,
                         Matrix         &A,
                         const Vector2	&b, bool randomsol=false)  const;


        template <class Matrix, class Perm, class Block> 
        Block& nullspacebasis(Block& x, 
                              unsigned long rank, 
                              const Matrix& U, 
                              const Perm& P)  const ;

        template <class Matrix, class Block> 
        Block& nullspacebasisin(Block& x, Matrix& A)  const;
        
        template <class Matrix, class Block> 
        Block& nullspacebasis(Block& x, const Matrix& A)  const;


            // Sparsest method
            //   erases elements while computing rank/det.
	template <class Matrix>
	unsigned long& InPlaceLinearPivoting(unsigned long &rank,
                                              Element& determinant,
                                              Matrix        &A,
                                              unsigned long Ni, 
                                              unsigned long Nj) const;

            // Same as the latter but keeps trace 
            //   of column permutations
            //   of remaining elements in the matrix
	template <class Matrix,class Perm>
	unsigned long& InPlaceLinearPivoting(unsigned long &rank,
                                              Element& determinant,
                                              Matrix        &A,
                                              Perm          &P,
                                              unsigned long Ni, 
                                              unsigned long Nj) const;


            /** \brief Sparse Gaussian elimination without reordering. 

                Gaussian elimination is done on a copy of the matrix.
                Using : SparseFindPivot
                eliminate

                Requirements : SLA is an array of sparse rows
                WARNING : NOT IN PLACE, THERE IS A COPY.
                Without reordering (Pivot is first non-zero in row)

            */
	template <class Matrix>
	unsigned long& NoReordering (unsigned long &rank, Element& determinant, Matrix &LigneA, unsigned long Ni, unsigned long Nj) const;

            /** \brief Dense in place LU factorization without reordering

                Using : FindPivot and LU
            */
	template <class Matrix>
	unsigned long &LUin (unsigned long &rank, Matrix &A) const;


            /** \brief Dense in place Gaussian elimination without reordering

                Using : FindPivot and LU
            */
	template <class Matrix>
	unsigned long &upperin (unsigned long &rank, Matrix &A) const;
        


    protected:
    
            //-----------------------------------------
            // Sparse elimination using a pivot row :
            // lc <-- lc - lc[k]/lp[0] * lp 
            // D is the number of elements per column
            //   it is updated and used for reordering
            // Vector is a vector of Pair (lin_pair.h)
            //-----------------------------------------
	template <class Vector, class D>
	void eliminate (Element             & headpivot,
                        Vector              &lignecourante,
			const Vector        &lignepivot,
			const unsigned long indcol,
			const long indpermut,
			const unsigned long npiv,
			D                   &columns) const;

	template <class Vector, class D>
	void eliminate (Vector              &lignecourante,
			const Vector        &lignepivot,
			const unsigned long &indcol,
			const long &indpermut,
			D                   &columns) const;

        template <class Vector>
        void permute (Vector              &lignecourante,
                      const unsigned long &indcol,
                      const long &indpermut) const;
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
			const long &indpermut) const;

            //-----------------------------------------
            // Dense elimination using a pivot row :
            // lc <-- lc - lc[k]/lp[0] * lp 
            // Computing only for k to n (and not 0 to n in LU)
            //-----------------------------------------
	template<class Vector>
	void Upper (Vector       &lignecur,
		    const Vector &lignepivot,
		    unsigned long indcol,
		    long indpermut) const;

            //-----------------------------------------
            // Dense elimination using a pivot row :
            // lc <-- lc - lc[k]/lp[0] * lp 
            //-----------------------------------------
	template <class Vector>
	void LU (Vector       &lignecur,
		 const Vector &lignepivot,
		 unsigned long indcol,
		 long indpermut) const;

            //------------------------------------------
            // Looking for a non-zero pivot in a row 
            // Using the column density for reordering
            // Pivot is chosen as to :
            // 1. Row density is minimum
            // 2. Column density is minimum for this row
            //------------------------------------------
	template <class Vector, class D>
	void SparseFindPivot (Vector &lignepivot, unsigned long &indcol, long &indpermut, D &columns, Element& determinant) const;
	
            //------------------------------------------
            // Looking for a non-zero pivot in a row 
            // No reordering
            //------------------------------------------
	template <class Vector>
	void SparseFindPivot (Vector &lignepivot, unsigned long &indcol, long &indpermut, Element& determinant) const;
	
            //------------------------------------------
            // Looking for a non-zero pivot in a row  
            // Dense search
            //------------------------------------------
	template <class Vector>
	void FindPivot (Vector &lignepivot, unsigned long &k, long &indpermut) const;

    };

} // namespace LinBox

#include "linbox/algorithms/gauss.inl"
#include "linbox/algorithms/gauss-pivot.inl"
#include "linbox/algorithms/gauss-elim.inl"
#include "linbox/algorithms/gauss-solve.inl"
#include "linbox/algorithms/gauss-nullspace.inl"
#include "linbox/algorithms/gauss-rank.inl"
#include "linbox/algorithms/gauss-det.inl"

#endif // __LINBOX_gauss_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
