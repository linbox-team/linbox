/* linbox/algorithms/gauss-gf2.h
 * Copyright (C) 2009 The LinBox group
 * Written by JG Dumas
 *
 * Time-stamp: <14 Jun 10 15:36:56 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 * SparseSeqMatrix is container< container< size_t > > 
 * as e.g. linbox/blackbox/zo-gf2.h
 */
// ========================================================================= //

#ifndef __LINBOX_gauss_gf2_H
#define __LINBOX_gauss_gf2_H

#include "linbox/util/debug.h"
#include "linbox/util/commentator.h"
#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/blackbox/zo-gf2.h"

namespace LinBox 
{

    template <>
    class GaussDomain<GF2> {
    public:
	typedef GF2 Field;
	typedef Field::Element Element;

	 // Preferred Matrix type
	typedef ZeroOne<GF2> Matrix;

    public:

            /** \brief The field parameter is the domain 
             * over which to perform computations
             */
	GaussDomain (const Field &) {}

            //Copy constructor
            /// 
	GaussDomain (const GaussDomain &)  {}

            /** accessor for the field of computation
             */
        const Field &field () const { return *(new GF2()); }

            /** @name rank
                Callers of the different rank routines\\
                -/ The "in" suffix indicates in place computation\\
                -/ Without Ni, Nj, the SparseSeqMatrix parameter must be a vector of sparse 
                row vectors, NOT storing any zero.\\
                -/ Calls {@link rankinLinearPivoting rankinLinearPivoting} (by default) or {@link rankinNoReordering rankinNoReordering}
            */
            //@{
            ///     
            ///
        template <class SparseSeqMatrix> unsigned long& rankin(unsigned long &rank,
                                                      SparseSeqMatrix        &A,
                                                      unsigned long  Ni,
                                                      unsigned long  Nj,
						      SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const ;

            ///        
    template <class SparseSeqMatrix> unsigned long& rankin(unsigned long &rank,
                                		  SparseSeqMatrix        &A,
                                		  SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;

            ///        
        template <class SparseSeqMatrix> unsigned long& rank(unsigned long &rk,
                                                    const SparseSeqMatrix        &A,
                                                    unsigned long  Ni,
                                                    unsigned long  Nj,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const ;

            ///        
	template <class SparseSeqMatrix> unsigned long& rank(unsigned long &rk,
                                                    const SparseSeqMatrix        &A,
                                                    SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const ;

            //@}

            /** @name det
                Callers of the different determinant routines\\
                -/ The "in" suffix indicates in place computation\\
                -/ Without Ni, Nj, the SparseSeqMatrix parameter must be a vector of sparse 
                row vectors, NOT storing any zero.\\
                -/ Calls {@link LinearPivoting } (by default) or {@link NoReordering}
            */
            //@{
            ///     
	template <class SparseSeqMatrix> Element& detin(Element &determinant,
                                               SparseSeqMatrix        &A,
                                               SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///
        template <class SparseSeqMatrix> Element& detin(Element &determinant,
                                               SparseSeqMatrix        &A,
                                               unsigned long  Ni,
                                               unsigned long  Nj,
                                               SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
	template <class SparseSeqMatrix> Element& det(Element &determinant,
                                             const SparseSeqMatrix        &A,
                                             SparseEliminationTraits::PivotStrategy   reord = SparseEliminationTraits::PIVOT_LINEAR) const;
            ///        
        template <class SparseSeqMatrix> Element& det(Element &determinant,
                                             const SparseSeqMatrix        &A,
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

                The SparseSeqMatrix parameter must meet the LinBox sparse matrix interface.
                [check details].
                The computedet indicates whether the algorithm must compute the determionant as it goes

                @ref [Jean-Guillaume Dumas and Gilles Villard, 
                Computing the rank of sparse matrices over finite fields.
                In Ganzha et~al. CASC'2002, pages 47--62.]
            */
	template <class SparseSeqMatrix, class Perm>
	unsigned long& QLUPin(unsigned long &rank,
                              Element& determinant,
                              Perm          &Q,
                              SparseSeqMatrix	    &L,
                              SparseSeqMatrix        &U,
                              Perm	    &P,
                              unsigned long Ni, 
                              unsigned long Nj) const;

        template <class SparseSeqMatrix, class Perm, class Vector1, class Vector2> 
        Vector1& solve(Vector1& x, unsigned long rank, const Perm& Q, const SparseSeqMatrix& L, const SparseSeqMatrix& U, const Perm& P, const Vector2& b, bool randomsol=false) const;
        

        template <class SparseSeqMatrix, class Perm, class Vector1, class Vector2> 
        Vector1& solve(Vector1& x, Vector1& w, unsigned long rank, const Perm& Q, const SparseSeqMatrix& L, const SparseSeqMatrix& U, const Perm& P, const Vector2& b) const;
        

	template <class SparseSeqMatrix, class Vector1, class Vector2>
	Vector1& solvein(Vector1& x,
                         SparseSeqMatrix        &A,
                         const Vector2& b, bool randomsol=false) const;


	template <class SparseSeqMatrix, class Perm>
	unsigned long& InPlaceLinearPivoting(unsigned long &rank,
                                              Element& determinant,
                                              SparseSeqMatrix        &A,
                                              Perm                   &P,
                                              unsigned long Ni, 
                                              unsigned long Nj) const;
	template <class SparseSeqMatrix>
        unsigned long& NoReordering (unsigned long & rank, Element& , SparseSeqMatrix &, unsigned long , unsigned long ) const {
		std::cerr << "Sparse elimination over GF2 without reordering not implemented" << std::endl;
		return rank;
	}


    protected:
    
            //-----------------------------------------
            // Sparse elimination using a pivot row :
            // lc <-- lc - lc[k]/lp[0] * lp 
            // D is the number of elements per column
            //   it is updated and used for reordering
            // Vector is a vector of Pair (lin_pair.h)
            //-----------------------------------------
	template <class Vector, class D>
	void eliminateBinary (Element             & headpivot,
                        Vector              &lignecourante,
			const Vector        &lignepivot,
			const unsigned long indcol,
			const long indpermut,
			const unsigned long npiv,
			D                   &columns) const;

       template <class Vector>
       void permuteBinary (Vector              &lignecourante,
                      const unsigned long &indcol,
                      const long &indpermut) const;
 
            //------------------------------------------
            // Looking for a non-zero pivot in a row 
            // Using the column density for reordering
            // Pivot is chosen as to :
            // 1. Row density is minimum
            // 2. Column density is minimum for this row
            //------------------------------------------
	template <class Vector, class D>
	void SparseFindPivotBinary (Vector &lignepivot, unsigned long &indcol, long &indpermut, D &columns, Element& determinant) const;
	
            //------------------------------------------
            // Looking for a non-zero pivot in a row 
            // No reordering
            //------------------------------------------
	template <class Vector>
	void SparseFindPivotBinary (Vector &lignepivot, unsigned long &indcol, long &indpermut, Element& determinant) const;
	
    };
} // namespace LinBox

#include "linbox/algorithms/gauss-gf2.inl"
#include "linbox/algorithms/gauss-pivot-gf2.inl"
#include "linbox/algorithms/gauss-elim-gf2.inl"
#include "linbox/algorithms/gauss-rank-gf2.inl"
#include "linbox/algorithms/gauss-solve-gf2.inl"

#endif // __LINBOX_gauss_gf2_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
