/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008 
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <21 Jun 10 14:43:11 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_gauss_nullspace_INL
#define __LINBOX_gauss_nullspace_INL

#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/vector/sparse.h"

namespace LinBox 
{

    

        // U is supposed full rank upper triangular
    template <class _Field>
    template <class Matrix, class Perm, class Block> inline Block& 
    GaussDomain<_Field>::nullspacebasis(Block& x, unsigned long rank, const Matrix& U, const Perm& P)  const {
        if (rank == 0) {
            for(size_t i=0; i<U.coldim(); ++i)
                x.setEntry(i,i,_F.one);
        } else {
            unsigned long nullity = U.coldim()-rank;
            if (nullity != 0) {
                    // compute U2T s.t. U = [ U1 | -U2T^T ]
                Matrix U2T(_F,nullity,rank);

                for(typename Matrix::ConstRawIndexedIterator uit=U.rawIndexedBegin();
                    uit != U.rawIndexedEnd(); ++uit) {
                    if (uit.colIndex() >= rank)
                        U2T.setEntry(uit.colIndex()-rank,uit.rowIndex(),uit.value());
                }
                for(typename Matrix::RawIterator u2it=U2T.rawBegin();
                    u2it != U2T.rawEnd(); ++u2it)
                    _F.negin(*u2it);

        
                    // Compute the basis vector by vector
                typedef Sparse_Vector< typename _Field::Element > SparseVect;
                for(size_t i=0; i<nullity; ++i) {
                    SparseVect W1Ti;
                        // Solve for upper part of basis
                    upperTriangularSparseSolve(W1Ti, rank, U, U2T[i]);
                        // Add identity for lower part
                    W1Ti.push_back( typename SparseVect::Element(rank+i, _F.one ) );

                    for(size_t j=0; j<W1Ti.size(); ++j) {
                            // P.applyTranspose(x[i],W1T[i]);
                            // Transposein(x)
                        x.setEntry( P.getStorage()[ W1Ti[j].first ], i, W1Ti[j].second );
                    }
                }
            }
        }
// x.write( std::cerr << "X:=", FORMAT_MAPLE ) << ';' << std::endl;
        return x;
    }

    template <class Matrix>
    inline bool nextnonzero(size_t& k, size_t Ni, const Matrix& A) {
        for(++k; k<Ni; ++k)
            if (A[k].size() > 0) return true;
        return false;
    }

        // Matrix A is upper triangularized
    template <class _Field>
    template <class Matrix, class Block> inline Block& 
    GaussDomain<_Field>::nullspacebasisin(Block& x, Matrix& A)  const {
        typename Field::Element det;
        unsigned long rank;
        size_t Ni(A.rowdim()),Nj(A.coldim());

        Permutation<Field> P(Nj,_F);

// A.write( std::cerr << "A:=", FORMAT_MAPLE ) << ';' << std::endl;
        this->InPlaceLinearPivoting(rank, det, A, P, Ni, Nj );

// P.write( std::cerr << "P:=", FORMAT_MAPLE ) << ';' << std::endl;
// A.write( std::cerr << "Ua:=", FORMAT_MAPLE ) << ';' << std::endl;

        for(size_t i=0; i< Ni; ++i) {
            if (A[i].size() == 0) {
                size_t j(i);
                if (nextnonzero(j,Ni,A)) {
                    A[i] = A[j];
                    A[j].resize(0);
                } else {
                    break;
                }
            }
        }

// A.write( std::cerr << "Ub:=", FORMAT_MAPLE ) << ';' << std::endl;

        return this->nullspacebasis(x, rank, A, P);
    }
    
    template <class _Field>
    template <class Matrix, class Block> inline Block& 
    GaussDomain<_Field>::nullspacebasis(Block& x, const Matrix& A)  const {
        SparseMatrix<Field, typename LinBox::Vector<Field>::SparseSeq> A1 (A); 
        return this->nullspacebasisin(x, A1);
    }
    
        
        

} // namespace LinBox

#endif // __LINBOX_gauss_nullspace_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
