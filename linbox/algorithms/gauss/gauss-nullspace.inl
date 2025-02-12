/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <26 Jan 24 16:05:55 Jean-Guillaume.Dumas@imag.fr>
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_gauss_nullspace_INL
#define __LINBOX_gauss_nullspace_INL

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/vector/sparse.h"

namespace LinBox
{



	// U is supposed full Rank upper triangular
	template <class _Field>
	template <class _Matrix, class Perm, class Block> inline Block&
	GaussDomain<_Field>::nullspacebasis(Block& x, size_t Rank, const _Matrix& U, const Perm& P)  const
	{
		if (Rank == 0) {
			for(size_t i=0; i<U.coldim(); ++i)
				x.setEntry(i,i,field().one);
		}
		else {
			size_t nullity = U.coldim()-Rank;
            x.resize(x.rowdim(),nullity);
			if (nullity != 0) {
				// compute U2T s.t. U = [ U1 | -U2T^T ]
				_Matrix U2T(field(),nullity,Rank);

				for(typename _Matrix::ConstIndexedIterator uit=U.IndexedBegin();
				    uit != U.IndexedEnd(); ++uit) {
					if (uit.colIndex() >= Rank)
						U2T.setEntry(uit.colIndex()-Rank,uit.rowIndex(),uit.value());
				}
				for(typename _Matrix::Iterator u2it=U2T.Begin();
				    u2it != U2T.End(); ++u2it)
					field().negin(*u2it);

				// Compute the basis vector by vector
				typedef Sparse_Vector< typename _Field::Element > SparseVect;
				for(size_t i=0; i<nullity; ++i) {
					SparseVect W1Ti;
					// Solve for upper part of basis
					upperTriangularSparseSolve(W1Ti, Rank, U, U2T[i]);
					// Add identity for lower part
					W1Ti.emplace_back((unsigned)(Rank+i), field().one );

					for(size_t j=0; j<W1Ti.size(); ++j) {
						// P.applyTranspose(x[i],W1T[i]);
						// Transposein(x)
						x.setEntry( (size_t)P.getStorage()[ (size_t)W1Ti[(size_t)j].first ], i, W1Ti[j].second );
					}
				}
			}
		}
		// x.write( std::cerr << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;
		return x;
	}

	template <class _Matrix>
	inline bool nextnonzero(size_t& k, size_t Ni, const _Matrix& A)
	{
		for(++k; k<Ni; ++k)
			if (A[k].size() > 0) return true;
		return false;
	}

	// _Matrix A is upper triangularized
	template <class _Field>
	template <class _Matrix, class Block> inline Block&
	GaussDomain<_Field>::nullspacebasisin(Block& x, _Matrix& A)  const
	{
		typename Field::Element Det;
		size_t Rank;
		size_t Ni(A.rowdim()),Nj(A.coldim());

		Permutation<Field> P(field(),(int)Nj);

// A.write( std::cerr << "A:=", Tag::FileFormat::Maple ) << ';' << std::endl;
		this->InPlaceLinearPivoting(Rank, Det, A, P, Ni, Nj );

// P.write( std::cerr << "P:=", Tag::FileFormat::Maple ) << ';' << std::endl;
// A.write( std::cerr << "Ua:=", Tag::FileFormat::Maple ) << ';' << std::endl;

		for(size_t i=0; i< Ni; ++i) {
			if (A[i].size() == 0) {
				size_t j(i);
				if (nextnonzero(j,Ni,A)) {
					A[i] = A[j];
					A[j].resize(0);
				}
				else {
					break;
				}
			}
		}

// A.write( std::cerr << "Ub:=", Tag::FileFormat::Maple ) << ';' << std::endl;

		return this->nullspacebasis(x, Rank, A, P);
	}

	template <class _Field>
	template <class _Matrix, class Block> inline Block&
	GaussDomain<_Field>::nullspacebasis(Block& x, const _Matrix& A)  const
	{
		Matrix A1 (A); // Must copy, then best to copy to preferred
		return this->nullspacebasisin(x, A1);
	}




} // namespace LinBox

#endif // __LINBOX_gauss_nullspace_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
