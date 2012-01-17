
/* linbox/blackbox/zero-one.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong, -bds
 * Time-stamp: <22 Jun 10 17:34:01 Jean-Guillaume.Dumas@imag.fr>
 *
 * ------------------------------------
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

#ifndef __LINBOX_zero_one_H
#define __LINBOX_zero_one_H

#include "linbox/integer.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/blackbox-interface.h"

// For STL pair in IndexIterator
#include <utility>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>
#include <iostream>

namespace LinBox
{

	/** \brief Time and space efficient representation of sparse {0,1}-matrices.
	 *
	 * A 0-1 matrix is a matrix with all 0's and 1's as entries.
	 * We're using a NAG-sparse format.
	 * Applies can be performed fast, using only additions.
	 * When initalizing this class, you only need to build 2 arrays of equal length:
	 * an array of the row indices for the non-zero (1's) entries, and an array of the column
	 * indices for the non-zero (1's) entries.

	 A {0, 1,-1} matrix can be effecively represented as the \ref Dif of two ZeroOne's.
	 \ingroup blackbox
	 */
	template<class _Field>
	class ZeroOne : public BlackboxInterface {
	protected:
		typedef size_t Index;
	public:
		typedef ZeroOne<_Field> Self_t;
		typedef _Field Field;
		typedef typename _Field::Element Element;

		// Default constructor, do nothing.
		ZeroOne(const Field& F);
		// The real constructor /todo give docs here
		ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort = false, bool colSort = false);
		// Destructor, once again do nothing
		~ZeroOne();

		/** apply.
		 *
		 * Uses one of the three
		 * private utility functions. It calls the generalized utility function
		 * _apply if there is no special ordering, _fyapply if there is C_ordering
		 * or _fxapply if there is fortran_ordering
		 */
		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const // y = Ax;
		{
			return applySpecialization(y,x,getType(_field));
		}

		/** applyTranspose.
		 *
		 * Uses one of the three
		 * private utility functions, in the manner described above.  Worthy of
		 * note is the fact that applyTranspose works by passing the column
		 * positions to the _apply functions as if they were rows, and row positions
		 * as if they were columns, as if the matrix had been transposed.
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const // y = ATx
		{
			return applyTransposeSpecialization(y,x,getType(_field));
		}

		size_t rowdim() const
		{
			return _rows;
		}

		size_t coldim() const
		{
			return _cols;
		}

		template<typename _Tp1>
		struct rebind {
			typedef ZeroOne<_Tp1> other;
			void operator() (other & Ap,
					 const Self_t& A,
					 const _Tp1& F)
			{
				// ZeroOne does not store any field element
			}
		};

		template<typename _Tp1>
		ZeroOne(const ZeroOne<_Tp1>& Z, const Field& F) :
			_field(F),
			_rows(Z.rowdim()), _cols(Z.coldim()), _nnz(Z.nnz()),
			_rowP(new Index[Z.nnz()]), _colP(new Index[Z.nnz()]),
			_rowSort(Z.isRowSorted()), _colSort(Z.isColSorted()),
			dynamic(true)
		{

			Index * rowit = _rowP;
			Index * colit = _colP;

			for(typename ZeroOne<_Tp1>::IndexIterator it = Z.indexBegin();
			    it != Z.indexEnd(); ++it,++rowit,++colit) {
				*rowit = (*it).first;
				*colit = (*it).second;
			}
		}

		/** Iterator class.
		 * Iterates straight through the values of the matrix
		 */
		class Iterator;

		Iterator Begin();
		Iterator End();
		const Iterator Begin() const;
		const Iterator End() const;

		/** IndexIterator.
		 * Iterates through the i and j of the current element
		 * and when accessed returns an STL pair containing the coordinates
		 */
		class IndexIterator;
		IndexIterator indexBegin();
		const IndexIterator indexBegin() const;
		IndexIterator indexEnd();
		const IndexIterator indexEnd() const;

		/** Read the matrix from a stream in the JGD's SMS format.
		 *  @param is Input stream from which to read the matrix
		 *  @return Reference to input stream
		 */
		std::istream &read (std::istream &is)
		{
			size_t i, j, k, m, n;

			char buf[80];
			buf[0]=0;
			is.getline (buf, 80);
			std::istringstream str (buf);
			str >> m >> n >> k;
			_rows = m;
			_cols = n;
			std::vector<size_t> rowP, colP;
			size_t x;
			while (is >> i >> j >> x) {
				if (i == 0 || i == (size_t) -1) break;
				if (x == 1UL) {
					rowP.push_back(i-1);
					colP.push_back(j-1);
				}
			}
			_nnz = rowP.size();
			_rowP = new size_t[_nnz];
			_colP = new size_t[_nnz];
			copy(rowP.begin(), rowP.end(), _rowP);
			copy(colP.begin(), colP.end(), _colP);
			return is;
		}

		std::ostream& write(std::ostream& out =std::cout)
		{
			size_t* i=_rowP;
			size_t* j=_colP;
			std::cout<<"Row dim: "<<rowdim()
			<<" Col dim: "<<coldim()
			<<" Total nnz: "<<nnz()<<"\n";
			for(;i<_rowP+nnz();++i,++j)
				std::cout<<*i<<" "<<*j<<"\n";
			return out;
		}

		const Field& field() const
		{
			return _field;
		}

		bool isRowSorted() const
		{
			return _rowSort;
		}
		bool isColSorted() const
		{
			return _colSort;
		}

		size_t nnz() const
		{
			return _nnz;
		};


	protected:


		Field _field; //!< @internal The field used by this class

		/*! @internal A temporary element used for initalization for the Begin() and
		 * End() methods of the ZeroOne class.  Is used to initalize a 1
		 * so that the Iterator returned stores a 1
		 */
		Element _tmp;

		Index _rows ;          //!<@internal number of rows of the Matrix
		Index _cols ;          //!<@internal number of columns
		Index _nnz;            //!<@internal Number of  Non-Zero elements in the Matrix.  It also happens to be the length of  the three NAGSparse arrays.
		mutable Index* _rowP ; //!<@internal pointer to an array of row indexes.
		mutable Index* _colP;  //!<@internal pointer to an array of column indexes. (\c _rowP and \c _colP are the other arrays of a  NAGSparse format Matrix.)
		mutable bool _rowSort ;
		mutable bool _colSort; //!<@internal status flags for sorting state
		bool dynamic;          // NO DOC

		/*! Tells the number of nonzero entries.
		 * Non blackbox function.
		 */
		void rowSort() const;
		void colSort() const;

		void _qsort(size_t start, size_t endp1, int &mode) const;   //!< @internal QuickSort function for when there is no sorting
		size_t _part( size_t start, size_t endp1, int &mode) const; //!< @internal Partition for quicksort

	private:

		class FieldType {};
		class NormField : public FieldType {};
		class Mod32Field : public FieldType {};

		template<class F>
		NormField getType(const F &  f) const
		{
			return NormField();
		}

		Mod32Field getType(const Modular<uint32_t> &) const
		{
			return Mod32Field();
		}

		template<class OutVector, class InVector>
		OutVector& applySpecialization(OutVector &, const InVector &,const NormField& ) const;
		template<class OutVector, class InVector>
		OutVector& applySpecialization(OutVector &, const InVector &, const Mod32Field& )const;
		template<class OutVector, class InVector>
		OutVector& applyTransposeSpecialization(OutVector &, const InVector &,const NormField& ) const;
		template<class OutVector, class InVector>
		OutVector& applyTransposeSpecialization(OutVector &, const InVector &, const Mod32Field& )const;

	}; //ZeroOne

} //LinBox

#include "linbox/blackbox/zero-one.inl"

#endif // __LINBOX_zero_one_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

