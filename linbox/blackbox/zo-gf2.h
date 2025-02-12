/* linbox/blackbox/zo-gf2.h
 * Copyright (C) 2009,2010 The LinBox group
 *
 * Time-stamp: <10 May 23 18:18:43 Jean-Guillaume.Dumas@imag.fr>
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
 *
 */
#ifndef __LINBOX_zo_gf2_H
#define __LINBOX_zo_gf2_H

#include <algorithm>
#include "linbox/blackbox/zero-one.h"
#include "linbox/field/gf2.h"
#include <givaro/zring.h>
#include "linbox/util/matrix-stream.h"
#include "linbox/vector/stream.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/vector/light_container.h"

namespace LinBox
{

	/** \brief Time and space efficient representation of sparse matrices over GF2.
	 * Representation if a full row array containing vector of non-zero locations
	 * A 0-1 matrix is a matrix with all 0's and 1's as entries.
	 \ingroup blackbox
	 */
	template<>
	class ZeroOne<GF2> : public LightContainer< LightContainer< size_t > > {
	public:
		typedef LightContainer< LightContainer< size_t > > Father_t;
		typedef LightContainer< size_t > Row_t;
		typedef GF2::Element Element;
		typedef size_t Index;
		typedef ZeroOne<GF2> Self_t;
		typedef GF2 Field;

		const GF2 *_field;

		ZeroOne(const GF2& ) :
			_nnz(0)
		{}
		ZeroOne(const GF2& , const size_t m) :
			Father_t(m), _rowdim(m), _coldim(m),_nnz(0)
		{}
		ZeroOne(const GF2& , const size_t m, const size_t n) :
			Father_t(m), _rowdim(m), _coldim(n),_nnz(0)
		{}

		ZeroOne():
			_nnz(0)
		{}
		ZeroOne(const size_t m) :
			Father_t(m), _rowdim(m), _coldim(m),_nnz(0)
		{}
		ZeroOne(const size_t m, const size_t n) :
			Father_t(m), _rowdim(m), _coldim(n),_nnz(0)
		{}

		ZeroOne(const GF2& , VectorStream<Row_t>& stream) :
			Father_t(stream.m()), _rowdim(stream.m()), _coldim(stream.n()), _nnz(0)
		{
			for (Father_t::iterator row=begin(); row != end(); ++row) {
				stream >> *row;
				_nnz += row->size();
			}
		}

		ZeroOne(const Self_t& A) :
			Father_t(static_cast<const Father_t&>(A)), _rowdim(A._rowdim), _coldim(A._coldim), _nnz(A._nnz)
		{ }

		ZeroOne(const GF2& , size_t* rowP, size_t* colP,
			const size_t m, const size_t n, const size_t Nnz, const bool ,const bool) :
			Father_t(m), _rowdim(m), _coldim(n), _nnz(Nnz)
		{
			for(size_t k=0; k<Nnz; ++k)
				this->operator[](rowP[k]).push_back(colP[k]);
		}


		size_t rowdim() const { return _rowdim; }
		size_t coldim() const { return _coldim; }


		const Element& setEntry(size_t i, size_t j, const Element& v) ;
		const Element& getEntry(size_t i, size_t j) const ;
		Element& getEntry(Element&, size_t i, size_t j) const ;

		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const; // y = A x

		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const; // y = A^T x

		/** Read the matrix from a stream in ANY format
		 *  entries are read as "long int" and set to 1 if they are odd,
		 *  0 otherwise
		 *  @param is Input stream from which to read the matrix
		 *  @return Reference to input stream
		 */
		std::istream &read (std::istream &is) ;
		std::ostream& write (std::ostream& out, Tag::FileFormat format=Tag::FileFormat::SMS) const ;

		const Field& field() const { return *_field; }

		template<typename _Tp1>
		struct rebind;


		template<typename _Tp1>
		ZeroOne(const ZeroOne<_Tp1>& A, const GF2 F2);

		size_t nnz() const { return _nnz; }
		bool isRowSorted() const { return true; }
		bool isColSorted() const { return true; }

		/** Iterator class.  Iterates straight through the values of the matrix
		*/
		class Iterator;

		Iterator Begin();
		Iterator End();
		const Iterator Begin() const;
		const Iterator End() const;

		/** IndexedIterator - Iterates through the i and j of the current element
		 * and when accessed returns an STL pair containing the coordinates
		 */
		class IndexedIterator;
		IndexedIterator indexBegin();
		const IndexedIterator indexBegin() const;
		IndexedIterator indexEnd();
		const IndexedIterator indexEnd() const;

    protected:
        // Merge A with self
        // Warning: respective supports must be disjoint
        template<typename _Tp1>
        void augment(const ZeroOne<_Tp1>&);

	private:
		size_t _rowdim, _coldim, _nnz;
	};

}

#include "linbox/blackbox/zo-gf2.inl"

#endif //__LINBOX_zo_gf2_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
