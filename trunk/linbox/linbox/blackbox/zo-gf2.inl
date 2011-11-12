/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/*! @file blackbox/zo-gf2.inl
 * @ingroup blackbox
 * @brief NO DOC
 */

#ifndef __LINBOX_zo_gf2_INL
#define __LINBOX_zo_gf2_INL

#include <givaro/givintfactor.h>

namespace LinBox
{
	// Dot product structure enabling std::transform call
	template<class Blackbox, class InVector>
	struct dotp {
		const typename Blackbox::Field& _field;
		const InVector& _x;
		dotp(const typename Blackbox::Field& F, const InVector& x) :
			_field(F), _x(x)
		{}

		bool operator()(const typename Blackbox::Row_t& row) const
		{
			bool tmp(false);
			for(typename Blackbox::Row_t::const_iterator loc = row.begin(); loc != row.end(); ++loc) {
				_field.addin(tmp,_x[*loc]);
			}
			return tmp;
		}
	};

#include <algorithm>
	template<class OutVector, class InVector>
	inline OutVector & ZeroOne<GF2>::apply(OutVector & y, const InVector & x) const
	{
		dotp<Self_t,InVector> mydp(this->_field, x);
		std::transform(this->begin(), this->end(), y.begin(), mydp );
		return y;
	}

#if 0
	template<class OutVector, class InVector>
	inline OutVector & ZeroOne<GF2>::apply(OutVector & y, const InVector & x) const
	{
		typename OutVector::iterator yit = y.begin();
		Self_t::const_iterator row = this->begin();
		for( ; row != this->end(); ++yit, ++row) {
			bool tmp(false);
			for(Row_t::const_iterator loc = row->begin();loc != row->end(); ++loc)
				_field.addin(tmp,x[*loc]);
			*yit = tmp;
		}
		return y;
	}
#endif

	template<class OutVector, class InVector>
	inline OutVector & ZeroOne<GF2>::applyTranspose(OutVector & y, const InVector & x) const
	{
		std::fill(y.begin(),y.end(),false);
		typename InVector::const_iterator xit = x.begin();
		Self_t::const_iterator row = this->begin();
		for( ; row != this->end(); ++row, ++xit) {
			for(typename Self_t::Row_t::const_iterator loc = row->begin(); loc != row->end(); ++loc) {
				_field.addin(y[*loc],*xit);
			}
		}
		return y;
	}


	inline void ZeroOne<GF2>::setEntry(size_t i, size_t j, const Element& v) {
		Row_t& rowi = this->operator[](i);
		Row_t::iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
		if (! _field.isZero(v) ) {
			if ( (there == rowi.end() ) || (*there != j) ) {
				rowi.insert(there, j);
				++_nnz;
			}
		}
		else {
			if ( (there != rowi.end() ) && (*there == j) ) {
				rowi.erase(there);
				--_nnz;
			}
		}
	}

	inline ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(Element& r, size_t i, size_t j) const
	{
		const Row_t& rowi = this->operator[](i);
		Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
		if (there != rowi.end() )
			return r=*there;
		else
			return r=_field.zero;
	}

	inline const ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(size_t i, size_t j) const
	{
		const Row_t& rowi = this->operator[](i);
		Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
		if (there != rowi.end() )
			return reinterpret_cast<const ZeroOne<GF2>::Element&>(*there);
		else
			return _field.zero;
	}

	inline std::istream &ZeroOne<GF2>::read (std::istream &is) {
		// Reads a long int and take it mod 2 afterwards (v&1)
		UnparametricField<long> Ints;
		MatrixStream<UnparametricField<long> > S(Ints, is);
		S.getDimensions( _rowdim, _coldim );
		this->resize(_rowdim);
		Index r, c;
		long v;
		_nnz = 0;
		while( S.nextTriple(r, c, v) ) {
			if (v&1) {
				this->operator[](r).push_back(c);
				++_nnz;
			}
		}
		for(Father_t::iterator i=this->begin(); i!=this->end(); ++i)
			std::sort(i->begin(),i->end());
		return is;
	}

	inline std::ostream& ZeroOne<GF2>::write (std::ostream& out, FileFormatTag format) const
	{
		if (format == FORMAT_GUILLAUME) {
			out << _rowdim << ' ' << _coldim << " M\n";
			for(size_t i=0; i<_rowdim; ++i) {
				const Row_t& rowi = this->operator[](i);
				for(Row_t::const_iterator it=rowi.begin(); it != rowi.end(); ++it)
					out << (i+1) << ' ' << (*it+1) << " 1\n";
			}
			return out << "0 0 0" << std::endl;
		}
		else if (format == FORMAT_MAPLE) {
			out << '[';
			bool firstrow=true;
			for (const_iterator i = begin (); i != end (); ++i) {
				if (firstrow) {
					out << '[';
					firstrow =false;
				}
				else
					out << ", [";

				Row_t::const_iterator j = i->begin ();
				for (long j_idx = 0; j_idx < static_cast<long>(_coldim); j_idx++) {
					if (j == i->end () || j_idx != static_cast<long>(*j) )
						out << '0';
					else {
						out << '1';
						++j;
					}
					if (j_idx < (static_cast<long>(_coldim)-1) )
						out << ',';
				}

				out << ']';
			}
			return out << ']';
		}
		else
			return out << "ZeroOne over GF(2), format other than SMS or Maple not implemented" << std::endl;
	}

	/*! Raw iterator.
	 * @ingroup iterators
	 */
	class ZeroOne<GF2>::Iterator {
	public:
		typedef Element value_type;

		Iterator(size_t pos, Element elem) :
			_elem(elem),_pos(pos)
		{}

		Iterator(const Iterator &In) :
			_elem(In._elem),_pos(In._pos)
		{}

		const Iterator& operator=(const Iterator& rhs)
		{
			_pos = rhs._pos;
			_elem = rhs._elem;
			return *this;
		}


		bool operator==(const Iterator &rhs)
		{
			return ( _pos == rhs._pos && _elem == rhs._elem);
		}

		bool operator!=(const Iterator &rhs)
		{
			return ( _pos != rhs._pos || _elem != rhs._elem );
		}

		Iterator & operator++()
		{
			++_pos;
			return *this;
		}

		Iterator operator++(int)
		{
			Iterator tmp = *this;
			_pos++;
			return tmp;
		}

		value_type operator*() { return _elem; }

		value_type operator*() const
		{
		       	return _elem;
		}

	private:
		value_type _elem;
		size_t _pos;
	};

	/* STL standard Begin and End functions.  Used to get
	 * the beginning and end of the data.  So that Iterator
	 * can be used in algorithms like a normal STL iterator.
	 */
	inline ZeroOne<GF2>::Iterator ZeroOne<GF2>::Begin()
	{
		return Iterator( 0, _field.one );
	}

	inline ZeroOne<GF2>::Iterator ZeroOne<GF2>::End()
	{
		return Iterator( _nnz, _field.one );
	}

	inline const ZeroOne<GF2>::Iterator ZeroOne<GF2>::Begin() const
	{
		return Iterator(0, _field.one );
	}

	inline const ZeroOne<GF2>::Iterator ZeroOne<GF2>::End() const
	{
		return Iterator(_nnz, _field.one );
	}

	/*! IndexIterator.
	 * @ingroup iterators
	 * Iterates through the i and j of the current element
	 * and when accessed returns an STL pair containing the coordinates
	 */
	class ZeroOne<GF2>::IndexIterator {
	public:
		typedef std::pair<size_t, size_t> value_type;

		IndexIterator() {}

		IndexIterator(size_t rowidx,
				 LightContainer<LightContainer<size_t> >::const_iterator rowbeg,
				 LightContainer<LightContainer<size_t> >::const_iterator rowend,
				 size_t colidx,
				 LightContainer<size_t>::const_iterator colbeg) :
			_rowbeg( LightContainer<LightContainer<size_t> >::const_iterator(rowbeg) ),
			_rowend( LightContainer<LightContainer<size_t> >::const_iterator(rowend) ),
			_colbeg( LightContainer<size_t>::const_iterator(colbeg) ),
			_row(rowidx),
			_col(colidx)
		{

			if( _rowbeg == _rowend ) return;

			while ( _colbeg == _rowbeg->end() ) {

				if (++_rowbeg == _rowend) return;

				_colbeg = _rowbeg->begin();

			}

		}

		IndexIterator(const IndexIterator &In) :
			_rowbeg(In._rowbeg), _rowend(In._rowend), _colbeg(In._colbeg), _row(In._row), _col(In._col)
		{}

		const IndexIterator &operator=(const IndexIterator &rhs)
		{
			_rowbeg = rhs._rowbeg;
			_rowend = rhs._rowend;
			_colbeg = rhs._colbeg;
			_row = rhs._row;
			_col = rhs._col;
			return *this;
		}

		bool operator==(const IndexIterator &rhs)
		{
			return _rowbeg == rhs._rowbeg && _colbeg == rhs._colbeg;
		}

		bool operator!=(const IndexIterator &rhs)
		{
			return _rowbeg != rhs._rowbeg || _colbeg != rhs._colbeg;
		}

		const IndexIterator& operator++() {



			++_colbeg;
			while(_colbeg == _rowbeg->end()) {
				if (++_rowbeg == _rowend) return *this;
				++_row;
				_colbeg = _rowbeg->begin();
			}
			_col = *_colbeg;


			return *this;
		}

		const IndexIterator operator++(int)
		{
			IndexIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		value_type operator*()
		{
			return std::pair<size_t,size_t>(_row, _col);
		}

		const value_type operator*() const
		{
			return std::pair<size_t,size_t>(_row, _col);
		}
	private:
		LightContainer<LightContainer<size_t> >::const_iterator _rowbeg, _rowend;
		LightContainer<size_t>::const_iterator _colbeg;
		size_t _row, _col;
	};

	inline ZeroOne<GF2>::IndexIterator ZeroOne<GF2>::indexBegin()
	{
		return IndexIterator(0, this->begin(), this->end(), 0, this->front().begin() );
	}

	inline const ZeroOne<GF2>::IndexIterator ZeroOne<GF2>::indexBegin() const
	{
		return IndexIterator(0, this->begin(), this->end(), 0, this->front().begin() );
	}

	inline ZeroOne<GF2>::IndexIterator ZeroOne<GF2>::indexEnd()
	{
		return IndexIterator(_rowdim, this->end(), this->end(), this->back().size(),this->back().end() );
	}

	inline const ZeroOne<GF2>::IndexIterator ZeroOne<GF2>::indexEnd() const
	{
		return IndexIterator(_rowdim, this->end(), this->end(), this->back().size(),this->back().end() );
	}





} // end of namespace LinBox


// Specialization of getentry
#include "linbox/solutions/getentry.h"
namespace LinBox
{
	template<> struct GetEntryCategory<ZeroOne<GF2> > { typedef GetEntryTags::Local Tag; };
	} // end of namespace LinBox

#endif //__LINBOX_zo_gf2_INL
