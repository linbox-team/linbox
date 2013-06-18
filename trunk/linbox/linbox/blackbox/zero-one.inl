
/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong Wan, -bds
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

/*! @file blackbox/zero-one.inl
 * @ingroup blackbox
 * @brief NO DOC
 */

#ifndef __LINBOX_bb_zero_one_INL
#define __LINBOX_bb_zero_one_INL

namespace LinBox
{
	/*! Raw iterator.
	 * @ingroup iterators
	 */
	template<class Field>
	class ZeroOne<Field>::Iterator {
	public:
		typedef Element value_type;

		Iterator(size_t pos, Element elem) :
			_pos(pos), _elem(elem) {}

		Iterator(const Iterator &In) :
			_pos(In._pos), _elem(In._elem) {}

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

		const value_type operator*() const { return _elem; }

	private:
		value_type _elem;
		size_t _pos;
	};

	/* STL standard Begin and End functions.  Used to get
	 * the beginning and end of the data.  So that Iterator
	 * can be used in algorithms like a normal STL iterator.
	 */
	template<class Field> typename
	ZeroOne<Field>::Iterator ZeroOne<Field>::Begin()
	{
	       	return Iterator( 0, field().init(_tmp, 1) );
	}

	template<class Field> typename
	ZeroOne<Field>::Iterator ZeroOne<Field>::End()
	{
	       	return Iterator( _nnz, field().init(_tmp, 1) );
	}

	template<class Field>
	const typename ZeroOne<Field>::Iterator ZeroOne<Field>::Begin() const
	{
	       	return Iterator(0, field().init(_tmp, 1) );
	}

	template<class Field>
	const typename ZeroOne<Field>::Iterator ZeroOne<Field>::End() const
	{
	       	return Iterator(_nnz, field().init(_tmp, 1) );
	}

	/*! IndexIterator.
	 * @ingroup iterators
	 * Iterates through the i and j of the current element
	 * and when accessed returns an STL pair containing the coordinates
	 */
	template<class Field>
	class ZeroOne<Field>::IndexIterator {
	public:
		typedef std::pair<size_t, size_t> value_type;

		IndexIterator() {}

		IndexIterator(size_t* row, size_t* col):
			_row(row), _col(col) {}

		IndexIterator(const IndexIterator &In):
			_row(In._row), _col(In._col)
		{}

		const IndexIterator &operator=(const IndexIterator &rhs)
		{
			_row = rhs._row;
			_col = rhs._col;
			return *this;
		}

		bool operator==(const IndexIterator &rhs)
		{
			return _row == rhs._row && _col == rhs._col;
		}

		bool operator!=(const IndexIterator &rhs)
		{
			return _row != rhs._row || _col != rhs._col;
		}

		const IndexIterator& operator++()
		{
			++_row; ++_col;
			return *this;
		}

		const IndexIterator operator++(int)
		{
			IndexIterator tmp = *this;
			++_row; ++_col;
			return tmp;
		}

		value_type operator*()
		{
			return std::pair<size_t,size_t>(*_row, *_col);
		}

		const value_type operator*() const
		{
			return std::pair<size_t,size_t>(*_row, *_col);
		}
	private:
		size_t* _row, *_col;
	};

	template<class Field> typename
	ZeroOne<Field>::IndexIterator ZeroOne<Field>::indexBegin()
	{
		return IndexIterator(_rowP, _colP);
	}

	template<class Field>
	const typename ZeroOne<Field>::IndexIterator ZeroOne<Field>::indexBegin() const
	{
		return IndexIterator(_rowP, _colP);
	}

	template<class Field> typename
	ZeroOne<Field>::IndexIterator ZeroOne<Field>::indexEnd()
	{
		return IndexIterator(_rowP + _nnz, _colP + _nnz);
	}

	template<class Field>
	const typename ZeroOne<Field>::IndexIterator ZeroOne<Field>::indexEnd() const
	{
		return IndexIterator(_rowP + _nnz, _colP + _nnz);
	}

	template<class Field>
	ZeroOne<Field>::ZeroOne(const Field& F) :
	       	_field(&F)
       	{
		srand((unsigned int) time(NULL) );
		dynamic = false;
	}

	template<class Field>
	ZeroOne<Field>::ZeroOne(Field& F, Index* rowP, Index* colP,
				Index rows, Index cols, Index NNz, bool RowSort, bool ColSort):
		_field(&F), _rows(rows), _cols(cols), _nnz(NNz), _rowP(rowP), _colP(colP), _rowSort(RowSort), _colSort(ColSort) , dynamic(false)
	{
	       	srand((unsigned)time(NULL));
	}

	template<class Field>
	ZeroOne<Field>::~ZeroOne()
	{
		if(dynamic) {
			delete [] _rowP;
			delete [] _colP;
		}
	}

	template<class Field>
	void ZeroOne<Field>::rowSort() const
	{
		if( _rowSort) return;  // Already sorted, we're done
		else {
			int mode = 0;
			_qsort( (size_t) 0, _nnz, mode);
		}
		_rowSort = true;
		return;
	}

	template<class Field>
	void ZeroOne<Field>::colSort() const
	{
		if( _colSort) return; // Already sorted, good to go
		else {
			int mode = 1;
			_qsort( (size_t) 0, _nnz, mode);
		}
		_colSort = true; _rowSort = false;
		return;
	}

	template<class Field>
	void ZeroOne<Field>::_qsort(size_t p, size_t e, int &mode) const
	{
		if( (e - p) <= 1) ;
		else
		{
			int i;
			i = 1 + (int)_part(p, e, mode);
			_qsort(p, (size_t)i, mode);
			_qsort((size_t)i, e, mode);
		}
	}

	template<class Field>
	size_t ZeroOne<Field>::_part(size_t p, size_t e, int &mode) const
	{
		size_t rtemp, ctemp, rowval, colval;
		int i = int(p +(size_t) rand() % (e - p)), j =(int) e;
		rtemp = _rowP[p];
		ctemp = _colP[p];
		_rowP[p] = _rowP[i];
		_colP[p] = _colP[i];
		_rowP[i] = rtemp;
		_colP[i] = ctemp;
		rowval = _rowP[p];
		colval = _colP[p];
		i = (int)p - 1;

		if(mode == 0)
		{ // Row mode, go by row order, then column
			while(true)
			{
				do j--; while( _rowP[j] > rowval || ( _rowP[j] == rowval && _colP[j] > colval ));
				do i++; while( _rowP[i] < rowval || ( _rowP[i] == rowval && _colP[i] < colval ));
				if( i < j)
				{
					rtemp = _rowP[j];
					ctemp = _colP[j];
					_rowP[j] = _rowP[i];
					_colP[j] = _colP[i];
					_rowP[i] = rtemp;
					_colP[i] = ctemp;
				}
				else return (size_t)j;
			}
		}
		else
		{ // Col mode, go by col order, then row
			while(true)
			{
				do j--; while( _colP[j] > colval || ( _colP[j] == colval && _rowP[j] > rowval ));
				do i++; while( _colP[i] < colval || ( _colP[i] == colval && _rowP[i] < rowval ));
				if( i < j) {
					rtemp = _rowP[j];
					ctemp = _colP[j];
					_rowP[j] = _rowP[i];
					_colP[j] = _colP[i];
					_rowP[i] = rtemp;
					_colP[i] = ctemp;
				}
				else return (size_t)j;
			}
		}
	}

	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::applySpecialization(OutVector & y, const InVector & x, const NormField& n) const
	{
		//std::cout<<"Call general case\n";
		linbox_check((y.size()==rowdim())&&(x.size()==coldim()));
		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		Index* ip, *jp;

		// 0 out y.  Note, this implementation assumes a dense vector.
		for(yp = y.begin(); yp != y.end(); ++yp)
			field().init(*yp , 0);

		rowSort();

		yp=y.begin();
		xp=x.begin();
		ip=_rowP;
		jp=_colP;
		size_t rowI =0;

		for(; ip <_rowP+nnz(); ++ip,++jp)
		{
			if( *ip == rowI)
				field().addin(*yp,*(xp +(ptrdiff_t) *jp));
			else
			{
				if((*ip-rowI)==1)
					++yp;
				else
					yp=y.begin()+(ptrdiff_t)*ip;

				rowI=*ip;
				field().addin(*yp,*(xp +(ptrdiff_t) *jp));
			}
		}
		return y;
	}


	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::applySpecialization(OutVector & y, const InVector & x, const Mod32Field& m) const
	{
		//std::cout<<"Called specialization\n";
		linbox_check((y.size()==rowdim())&&(x.size()==coldim()));

		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		typedef typename OutVector::value_type val_t ;
		Index* ip, *jp;

		for(yp = y.begin(); yp != y.end(); ++yp)
			field().init(*yp , 0);

		rowSort();

		yp=y.begin();
		xp=x.begin();
		ip=_rowP;
		jp=_colP;
		size_t rowI =0;
		integer _prime;

		field().characteristic(_prime);

		uint32_t prime = static_cast<uint32_t>(_prime);

		uint64_t accum =0;

		for(; ip <_rowP+nnz(); ++ip,++jp)
		{
			if( *ip == rowI)
				accum=accum+*(xp +(ptrdiff_t) *jp);
			else
			{
				*yp= (val_t)(accum % prime);
				if((*ip-rowI)==1)
					++yp;
				else
					yp=y.begin()+(ptrdiff_t)*ip;

				rowI=*ip;
				accum=*(xp+(ptrdiff_t)*jp);
			}
		}
		if(rowI)
			*yp= val_t(accum % prime);

		return y;
	}


	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::applyTransposeSpecialization(OutVector & y, const InVector & x, const NormField& n) const
	{
		//std::cout<<"Call general case\n";
		linbox_check((y.size()==coldim())&&(x.size()==rowdim()));
		typename OutVector::iterator yp;
		typename InVector::const_iterator xp;
		Index* ip, *jp;

		// 0 out y.  Note, this implementation assumes a dense vector.
		for(yp = y.begin(); yp != y.end(); ++yp)
			field().init(*yp , 0);

		rowSort();

		yp=y.begin();
		xp=x.begin();
		ip=_rowP;
		jp=_colP;
		size_t rowI =0;

		for(; ip <_rowP+nnz(); ++ip,++jp)
		{
			if( *ip == rowI)
				field().addin(*(yp+(ptrdiff_t)*jp),*xp);
			else
			{
				if((*ip-rowI)==1)
					++xp;
				else
					xp=x.begin()+(ptrdiff_t)*ip;

				rowI=*ip;
				field().addin(*(yp+(ptrdiff_t)*jp),*xp);
			}
		}

		return y;
	}


	template<class Field>
	template<class OutVector, class InVector>
	OutVector & ZeroOne<Field>::applyTransposeSpecialization(OutVector & y, const InVector & x, const Mod32Field& m) const
	{
		//std::cout<<"Called specialization\n";
		linbox_check((y.size()==coldim())&&(x.size()==rowdim()));

		std::vector<uint64_t> y_c (y.size(),0);

		typename OutVector::iterator         yp;
		typename InVector::const_iterator    xp;
		typedef typename OutVector::value_type   val_t ;
		Index* ip, *jp;

		rowSort();

		xp=x.begin();
		ip=_rowP;
		jp=_colP;
		size_t rowI =0;
		std::vector<uint64_t>::iterator y_cp;
		y_cp=y_c.begin();

		for(; ip <_rowP+nnz(); ++ip,++jp)
		{
			if( *ip == rowI) {
				*(y_cp+(ptrdiff_t)*jp) += *xp;
			}
			else {
				if((*ip-rowI)==1)
					++xp;
				else
					xp=x.begin()+(ptrdiff_t)*ip;

				rowI=*ip;
				*(y_cp+(ptrdiff_t)*jp) += *xp;
			}
		}

		integer _prime;
		field().characteristic(_prime);
		uint32_t prime = static_cast<uint32_t>(_prime);

		yp=y.begin();
		y_cp=y_c.begin();
		for(;yp!=y.end();++yp,++y_cp)
			*yp =(val_t)( (*y_cp) % prime );

		return y;
	}

}//End of LinBox

#endif // __LINBOX_bb_zero_one_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

