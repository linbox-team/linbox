/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@lirmm.fr
 *               Clément Pernet clement.pernet@imag.fr
 *               Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file matrix/densematrix/new-blas-submatrix.inl
 * @ingroup densematrix
 */


namespace LinBox {

    //////////////////
    // CONSTRUCTORS //
    //////////////////

    /*  constructors */
    
    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (typename BlasSubmatrix<_Matrix>::matrixType &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
		_stride(M.getStride()),
        _field(M.field()),
		_AD(M.field())
	{
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (BlasSubmatrix<_Matrix>  &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim) :
        _ptr(M.getPointer()+rowbeg*M.getStride()+colbeg),
        _row (Rowdim),
        _col(Coldim),
		_stride(M.getStride()),
        _field(M.field()),
		_AD(M.field())
	{
        linbox_check ( rowbeg  <= M.rowdim() ); // allow for NULL matrix
		linbox_check ( colbeg  <= M.coldim() );
        linbox_check ( rowbeg+Rowdim <= M.rowdim() );
        linbox_check ( colbeg+Coldim <= M.coldim() );
    }

    template < class _Matrix >
    BlasSubmatrix<_Matrix>::BlasSubmatrix (typename BlasSubmatrix<_Matrix>::matrixType &M):        
        _ptr(M.getPointer()),
        _row (M.rowdim()),
        _col (M.coldim()),
		_stride(M.getStride()),
        _field(M.field()),
		_AD(M.field())
    {

    }

  

    //////////////////
    //   ELEMENTS   //
    //////////////////

    template < class _Matrix >
    const typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::setEntry (size_t i, size_t j, const Element &a_ij) {
        return field().assign(a_ij, _ptr[i*_stride+j]);
    }

    template < class _Matrix >
    typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::refEntry (size_t i, size_t j) {
        return _ptr[i*_stride+j];
    }

    template < class _Matrix >
    const typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::getEntry (size_t i, size_t j) const {
        return _ptr[i*_stride+j];
    }

    template < class _Matrix >
    typename BlasSubmatrix<_Matrix>::Element& BlasSubmatrix<_Matrix>::getEntry (Element &x, size_t i, size_t j) const {
        return field().assign(x,_ptr[i*_stride+j]);
    }




    ///////////////////
    //      I/O      //
    ///////////////////

    template < class _Matrix >
    std::istream& BlasSubmatrix<_Matrix>::read (std::istream &file)
    {
        // must be improved maybe to allow any matrix format
		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;        		
        linbox_check (m == _row);
        linbox_check (m == _col);
		
		if ((c != 'M') && (c != 'm')) {
            for (p = Begin (); p != End (); ++p) {
                field().read (file, *p);
			}
		}
		else { // sparse file format
			int i=0, j=0;
			while (true)
                {
                    file >> i >> j;
                    if (i+j <= 0) break;
                    field().read (file, _ptr+ (i-1)*_stride+(j-1));
                }
		}

		return file;
    }

    template < class _Matrix >
    std::ostream& BlasSubmatrix<_Matrix>::write (std::ostream &os,LINBOX_enum (Tag::FileFormat) f) const {
        
		ConstRowIterator p;
		switch(f) {
        case (Tag::FileFormat::MatrixMarket ) : /* Matrix Market */
            {
                writeMMArray(os, *this, "BlasSubmatrix");
            }
            break;
        case (Tag::FileFormat::Plain) : /*  raw output */
            {
                integer c;
                int wid;

                field().cardinality (c);

                if (c >0)
                    wid = (int) ceil (log ((double) c) / M_LN10);
                else {
                    wid=1000;
                }
                    
                for (p = rowBegin (); p != rowEnd ();++p) {
                    typename ConstRow::const_iterator pe;

                    os << "  [ ";

                    for (pe = p->begin (); pe != p->end (); ++pe) {
                        os.width (wid);
                        field().write (os, *pe);
                        os << " ";
                    }
                    os << "]" << std::endl;
                }
            }
            break;
        case (Tag::FileFormat::Maple) : /*  maple format */
            {

                os << "Matrix( " << rowdim() << ',' << coldim() << ",\n[" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;
                    if (p!=rowBegin()) os << ' ';
                    os << "[ ";
                        
                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os, *pe);
                        ++pe ;
                        if (pe != p->end())
                            os << ", ";
                    }

                    os << "]" ;
                    ++p ;
                    if (p != rowEnd() )
                        os << ',' << std::endl;;

                }
                os << "])" ;
            }
            break;
        case (Tag::FileFormat::HTML) : /*  HTML format */
            {

                os << "<table border=\"1\">" ;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    os << "<tr>";

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os<< "<td>", *pe)<<"</td>";
                        ++pe ;
                    }

                    os << "</tr>" << std::endl;
                    ++p ;

                }
                os << "</table>" ;
            }
            break;
        case (Tag::FileFormat::LaTeX) : /*  LaTex format */
            {
                os << "\\begin{pmatrix} " << std::endl;
                for (p = rowBegin (); p != rowEnd (); ) {
                    typename ConstRow::const_iterator pe;

                    for (pe = p->begin (); pe != p->end (); ) {
                        field().write (os, *pe);
                        ++pe ;
                        if (pe != p->end())
                            os << "& ";
                    }

                    os << "\\\\" << std::endl;
                    ++p ;

                }
                os << "\\end{pmatrix}" ;
            }
            break;
        default : /*  this is an error */
            {
                throw LinBoxError("unknown format to write matrix in");
            }
		}

		return os;
	}

    


    ///////////////////
    //   ITERATORS   //
    ///////////////////



	/*! Raw Iterators.
	 * @ingroup iterators
	 *
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::Iterator {
	public:
		Iterator (){}

		/*! @internal
		 * @brief NO DOC
		 */
		Iterator (const Element_ptr& cur,
                  const size_t c_dim,
                  const size_t stride,
                  const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
		{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		Iterator& operator = (const Iterator& r)
		{
			_cur    = r._cur;
			_c_dim  = r._c_dim;
			_stride = r._stride;
			_c_idx  = r._c_idx;
			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		Iterator& operator ++()
		{
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
			}

			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		Iterator& operator++ (int)
		{
			return this->operator++ ();
		}


		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const Iterator& r) const
		{
			return (_cur != r._cur || _c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		typename _Matrix::Element& operator * ()
		{
			return *_cur;
		}

		//! @internal operator *.
		const typename _Matrix::Element& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<typename _Matrix::Field,typename _Matrix::Rep>::Iterator _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};

	/*! Raw Iterators (const version).
	 * @ingroup iterators
	 * The raw iterator is a method for accessing all entries in the matrix
	 * in some unspecified order. This can be used, e.g. to reduce all
	 * matrix entries modulo a prime before passing the matrix into an
	 * algorithm.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::ConstIterator {
	public:
		//! @internal Null constructor
		ConstIterator (){}


		/*! @internal
		 * @brief NO DOC
		 */
		ConstIterator (const typename BlasMatrix<typename _Matrix::Field, typename _Matrix::Rep>::ConstIterator& cur,
                       const size_t c_dim,
                       const size_t stride,
                       const size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride(stride), _c_idx (c_idx)
		{}

		/*! @internal
		 * @brief copy operator.
		 * @param r Iterator to copy.
		 */
		ConstIterator& operator = (const ConstIterator& r)
		{
			_cur = r._cur;
			_c_dim = r._c_dim;
			_stride = r._stride;
			_c_idx = r._c_idx;
			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator ++()
		{
			if (_c_idx < _c_dim - 1){
				++_cur; ++_c_idx;
			}
			else {
				linbox_check(_stride > _c_dim);
				_cur = _cur + (ptrdiff_t)(_stride - _c_dim + 1);
				_c_idx = 0;
			}

			return *this;
		}

		/*! @internal
		 * increment.
		 * ??
		 */
		ConstIterator& operator++ (int)
		{
			return this->operator++ ();
		}

		/*! @internal
		 * @brief  operator !=.
		 * @param r Iterator to test inequaltity from.
		 */
		bool operator != (const ConstIterator& r) const
		{
			return (_cur != r._cur) || (_c_dim != r._c_dim) || (_stride != r._stride) || (_c_idx != r._c_idx);
		}

		//! @internal operator *.
		const typename BlasSubmatrix<_Matrix>::Element& operator * () const
		{
			return *_cur;
		}

	protected:
		typename BlasMatrix<typename _Matrix::Field, typename _Matrix::Rep>::ConstIterator _cur;
		size_t _c_dim;
		size_t _stride;
		size_t _c_idx;
	};


	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::Iterator BlasSubmatrix<_Matrix>::Begin ()
	{
		return Iterator (_Mat.Begin () + (ptrdiff_t)( _off ),
                         _col, _stride, 0);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::Iterator BlasSubmatrix<_Matrix>::End ()
	{
		return Iterator (_Mat.Begin () +(ptrdiff_t) ( (_row) * _stride + _off ),
                         _col, _stride, 0);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstIterator BlasSubmatrix<_Matrix>::Begin () const
	{
		return ConstIterator (_Mat.Begin () +(ptrdiff_t) ( _off ),
                              _col, _stride, 0);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstIterator BlasSubmatrix<_Matrix>::End () const
	{
		return ConstIterator (_Mat.Begin () +(ptrdiff_t) ( (_row) * _stride + _off ),
                              _col, _stride, 0);
	}

	/*! Raw Indexed Iterator.
	 * @ingroup iterators
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's \c rowIndex() and \c colIndex() functions.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::IndexedIterator {
	public:
		IndexedIterator (){}

		IndexedIterator (const typename BlasMatrix<typename _Matrix::Field,typename _Matrix::Rep>::Iterator& cur,
                         size_t c_dim,
                         size_t stride,
                         size_t r_idx,
                         size_t c_idx) :
			_cur (cur), _c_dim (c_dim), _stride (stride), _r_idx(r_idx), _c_idx (c_idx)
		{}

		IndexedIterator& operator = (const IndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		IndexedIterator& operator++()
		{
			if (_c_idx < _c_dim - 1){
				++_c_idx;
				++_cur;
			}
			else {
				_cur = _cur + _stride - _c_dim + 1;
				_c_idx = 0;
				++_r_idx;
			}
			return *this;
		}

		IndexedIterator& operator--()
		{
			if (_c_idx > 0){
				--_c_idx;
				--_cur;
			}
			else {
				_cur = _cur - _stride + _c_dim -1;
				_c_idx = 0;
				--_r_idx;
			}
			return *this;
		}

		IndexedIterator operator++(int)
		{
			IndexedIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		IndexedIterator operator--(int)
		{
			IndexedIterator tmp = *this;
			this->operator--();
			return tmp;
		}

		bool operator != (const IndexedIterator& r) const
		{
			return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
		}

		const typename _Matrix::Field& operator*() const {return *_cur;}

		typename _Matrix::Element& operator*() {return *_cur;}

		size_t rowIndex () const { return _r_idx; }

		size_t colIndex () const { return _c_idx; }

		const typename _Matrix::Element& value () const {return *_cur;}

	protected:
		typename BlasMatrix<typename _Matrix::Field, typename _Matrix::Rep>::Iterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::IndexedIterator BlasSubmatrix<_Matrix>::IndexedBegin ()
	{
		return IndexedIterator (_Mat.Begin () +(ptrdiff_t) ( (_off) ),
                                _col , _stride, 0, 0);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::IndexedIterator BlasSubmatrix<_Matrix>::IndexedEnd ()
	{
		return IndexedIterator (_Mat.Begin () +(ptrdiff_t) ( (_row) * _stride + (_col+_off) ),
                                _col, _stride, _row-1, _col-1);
	}

	/*! Raw Indexed Iterator (const version).
	 * @ingroup iterators
	 *
	 * Like the raw iterator, the indexed iterator is a method for
	 * accessing all entries in the matrix in some unspecified order.
	 * At each position of the the indexed iterator, it also provides
	 * the row and column indices of the currently referenced entry.
	 * This is provided through it's \c rowIndex() and \c colIndex() functions.
	 */
	template < class _Matrix >
	class BlasSubmatrix<_Matrix>::ConstIndexedIterator {
	public:
		ConstIndexedIterator (){}

		ConstIndexedIterator (const typename BlasMatrix<typename _Matrix::Field, typename _Matrix::Rep>::ConstIterator& cur,
                              size_t c_dim,
                              size_t stride,
                              size_t r_idx,
                              size_t c_idx) :
			_cur (cur), _stride (stride), _c_dim (c_dim), _r_idx(r_idx), _c_idx (c_idx)
		{}

		ConstIndexedIterator& operator = (const IndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		ConstIndexedIterator& operator = (const ConstIndexedIterator& r)
		{
			_cur = r._cur;
			_stride = r._stride;
			_c_dim = r._c_dim;
			_r_idx = r._r_idx;
			_c_idx = r._c_idx;
			return *this;
		}

		ConstIndexedIterator& operator++()
		{
			if (_c_idx < _c_dim - 1){
				++_c_idx;
				++_cur;
			}
			else {
				_cur = _cur + _stride - _c_dim +1;
				_c_idx = 0;
				++_r_idx;
			}
			return *this;
		}

		IndexedIterator& operator--()
		{
			if (_c_idx > 0){
				--_c_idx;
				--_cur;
			}
			else {
				_cur = _cur - _stride + _c_dim -1;
				_c_idx = 0;
				--_r_idx;
			}
			return *this;
		}

		ConstIndexedIterator operator++(int)
		{
			ConstIndexedIterator tmp = *this;
			this->operator++();
			return tmp;
		}

		ConstIndexedIterator operator--(int)
		{
			ConstIndexedIterator tmp = *this;
			this->operator--();
			return tmp;
		}

		size_t rowIndex () const
		{
			return _r_idx;
		}

		size_t colIndex () const
		{
			return _c_idx;
		}

		bool operator != (const ConstIndexedIterator& r) const
		{
			return ((_c_idx != r._c_idx) || (_r_idx != r._r_idx) ||(_stride != r._stride) || (_c_dim != r._c_dim) );
		}

		const typename _Matrix::Element& operator*() const
		{
			return *_cur;
		}


		friend std::ostream& operator<<(std::ostream& out, const ConstIndexedIterator m)
		{
			return out /* << m._cur << ' ' */
                << m._stride << ' '
                << m._c_dim << ' '
                << m._r_idx << ' '
                << m._c_idx;
		}

		const typename _Matrix::Element & value() const
		{
			return this->operator*();

		}

	protected:
		typename BlasMatrix<typename _Matrix::Field, typename _Matrix::Rep>::ConstIterator _cur;
		size_t _stride;
		size_t _c_dim;
		size_t _r_idx;
		size_t _c_idx;
	};

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstIndexedIterator BlasSubmatrix<_Matrix>::IndexedBegin () const
	{
		return ConstIndexedIterator (_Mat.Begin () +(ptrdiff_t) ( _off ),
                                     _row, _stride, 0, 0);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstIndexedIterator BlasSubmatrix<_Matrix>::IndexedEnd () const
	{
		return ConstIndexedIterator (_Mat.Begin () +(ptrdiff_t) ( (_row) * _stride + (_off+_col) ),
                                     _col, _stride, _row-1, _col-1);
	}

	////////
	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::RowIterator BlasSubmatrix<_Matrix>::rowBegin ()
	{
		return RowIterator (_Mat.Begin () +(ptrdiff_t) ( _off  ),
                            _col, _stride);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::RowIterator BlasSubmatrix<_Matrix>::rowEnd ()
	{
		return RowIterator (_Mat.Begin () +(ptrdiff_t) ( (_row) * _stride + _off ),
                            _col, _stride);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstRowIterator BlasSubmatrix<_Matrix>::rowBegin () const
	{
		return ConstRowIterator (_Mat.Begin () + (ptrdiff_t)( _off ),
                                 _col, _stride);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstRowIterator BlasSubmatrix<_Matrix>::rowEnd () const
	{
		return ConstRowIterator (_Mat.Begin () + (ptrdiff_t)( (_row) * _stride + _off ),
                                 _col, _stride);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ColIterator BlasSubmatrix<_Matrix>::colBegin ()
	{
		return ColIterator (_Mat.Begin () + (ptrdiff_t)( _off ),
                            _stride, _row);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ColIterator BlasSubmatrix<_Matrix>::colEnd ()
	{
		return ColIterator (_Mat.Begin () + (ptrdiff_t)( (_col) + _off ),
                            _stride, _row);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstColIterator BlasSubmatrix<_Matrix>::colBegin () const
	{
		return ConstColIterator (_Mat.Begin () + (ptrdiff_t)( _off ),
                                 _stride, _row);
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstColIterator BlasSubmatrix<_Matrix>::colEnd () const
	{
		return ConstColIterator (_Mat.Begin () + (ptrdiff_t)( (_col) + _off ),
                                 _stride, _row);
	}

	/*  operators */
	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::Row BlasSubmatrix<_Matrix>::operator[] (size_t i)
	{
		return Row (_Mat.Begin () +(ptrdiff_t) (_r0+i) * _stride, _Mat.Begin () + (ptrdiff_t)((_r0+i) * _stride + _stride) );
	}

	template < class _Matrix >
	typename BlasSubmatrix<_Matrix>::ConstRow BlasSubmatrix<_Matrix>::operator[] (size_t i) const
	{
		return Row (_Mat.Begin () + (ptrdiff_t)(_r0+i) * _stride, _Mat.Begin () + (ptrdiff_t)((_r0+i) * _stride + _stride) );
	}

    
    
}



/* // Local Variables: */
/* // mode: C++ */
/* // tab-width: 4 */
/* // indent-tabs-mode: nil */
/* // c-basic-offset: 4 */
/* // End: */
/* // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s */

