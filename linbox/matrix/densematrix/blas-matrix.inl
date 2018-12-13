/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@lirmm.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
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

/*!@internal
 * @file matrix/densematrix/blas-matrix.inl
 * @ingroup densematrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
 */

#ifndef __LINBOX_densematrix_blas_matrix_INL
#define __LINBOX_densematrix_blas_matrix_INL


#include "linbox/field/hom.h"
#include "linbox/matrix/matrixdomain/apply-domain.h"
#include "fflas-ffpack/fflas/fflas.h"

///////////////////
//   PROTECTED   //
///////////////////

namespace LinBox
{

	template<class _Field, class _Storage>
    template <class constIterator>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const constIterator& v)
	{
		constIterator v_end = v+(_col*_row) ;
		Element_ptr iter_addr = _ptr ;
		for (; v != v_end ; ++v, ++iter_addr)
            {
                field().init(*iter_addr);
                field().assign(*iter_addr,*v);
            }
	}

	template<class _Field, class _Storage>
	template <class _Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const _Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::Container)
	{
        linbox_check( areFieldEqual(A.field(),field() ) );

		typename _Matrix::ConstIterator         iter_value = A.Begin();
		typename _Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

        // PG -> BUG if we use iter_value !=A.End() 
		for (;iter_index != A.IndexedEnd(); ++iter_value,++iter_index){
			int64_t i,j;
			i=(int64_t)iter_index.rowIndex()-(int64_t)i0;
			j=(int64_t)iter_index.colIndex()-(int64_t)j0;
			if ( (i>=0) && (j>=0) && (i< (int64_t)m) && (j < (int64_t)n))
				setEntry((size_t)i, (size_t)j, *iter_value);
		}
	}

    template<class _Field, class _Storage>
	template <class _Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const _Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::BlasContainer)
	{
        linbox_check( areFieldEqual(A.field(),field() ) );
        FFLAS::fassign(field(), m, n, A.getPointer(), A.getStride(), _ptr, getStride());
	}

	template<class _Field, class _Storage>
	template <class Matrix>
	void BlasMatrix< _Field, _Storage >::createBlasMatrix (const Matrix& A,
                                                           const size_t i0,const size_t j0,
                                                           const size_t m, const size_t n,
                                                           MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );

        BlasVector<Field> e(A.field(),A.coldim(), field().zero), tmp(A.field(),A.rowdim());
		ColIterator col_p;
        
		typename BlasMatrix< _Field, _Storage >::Col::iterator elt_p;
		typename BlasVector<Field>::iterator e_p, tmp_p;


		for (col_p = colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p)
            {

                field().assign(*e_p, field().one);

                A.apply (tmp, e);

                for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
                     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

                    field().assign(*elt_p, *tmp_p);

                field().assign(*e_p, field().zero);
            }
	}
            	
    ////////////
    // MEMORY //
    ////////////

    
	template < class _Field, class _Storage >
	void BlasMatrix< _Field, _Storage >::init(const size_t & m, const size_t & n)
	{
        _row = m; _col = n;
		_rep.resize(m*n, field().zero);
        //FFLAS::finit(field(), _row, _col, _ptr, _col); not yet in FFLAS
        for (size_t i=0;i<_row*_col;i++)
            field().init(_rep[i]);
		_ptr = _rep.data();
	}

    template < class _Field, class _Storage >
	void BlasMatrix< _Field, _Storage >::resize (const size_t & m, const size_t & n, const Element& val )
	{
#ifndef NDEBUG
		if (_col > 0 && _col != n)
			std::cerr << " ***Warning*** you are resizing a matrix, possibly loosing data. " << std::endl;
#endif
		_rep.resize (m * n, val);
		_ptr = _rep.data();
        if ( m*n > _row*_col){
            //FFLAS::finit(field(), m*n-_row*_col, _ptr+_row*_col); not yet in FFLAS
            for (size_t i=0;i<m*n-_row*_col;i++)
                field().init((_ptr+_row*_col)[i]);
        }
		_row = m;
		_col = n;
	}


    
    //////////////////
    // CONSTRUCTORS //
    //////////////////

    
	template < class _Field, class _Storage >
	BlasMatrix< _Field, _Storage >::BlasMatrix (const _Field &F) :
		_row(0),_col(0),_rep(0),_ptr(NULL),_field(F){}


	template < class _Field, class _Storage >
	BlasMatrix< _Field, _Storage >::BlasMatrix ( const _Field &F, const size_t & m, const size_t & n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(_rep.data()),_field(F)
    {
        // FFLAS::finit(field(), _row, _col, _ptr, getStride()); NO YET IN FFLAS
        for (size_t i=0;i<_row*_col;i++)
            field().init(_ptr[i]);
    }


	template < class _Field, class _Storage >
	template <class Matrix>
	BlasMatrix< _Field, _Storage >::BlasMatrix (const Matrix &A) : 
		_row(A.rowdim()),_col(A.coldim()),_rep(_row*_col),_ptr(_rep.data()),_field(A.field())
    {
        createBlasMatrix(A,0,0,_row,_col,typename MatrixContainerTrait<Matrix>::Type());
    }

    template < class _Field, class _Storage >
    template <class Matrix>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const Matrix& A,const size_t &i0, const size_t &j0,const size_t &m,  const size_t &n) :
        _row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0]),_field(A.field())
    {
                   
        createBlasMatrix(A, i0, j0, m, n,typename MatrixContainerTrait<Matrix>::Type());
    }

    template < class _Field, class _Storage >
    template<class constIterator>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const _Field &F,
                                                const size_t & m, const size_t & n,
                                                const constIterator& it) :
        _row(m), _col(n),_rep(_row*_col),_ptr(_rep.data()),
        _field(F)
    {
        createBlasMatrix(it);
    }

    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >::BlasMatrix(MatrixStream<_Field>& ms) :
        _row(0),_col(0),_rep(0),
        _field((ms.getField()))
    {
        if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
            throw ms.reportError(__FUNCTION__,__LINE__);
        _ptr = _rep.data();
    }

    template < class _Field, class _Storage >
    template <class StreamVector>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) :
        _row(stream.size ()), _col(stream.dim ()), _rep(_row*_col), _ptr(_rep.data()),
        _field (F)
    {
        VectorDomain<Field>    VD(F);
        StreamVector tmp(F);
        typename BlasMatrix<Field,_Storage>::RowIterator p;
        VectorWrapper::ensureDim (tmp, stream.dim ());
        for (p = BlasMatrix<Field,_Storage>::rowBegin (); p != BlasMatrix<Field,_Storage>::rowEnd (); ++p) {
            stream >> tmp;

            VD.copy (*p, tmp);
        }
    }

    template < class _Field, class _Storage >
    template<class OtherMatrix>
    BlasMatrix< _Field, _Storage >::BlasMatrix (const OtherMatrix &A, const _Field &F) :
        _row(A.rowdim()), _col(A.coldim()), _rep(_row*_col), _ptr(_rep.data()), _field(F)
    {
		typename OtherMatrix::template rebind<_Field>()(*this,A);        
    }

    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >::BlasMatrix (const  BlasMatrix< _Field, _Storage>  &A) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(_rep.data()),_field(A.field())
	{
        createBlasMatrix(A,0,0,_row,_col,MatrixContainerCategory::BlasContainer());
	}

    
    
    template < class _Field, class _Storage >
    BlasMatrix< _Field, _Storage >& BlasMatrix< _Field, _Storage >::operator= (const BlasMatrix< _Field, _Storage >& A)
    {
        if ( &A == this)
            return *this;

        _col = A.coldim();
        _row = A.rowdim();
        _rep = Storage(_row*_col);
        _ptr = _rep.data() ;
        createBlasMatrix(A,0,0,_row,_col,MatrixContainerCategory::BlasContainer());
        return *this;
    }
   
    template < class _Field, class _Storage >
    template < class Matrix>
    BlasMatrix< _Field, _Storage >& BlasMatrix< _Field, _Storage >::operator= (const Matrix& A)
    {
        _col = A.coldim();
        _row = A.rowdim();
        _rep = Storage(_row*_col);
        _ptr = _rep.data() ;
        createBlasMatrix(A,0,0,_row,_col,typename MatrixContainerTrait<Matrix>::Type());
        return *this;
    }


    //! @bug other Storage
    template < class _Field, class _Storage >
    template<typename _Tp1>
    struct BlasMatrix< _Field, _Storage >::rebind {
		typedef BlasMatrix<_Tp1,typename Vector<_Tp1>::Dense> other;

		void operator() (other & Ap, const Self_t& A) {
			// typedef Self_t::ConstIterator ConstSelfIterator ;
			typedef typename BlasMatrix< _Field, _Storage >::ConstIterator ConstSelfIterator ;
			typedef typename other::Iterator OtherIterator ;
			OtherIterator    Ap_i = Ap.Begin();
			ConstSelfIterator A_i = A.Begin();
			Hom<Field, _Tp1> hom(A. field(), Ap. field()) ;
			for ( ; A_i != A. End(); ++ A_i, ++ Ap_i)
				hom.image (*Ap_i, *A_i);
		}
    };    

    ///////////////////
    //      I/O      //
    ///////////////////



    template < class _Field, class _Storage >
    std::istream &BlasMatrix< _Field, _Storage >::read (std::istream &file)
    {
        MatrixStream<Field> ms(_field, file);
        if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) ){
            throw ms.reportError(__FUNCTION__,__LINE__);
        }
        _ptr = _rep.data();
        return file;
    }

    template < class _Field, class _Storage >
    std::ostream &BlasMatrix< _Field, _Storage >::write (std::ostream &os, LINBOX_enum (Tag::FileFormat) f) const
    {
        constSubMatrixType B(*this);
        return B.write(os,f);
    }
    


    ///////////////////
    //   ITERATORS   //
    ///////////////////

    template < class _Field, class _Storage >
    class BlasMatrix< _Field, _Storage >::ConstRowIterator {
    public:
        ConstRowIterator (const typename Storage::const_iterator& p, size_t len, size_t d) :
            _row (p, p + (ptrdiff_t)len), _dis (d)
        {}

        ConstRowIterator () {}

        ConstRowIterator (const ConstRowIterator& colp) :
            _row (colp._row), _dis (colp._dis)
        {}

        ConstRowIterator& operator = (const ConstRowIterator& colp)
        {
            _row = colp._row;
            _dis = colp._dis;
            return *this;
        }

        ConstRowIterator& operator --()
        {
            _row = ConstRow (_row.begin () - (ptrdiff_t)_dis, _row.end () - (ptrdiff_t)_dis);
            return *this;
        }

        ConstRowIterator  operator-- (int)
        {
            ConstRowIterator tmp (*this);
            --*this;
            return tmp;
        }


        ConstRowIterator& operator++ ()
        {
            _row = ConstRow (_row.begin () + (ptrdiff_t)_dis, _row.end () + (ptrdiff_t) _dis);
            return *this;
        }

        ConstRowIterator  operator++ (int)
        {
            ConstRowIterator tmp (*this);
            ++*this;
            return tmp;
        }

        ConstRowIterator operator+ (int i)
        {
            return ConstRowIterator (_row.begin () + (ptrdiff_t)((int)_dis * i), _row.size (), _dis);
        }

        ConstRowIterator& operator += (int i)
        {
            _row = ConstRow (_row.begin () + (ptrdiff_t)((int)_dis * i), _row.end () + (ptrdiff_t)((int)_dis * i));
            return *this;
        }

        ConstRow operator[] (int i) const
        {
            return ConstRow (_row.begin () + (ptrdiff_t)((int)_dis * i), _row.end () + (ptrdiff_t)((int)_dis * i));
        }

        ConstRow* operator-> ()
        {
            return &_row;
        }

        ConstRow& operator* ()
        {
            return _row;
        }

        bool operator!= (const ConstRowIterator& c) const
        {
            return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis);
        }

    private:
        ConstRow _row;
        size_t _dis;
    };

    template < class _Field, class _Storage >
    class BlasMatrix< _Field, _Storage >::RowIterator {
    public:
        RowIterator (const typename Storage::iterator& p, size_t len, size_t d) :
            _row (p, p + (ptrdiff_t)len), _dis (d)
        {}

        RowIterator () {}

        RowIterator (const RowIterator& colp) :
            _row (colp._row), _dis (colp._dis)
        {}

        RowIterator& operator = (const RowIterator& colp)
        {
            _row = colp._row;
            _dis = colp._dis;
            return *this;
        }

        RowIterator& operator ++ ()
        {
            _row = Row (_row.begin () + (ptrdiff_t)_dis, _row.end () + (ptrdiff_t)_dis);
            return *this;
        }

        RowIterator  operator ++ (int)
        {
            RowIterator tmp (*this);
            ++*this;
            return tmp;
        }

        RowIterator& operator -- ()
        {
            _row = Row (_row.begin () - (ptrdiff_t)_dis, _row.end () - (ptrdiff_t)_dis);
            return *this;
        }

        RowIterator  operator -- (int)
        {
            RowIterator tmp (*this);
            --*this;
            return tmp;
        }

        RowIterator operator + (int i)
        {
            return RowIterator (_row.begin () + (ptrdiff_t)((int)_dis * i), _row.size (), _dis);
        }

        RowIterator& operator += (int i)
        {
            _row = Row (_row.begin () + (ptrdiff_t)((int)_dis * i), _row.end () + (ptrdiff_t)((int)_dis * i));
            return *this;
        }

        Row operator[] (int i) const
        {
            return Row (const_cast<Row&> (_row).begin () + (ptrdiff_t)((int)_dis * i),
                        const_cast<Row&> (_row).end () + (ptrdiff_t)((int)_dis * i));
        }

        Row* operator-> ()
        {
            return &(this->_row);
        }

        Row& operator* ()
        {
            return _row;
        }

        bool operator!= (const RowIterator& c) const
        {
            return (_row.begin () != c._row.begin ()) || (_row.end () != c._row.end ()) || (_dis != c._dis);
        }

        operator ConstRowIterator ()
        {
            return ConstRowIterator (_row.begin (), _row.size (), _dis);
        }

    private:
        Row _row;
        size_t _dis;
    };


    
	template < class _Field, class _Storage >
	class BlasMatrix< _Field, _Storage >::ConstColIterator {
	public:
		ConstColIterator (typename _Storage::const_iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename _Storage::const_iterator> (p, (ptrdiff_t)stride),
			      Subiterator<typename _Storage::const_iterator> (p + (ptrdiff_t)(len * stride), (ptrdiff_t)stride)),
			_stride (stride)
		{}

		ConstColIterator (const ConstCol& col, size_t stride) :
			_col (col),
			_stride (stride)
		{}

		ConstColIterator () {}

		ConstColIterator (const ConstColIterator& rowp) :
			_col (rowp._col),
			_stride (rowp._stride)
		{}

		ConstColIterator& operator= (const ConstColIterator& rowp)
		{
			_col = rowp._col;
			_stride = rowp._stride;
			return *this;
		}

		ConstColIterator& operator++ ()
		{
			_col = ConstCol (Subiterator<typename _Storage::const_iterator> (_col.begin ().operator-> () + 1, (ptrdiff_t)_stride),
					 Subiterator<typename _Storage::const_iterator> (_col.end ().operator-> () + 1, (ptrdiff_t)_stride));
			return *this;
		}

		ConstColIterator  operator++ (int)
		{
			ConstColIterator old(*this);
			this->operator++ ();
			return old;
		}

		ConstColIterator operator + (int i)
		{
			return ConstColIterator (_col.begin ().operator-> () + i, _stride, _col.size ());
		}

		ConstColIterator& operator += (int i)
		{
			_col = ConstCol (Subiterator<typename _Storage::const_iterator> (_col.begin ().operator-> () + i, _stride),
					 Subiterator<typename _Storage::const_iterator> (_col.end ().operator-> () + i, _stride));
			return *this;
		}

		ConstCol operator[] (int i) const
		{
			return ConstCol (Subiterator<typename _Storage::const_iterator> (_col.begin ().operator-> () + i, _stride),
					 Subiterator<typename _Storage::const_iterator> (_col.end ().operator-> () + i, _stride));
		}

		ConstCol* operator-> ()
		{
			return &_col;
		}

		ConstCol& operator* ()
		{
			return _col;
		}

		bool operator!= (const ConstColIterator& c) const
		{
			return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ());
		}

	private:
		ConstCol _col;
		size_t _stride;
	};

	template < class _Field, class _Storage >
	class BlasMatrix< _Field, _Storage >::ColIterator {
	public:
		ColIterator (typename _Storage::iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename _Storage::iterator> (p, (long)stride),
			      Subiterator<typename _Storage::iterator> (p + (ptrdiff_t)(len * stride),(long) stride)), _stride (stride)
		{}

		ColIterator () {}

		ColIterator (const ColIterator& rowp) :
			_col (rowp._col)
		{}

		ColIterator& operator= (const ColIterator& rowp)
		{
			_col = rowp._col;
			_stride = rowp._stride;
			return *this;
		}

		const ColIterator& operator= (const ColIterator& rowp) const
		{
			const_cast<ColIterator*> (this)->_col = rowp._col;
			return *this;
		}

		ColIterator& operator++ ()
		{
			_col = Col (Subiterator<typename _Storage::iterator> (_col.begin ().operator-> () + 1, (const long)_stride),
				    Subiterator<typename _Storage::iterator> (_col.end ().operator-> () + 1,(const long) _stride));
			return *this;
		}

		ColIterator  operator++ (int)
		{
			Col tmp (_col);
			this->operator++ ();
			return tmp;
		}

		ColIterator operator + (int i)
		{
			return ColIterator (_col.begin ().operator-> () + i, _stride, _col.size ());
		}

		ColIterator& operator += (int i)
		{
			_col = Col (Subiterator<typename _Storage::iterator> (_col.begin ().operator-> () + i, _stride),
				    Subiterator<typename _Storage::iterator> (_col.end ().operator-> () + i, _stride));
			return *this;
		}

		Col operator[] (int i) const
		{
			return Col (Subiterator<typename _Storage::iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride),
				    Subiterator<typename _Storage::iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride));
		}

		Col* operator-> ()
		{
			return &(this->_col);
		}

		Col& operator* ()
		{
			return _col;
		}

		bool operator!= (const ColIterator& c) const
		{
			return (_col.begin () != c._col.begin ()) || (_col.end () != c._col.end ());
		}

		operator ConstColIterator ()
		{
			return ConstColIterator (reinterpret_cast<ConstCol&> (_col) , _stride);
		}

	private:

		Col _col;
		size_t _stride;
	};

   
    /*!   Indexed Iterator.
     * @ingroup iterators
     * @brief NO DOC
     */
    template < class _Field, class _Storage >
    class BlasMatrix< _Field, _Storage >::IndexedIterator {
        size_t _r_index;
        size_t _c_index;
        size_t _dim;
        typename Storage::iterator _begin;
        typedef typename _Field::Element value_type;

    public:
        IndexedIterator (const size_t  &dim,
                         const size_t  &r_index,
                         const size_t  &c_index,
                         const typename Storage::iterator &begin) :
            _r_index (r_index), _c_index (c_index), _dim (dim), _begin (begin)
        {}

        IndexedIterator () :
            _r_index (0), _c_index (0), _dim (1), _begin (0)
        {}

       
		IndexedIterator (const IndexedIterator& r) :
			_r_index (r._r_index), _c_index (r._c_index), _dim (r._dim), _begin (r._begin)
		{}

		IndexedIterator& operator = (const IndexedIterator &iter)
		{
			_r_index = iter._r_index;
			_c_index = iter._c_index;
			_dim = iter._dim;
			_begin = iter._begin;
			return *this;
		}

		bool operator == (const IndexedIterator &iter) const
		{
			return (_r_index == iter._r_index) &&
                (_c_index == iter._c_index) &&
                (_dim == iter._dim) &&
                (_begin==iter._begin);
		}

		bool operator != (const IndexedIterator& iter) const
		{
			return (_r_index != iter._r_index) ||
                (_c_index != iter._c_index) ||
                (_dim != iter._dim) ||
                (_begin!=iter._begin);
		}

		IndexedIterator &operator ++ ()
		{
			++_c_index;

			if (_c_index == _dim) {
				_c_index = 0;
				++_r_index;
			}

			return *this;
		}


		IndexedIterator operator ++ (int)
		{
			IndexedIterator tmp = *this;
			++(*this);
			return tmp;
		}

		IndexedIterator &operator -- ()
		{
			if (_c_index)
				--_c_index;
			else {
				--_r_index;
				_c_index = _dim - 1;
			}

			return *this;
		}


		IndexedIterator operator -- (int)
		{
			IndexedIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * () const
		{
			return *(_begin +(ptrdiff_t) (_r_index * _dim + _c_index));
		}


		value_type * operator -> () const
		{
			return _begin + (ptrdiff_t)(_r_index * _dim + _c_index);
		}


		size_t rowIndex () const
		{
			return _r_index;
		}

		size_t colIndex () const
		{
			return _c_index;
		}

		const value_type &value () const
		{
			return *(_begin + (ptrdiff_t)(_r_index * _dim + _c_index));
		}


	};
    
    template < class _Field, class _Storage >
    class BlasMatrix< _Field, _Storage >::ConstIndexedIterator {
        size_t _r_index;
        size_t _c_index;
        size_t _dim;
        typedef typename _Field::Element value_type;
        typename Storage::const_iterator _begin;

    public:
        ConstIndexedIterator (const size_t  &my_dim,
                              const size_t  &r_index,
                              const size_t  &c_index,
                              const typename Storage::const_iterator &begin) :
            _r_index (r_index), _c_index (c_index), _dim (my_dim), _begin (begin)
        {}

        ConstIndexedIterator () :
            _r_index (0), _c_index (0), _dim (1), _begin (0)
        {}

        ConstIndexedIterator (const ConstIndexedIterator& r) :
            _r_index (r._r_index), _c_index (r._c_index), _dim (r._dim), _begin (r._begin)
        {}

        ConstIndexedIterator& operator = (const ConstIndexedIterator &iter)
        {
            _r_index = iter._r_index;
            _c_index = iter._c_index;
            _dim = iter._dim;
            _begin = iter._begin;
            return *this;
        }

        bool operator == (const ConstIndexedIterator &iter) const
        {
            return (_r_index == iter._r_index) &&
                (_c_index == iter._c_index) &&
                (_dim == iter._dim) &&
                (_begin==iter._begin);
        }

        bool operator != (const ConstIndexedIterator& iter) const
        {
            return (_r_index != iter._r_index) ||
                (_c_index != iter._c_index) ||
                (_dim != iter._dim) ||
                (_begin!=iter._begin);
        }

        ConstIndexedIterator &operator ++ ()
        {
            ++_c_index;

            if (_c_index == _dim) {
                _c_index = 0;
                ++_r_index;
            }

            return *this;
        }


        ConstIndexedIterator operator ++ (int)
        {
            ConstIndexedIterator tmp = *this;
            ++(*this);
            return tmp;
        }

        ConstIndexedIterator &operator -- ()
        {
            if (_c_index)
                --_c_index;
            else {
                --_r_index;
                _c_index = _dim - 1;
            }

            return *this;
        }

        ConstIndexedIterator operator -- (int)
        {
            ConstIndexedIterator tmp = *this;
            --(*this);
            return tmp;
        }

        const value_type &operator * () const
        {
            return *(_begin + (ptrdiff_t)(_r_index * _dim + _c_index));
        }

        const value_type *operator -> () const
        {
            return _begin + (ptrdiff_t)(_r_index * _dim + _c_index);
        }

        size_t rowIndex () const
        {
            return _r_index;
        }

        size_t colIndex () const
        {
            return _c_index;
        }

        const value_type &value() const
        {
            return *(_begin + (ptrdiff_t)(_r_index * _dim + _c_index));
        }
    };


   

    /*   */

    // Entry access  view.  Size m*n vector in C (row major) order.
    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::Iterator BlasMatrix< _Field, _Storage >::Begin ()
    {
        return _rep.begin ();
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::Iterator BlasMatrix< _Field, _Storage >::End ()
    {
        return _rep.end ();
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstIterator BlasMatrix< _Field, _Storage >::Begin () const
    {
        return _rep.begin ();
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstIterator BlasMatrix< _Field, _Storage >::End () const
    {
        return _rep.end ();
    }

    /*   Indexed  */

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::IndexedIterator BlasMatrix< _Field, _Storage >::IndexedBegin ()
    {
        return IndexedIterator (coldim (), 0, 0, _rep.begin ());
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::IndexedIterator BlasMatrix< _Field, _Storage >::IndexedEnd ()
    {
        return IndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstIndexedIterator BlasMatrix< _Field, _Storage >::IndexedBegin () const
    {
        return ConstIndexedIterator (coldim (), 0, 0, _rep.begin ());
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstIndexedIterator BlasMatrix< _Field, _Storage >::IndexedEnd () const
    {
        return ConstIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
    }

    /*  Row  */

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::RowIterator BlasMatrix< _Field, _Storage >::rowBegin ()
    {
        return RowIterator (_rep.begin (), _col, _col);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::RowIterator BlasMatrix< _Field, _Storage >::rowEnd ()
    {
        return RowIterator (_rep.end (), _col, _col);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstRowIterator BlasMatrix< _Field, _Storage >::rowBegin () const
    {
        return ConstRowIterator (_rep.begin (), _col, _col);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstRowIterator BlasMatrix< _Field, _Storage >::rowEnd () const
    {
        return ConstRowIterator (_rep.end (), _col, _col);
    }

    /*  Col */

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ColIterator BlasMatrix< _Field, _Storage >::colBegin ()
    {
        return  typename BlasMatrix< _Field, _Storage >::ColIterator (_rep.begin (), _col, _row);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ColIterator BlasMatrix< _Field, _Storage >::colEnd ()
    {
        return  typename BlasMatrix< _Field, _Storage >::ColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstColIterator BlasMatrix< _Field, _Storage >::colBegin () const
    {
        return  typename BlasMatrix< _Field, _Storage >::ConstColIterator (_rep.begin (), _col, _row);
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstColIterator BlasMatrix< _Field, _Storage >::colEnd () const
    {
        return  typename BlasMatrix< _Field, _Storage >::ConstColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
    }

    /*  operators */
    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::Row BlasMatrix< _Field, _Storage >::operator[] (size_t i)
    {
        return Row (_rep.begin () +(ptrdiff_t)( i * _col), _rep.begin () + (ptrdiff_t)(i * _col +_col));
    }

    template < class _Field, class _Storage >
    typename BlasMatrix< _Field, _Storage >::ConstRow BlasMatrix< _Field, _Storage >::operator[] (size_t i) const
    {
        return Row (_rep.begin () +(ptrdiff_t) (i * _col), _rep.begin () + (ptrdiff_t)( i * _col + _col));
    }



} // end of LinBox namespace


#endif // __LINBOX_densematrix_blas_matrix_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
