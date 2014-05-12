/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
 *               Brice Boyer    <bboyer@imag.fr>
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
 * @file matrix/DenseMatrix/blas-matrix.inl
 * @ingroup densematrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
 */

#ifndef __LINBOX_densematrix_blas_matrix_INL
#define __LINBOX_densematrix_blas_matrix_INL

/////////////////
//   PRIVATE   //
/////////////////

namespace LinBox
{
	template<class _Field, class _Rep>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const BlasMatrix< _Field, _Rep >& A)
	{
#ifndef NDEBUG
		if (!areFieldEqual(A.field(),field())) {
			A.field().write(std::cout) << "!=" ;
			field().write(std::cout) << std::endl;
		}
#endif
		linbox_check( areFieldEqual(A.field(),field()));
		createBlasMatrix(A,MatrixContainerCategory::BlasContainer());
	}

	template<class _Field, class _Rep>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const Element * v)
	{
		// Element * iter_v = const_cast<Element*>(v) ;
		Element * v_end = const_cast<Element*>(v+(_col*_row)) ;
		// Iterator  iter_addr = this->Begin();
		Element * iter_addr = _ptr ;
		for (; v != v_end ; ++v, ++iter_addr)
		{
			field().init(*iter_addr);
			field().assign(*iter_addr,*v);
		}
	}

	template<class _Field, class _Rep>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const std::vector<Element> & v)
	{
		typename std::vector< Element>::const_iterator iter_value = v.begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != v.end(); ++iter_value,++iter_addr)
		{
			field().init(*iter_addr);
			field().assign(*iter_addr,*iter_value);
		}
	}


	template<class _Field, class _Rep>
	template <class _Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const _Matrix& A,
							   MatrixContainerCategory::BlasContainer)
	{
		linbox_check( areFieldEqual(A.field(),field()));
		typename _Matrix::ConstIterator         iter_value = A.Begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != A.End(); ++iter_value,++iter_addr)
		{
			field().init(*iter_addr);
			field().assign(*iter_addr, *iter_value);
		}
	}

	template<class _Field, class _Rep>
	template <class Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const Matrix& A,
							   MatrixContainerCategory::Container)
	{
		linbox_check( areFieldEqual(A.field(),field()));
		// const Field & F = A.field();
		//!@bug With both iterators, it is Segfaulting !!!!
		typename Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();
		for (;iter_index != A.IndexedEnd(); ++iter_index)
			setEntry( iter_index.rowIndex(),
				  iter_index.colIndex(),
				  A.getEntry(iter_index.rowIndex(),iter_index.colIndex())
				);
	}

	template<class _Field, class _Rep>
	template <class Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const Matrix& A,
							   MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field()) );

		BlasVector<Field> e(A.field(),A.coldim(), field().zero), tmp(A.field(),A.rowdim());
		ColIterator col_p;

		typename BlasMatrix< _Field, _Rep >::Col::iterator elt_p;
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

	template<class _Field, class _Rep>
	template <class _Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const _Matrix& A,
							   const size_t i0,const size_t j0,
							   const size_t m, const size_t n,
							   MatrixContainerCategory::Container)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );

		typename _Matrix::ConstIterator         iter_value = A.Begin();
		typename _Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

		for (;iter_value != A.End(); ++iter_value,++iter_index){
			int i,j;
			i=(int)iter_index.rowIndex()-(int)i0;
			j=(int)iter_index.colIndex()-(int)j0;
			if (( i >= 0)  && (j >= 0) && (i< (int)m) && (j < (int)n))
				setEntry(i, j, *iter_value);
		}
	}

	template<class _Field, class _Rep>
	template <class Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const Matrix& A,
							   const size_t i0,const size_t j0,
							   const size_t m, const size_t n,
							   MatrixContainerCategory::BlasContainer)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );

		typename Matrix::ConstIterator         iter_value = A.Begin();
		typename Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

		for (;iter_value != A.End(); ++iter_value,++iter_index){
			int i,j;
			i=(int)iter_index.rowIndex()-(int)i0;
			j=(int)iter_index.colIndex()-(int)j0;
			if ( (i>=0) && (j>=0) && (i< (int)m) && (j < (int)n))
				setEntry((size_t)i, (size_t)j, *iter_value);
		}
	}

	template<class _Field, class _Rep>
	template <class Matrix>
	void BlasMatrix< _Field, _Rep >::createBlasMatrix (const Matrix& A,
							   const size_t i0,const size_t j0,
							   const size_t m, const size_t n,
							   MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );


		BlasVector<Field> e(A.field(),A.coldim(), field().zero), tmp(A.field(),A.rowdim());
		ColIterator col_p;

		typename BlasMatrix< _Field, _Rep >::Col::iterator elt_p;
		typename BlasVector<Element>::iterator e_p, tmp_p;


		for (col_p = colBegin(), e_p = e.begin()+(ptrdiff_t)j0;
		     e_p != e.begin()+(ptrdiff_t)(j0+n); ++ col_p, ++ e_p) {

			field().assign(*e_p, field().one);

			A.apply (tmp, e);

			for (tmp_p = tmp.begin()+(ptrdiff_t)i0, elt_p = col_p -> begin();
			     elt_p != col_p.end(); ++ tmp_p, ++ elt_p) {
				field().assign(*elt_p, *tmp_p);
			}

			field().assign(*e_p, field().zero);
		}
	}

} // LinBox

//////////////////
// CONSTRUCTORS //
//////////////////

namespace LinBox
{


	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix (const _Field &F) :
		_row(0),_col(0),_rep(0)
		// ,_use_fflas(false)
		,_ptr(NULL)
		,_field(&F),_MD(F),_VD(F)
		// ,_AD(F)
	{
		// std::cout << "cstor 1 called" << std::endl;
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}

#if 0
	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix () :
		_row(0),_col(0),_rep(0),_ptr(NULL),
		_field(Field()),_MD(_field ),_VD(_field )
	{}
#endif

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::init(const _Field &F, const size_t & r, const size_t & c)
	{
		_field = &F; _row = r; _col = c;
		_rep.resize(r*c, F.zero);
		_ptr = &_rep[0];
		_VD.init(F); _MD.init(F);
		// _AD.init(F);
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix ( const _Field &F, const size_t & m, const size_t & n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
		// ,_AD(F)
	{
		// std::cout << "cstor 2 called" << std::endl;
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}


	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix(MatrixStream<_Field>& ms) :
		_row(0),_col(0),_rep(0),
		_field(&(ms.getField())),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 3 called" << std::endl;
		if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
			throw ms.reportError(__FUNCTION__,__LINE__);
		_ptr = &_rep[0];
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}

	template < class _Field, class _Rep >
	template <class StreamVector>
	BlasMatrix< _Field, _Rep >::BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) :
		_row(stream.size ()), _col(stream.dim ()), _rep(_row*_col), _ptr(&_rep[0]),
		_field (&F), _MD (F), _VD(F)
		// ,_AD(F)
	{
		// std::cout << "cstor 4 called" << std::endl;
		StreamVector tmp(F);
		typename BlasMatrix<Field,_Rep>::RowIterator p;

		VectorWrapper::ensureDim (tmp, stream.dim ());

		for (p = BlasMatrix<Field,_Rep>::rowBegin (); p != BlasMatrix<Field,_Rep>::rowEnd (); ++p) {
			stream >> tmp;
			_VD.copy (*p, tmp);
		}
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}


	template < class _Field, class _Rep >
	template <class Matrix>
	BlasMatrix< _Field, _Rep >::BlasMatrix (const Matrix &A) :
		_row(A.rowdim()),_col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 5 called" << std::endl;
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(A, typename MatrixContainerTrait<Matrix>::Type());
	}

	template < class _Field, class _Rep >
	template <class Matrix>
	BlasMatrix< _Field, _Rep >::BlasMatrix (const Matrix& A,
						const size_t &i0, const size_t &j0,
						const size_t &m,  const size_t &n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 6 called" << std::endl;
		_use_fflas = Protected::checkBlasApply(field(), _col);
		// makePointer();
		createBlasMatrix(A, i0, j0, m, n,
				 typename MatrixContainerTrait<Matrix>::Type());
	}

	template < class _Field, class _Rep >
	template<class _Matrix>
	BlasMatrix< _Field, _Rep >::BlasMatrix (const _Matrix &A,  const _Field &F) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 7 called" << std::endl;
		_use_fflas = Protected::checkBlasApply(field(), _col);
		typename _Matrix::template rebind<_Field>()(*this,A);
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix (const BlasMatrix< _Field, _Rep >& A) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 8 called" << std::endl;
		_use_fflas = Protected::checkBlasApply(field(), _col);
		// makePointer();
		createBlasMatrix(A);
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix (const _Field &F,
						const std::vector<typename _Field::Element>& v,
						const size_t & m, const size_t & n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F),_MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 9 called" << std::endl;
		linbox_check(v.size() == m*n);
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(v);
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::BlasMatrix (const _Field &F,
						const typename _Field::Element * v,
						const size_t & m, const size_t & n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F), _MD(field() ),_VD(field() )
		// ,_AD(field())
	{
		// std::cout << "cstor 10 called" << std::endl;
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(v);
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >::~BlasMatrix ()
	{
		// if (_ptr)
		// free(_ptr);
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::zero()
	{
		typename _Field::Element x; field().init(x);
		for (size_t i = 0; i < rowdim(); ++i)
			for (size_t j = 0; j < coldim(); ++j)
				setEntry(i, j, field().zero);
	}

} // LinBox

///////////////////
//      I/O      //
///////////////////

namespace LinBox
{

	template < class _Field, class _Rep >
	std::istream &BlasMatrix< _Field, _Rep >::read (std::istream &file)
	{
		MatrixStream<Field> ms(field(), file);
		if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
			throw ms.reportError(__FUNCTION__,__LINE__);
		_ptr = &_rep[0];
		_use_fflas = Protected::checkBlasApply(field(), _col);
		return file;
	}

	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep >& BlasMatrix< _Field, _Rep >::operator= (const BlasMatrix< _Field, _Rep >& A)
	{
		if ( &A == this)
			return *this;

		_col = A.coldim();
		_row = A.rowdim();
		_rep = Rep(_row*_col);
		_ptr = &_rep[0] ;
		// const_cast<_Field&>(_field ) = A.field();
		// changeField( A.field() );
		createBlasMatrix(A);

		return *this;
	}

#if 0 /*  loop rebind */
	template < class _Field, class _Rep >
	template<typename _Tp1>
	struct BlasMatrix< _Field, _Rep >::rebind {
		typedef BlasMatrix<_Tp1> other;

		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
			// typedef typename BlasMatrix< _Field, _Rep >::ConstIndexedIterator ConstIndexedIterator ;
			// typedef typename BlasMatrix< _Field, _Rep >::ConstIterator ConstIterator ;
			ConstIterator         iter_value = A.Begin();
			ConstIndexedIterator  iter_index = A.IndexedBegin();
			typename _Tp1::Element tmp;
			for (;iter_value != A.End(); ++iter_value,++iter_index){
				Ap.field().init(tmp);
				Ap.field().assign(tmp, *iter_value);
				Ap.setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
			}
		}
	};
#endif

#if 1 /*  HOM */
	//! @bug other rep
	template < class _Field, class _Rep >
	template<typename _Tp1>
	struct BlasMatrix< _Field, _Rep >::rebind {
		typedef BlasMatrix<_Tp1,typename Vector<_Tp1>::Dense> other;

		void operator() (other & Ap, const Self_t& A) {
			// typedef Self_t::ConstIterator ConstSelfIterator ;
			typedef typename BlasMatrix< _Field, _Rep >::ConstIterator ConstSelfIterator ;
			typedef typename other::Iterator OtherIterator ;
			OtherIterator    Ap_i = Ap.Begin();
			ConstSelfIterator A_i = A.Begin();
			Hom<Field, _Tp1> hom(A. field(), Ap. field()) ;
			for ( ; A_i != A. End(); ++ A_i, ++ Ap_i)
				hom.image (*Ap_i, *A_i);
		}
	};
#endif


} // LinBox

//////////////////
//  DIMENSIONS  //
//////////////////

namespace LinBox
{

	template < class _Field, class _Rep >
	size_t BlasMatrix< _Field, _Rep >::rowdim() const
	{
		return _row ;
	}

	template < class _Field, class _Rep >
	size_t BlasMatrix< _Field, _Rep >::coldim() const
	{
		return _col ;
	}

	template < class _Field, class _Rep >
	size_t  BlasMatrix< _Field, _Rep >::getStride() const
	{
		return _col;
	}

	template < class _Field, class _Rep >
	size_t&  BlasMatrix< _Field, _Rep >::getWriteStride()
	{
		return _col;
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::resize (const size_t & m, const size_t & n, const Element& val )
	{
#ifndef NDEBUG
		if (_col > 0 && _col != n)
			std::cerr << " ***Warning*** you are resizing a matrix, possibly loosing data. " << std::endl;
#endif
		_row = m;
		_col = n;
		_rep.resize (m * n, val);
		_ptr = &_rep[0];
#if 0
		if (_ptr) {
			if (m && n)
				realloc(_ptr,m*n*sizeof(Element));
			else {
				free(_ptr);
				_ptr=NULL ;
			}
		}
		else
			makePointer();
#endif
	}

} // LinBox

//////////////////
//   ELEMENTS   //
//////////////////

namespace LinBox
{


	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::pointer
	BlasMatrix< _Field, _Rep >::getPointer() const
	{
		return _ptr;
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::const_pointer &
	BlasMatrix< _Field, _Rep >::getConstPointer() const
	{
		return (const_pointer)_ptr;
	}


	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::pointer&
	BlasMatrix< _Field, _Rep >::getWritePointer()
	{
		return _ptr;
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::setEntry (size_t i, size_t j, const Element &a_ij)
	{
		_ptr[i * _col + j] = a_ij;
	}

	template < class _Field, class _Rep >
	typename _Field::Element & BlasMatrix< _Field, _Rep >::refEntry (size_t i, size_t j)
	{
		return _ptr[i * _col + j];
	}

	template < class _Field, class _Rep >
	const typename _Field::Element & BlasMatrix< _Field, _Rep >::getEntry (size_t i, size_t j) const
	{
		return _ptr[i * _col + j];
	}

	template < class _Field, class _Rep >
	typename _Field::Element & BlasMatrix< _Field, _Rep >::getEntry (Element &x, size_t i, size_t j) const
	{
		x = _ptr[i * _col + j];
		return x;
	}

} // LinBox

///////////////////
// TRANSPOSE &AL //
///////////////////

namespace LinBox
{
	template < class _Field, class _Rep >
	BlasMatrix< _Field, _Rep > BlasMatrix< _Field, _Rep >::transpose(BlasMatrix< _Field, _Rep > & tM) const
	{
		size_t r = this->rowdim() ;
		size_t c = this->coldim() ;
		linbox_check(tM.coldim() == r );
		linbox_check(tM.rowdim() == c);
		for (size_t i = 0 ; i < r ; ++i)
			for (size_t j = 0 ; j < c ; ++j)
				tM.setEntry(j,i,this->getEntry(i,j));
		return tM;
	}

	namespace Protected
	{
		/*! @internal
		 *  @brief In-Place Tranpose.
		 * Adapted from the Wikipedia article.
		 * @todo make it for strides and Submatrices
		 * @todo use specialized versions when available (eg dgetmi)
		 * @todo make transpose have an inplace=true default parameter
		 * (execpt maybe when below a threshold).
		 * @param m pointer to the beginning of a row-major matrix vector
		 * @param w rows in the matrix
		 * @param h cols in the matrix
		 */
		template<class T>
		void transposeIP(T *m, size_t h, size_t w)
		{
			size_t start;
			T tmp;

			for (start = 0; start <= w * h - 1; ++start) {
				size_t next, i ;
				next = start;
				i = 0;
				do {
					++i;
					next = (next % h) * w + next / h;
				} while (next > start);
				if (next < start || i == 1)
					continue;

				tmp = m[next = start];
				do {
					i = (next % h) * w + next / h;
					m[next] = (i == start) ? tmp : m[i];
					next = i;
				} while (next > start);
			}
		}
	} // Protected

	template < class _Field, class _Rep >
	template<bool _IP>
	void BlasMatrix< _Field, _Rep >::transpose()
	{
		size_t r = this->rowdim() ;
		size_t c = this->coldim() ;

		if ( r == c) {
			for (size_t i = 0 ; i < r ; ++i)
				for (size_t j = i+1 ; j < c ; ++j)
					std::swap(this->refEntry(i,j),this->refEntry(j,i));
		}
		else {
			if (!_IP) {
				BlasMatrix< _Field, _Rep > MM(*this);
				std::swap(_row,_col);
				// iterating row first seems faster.
#ifdef _BLOCKIT
				size_t B ;
				B =  1024;

				for (size_t bi = 0 ; bi < r/B ; ++bi) {
					for (size_t bj = 0 ; bj < c/B ; ++bj){
						for (size_t i = 0 ; i < B ; ++i)
							for (size_t j = 0 ; j < B ; ++j)
								this->setEntry(bj*B+j,bi*B+i,
									       MM.getEntry(bi*B+i,bj*B+j));
					}
				}
				for (size_t i = r/B*B ; i < r ; ++i)
					for (size_t j = c/B*B ; j < c ; ++j)
						this->setEntry(j,i,
							       MM.getEntry(i,j));
#else
				for (size_t i = 0 ; i < r ; ++i)
					for (size_t j = 0 ; j < c ; ++j)
						this->setEntry(j,i,
							       MM.getEntry(i,j));
#endif

			}
			else {
				Protected::transposeIP<Element>(_ptr,_row,_col);
				std::swap(_row,_col);
			}
		}
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::transpose()
	{
		this->transpose<false>();
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::reverseRows()
	{
		size_t r = this->rowdim()/2 ;
		for (size_t i = 0 ; i <  r ; ++i) {
			_VD.swap( this->rowBegin()+i,
				  this->rowBegin()+(r-1-i) );
		}

	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::reverseCols()
	{
		size_t r = this->rowdim() ;
		size_t c = this->coldim()/2 ;
		for (size_t j = 0 ; j < c ; ++j) {
			for (size_t i = 0 ; i < r ; ++i) {
				std::swap(this->refEntry(i,j),
					  this->refEntry(i,c-j-1));

			}
		}
	}

	template < class _Field, class _Rep >
	void BlasMatrix< _Field, _Rep >::reverse()
	{
		size_t r = this->rowdim() ;
		size_t c = this->coldim() ;
		for (size_t j = 0 ; j < c ; ++j) {
			for (size_t i = 0 ; i < r ; ++i) {
				std::swap(this->refEntry(i,j),
					  this->refEntry(r-i-1,c-j-1));

			}
		}
	}
} // LinBox

///////////////////
//   ITERATORS   //
///////////////////

namespace LinBox
{

	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::ConstRowIterator {
	public:
		ConstRowIterator (const typename Rep::const_iterator& p, size_t len, size_t d) :
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

	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::RowIterator {
	public:
		RowIterator (const typename Rep::iterator& p, size_t len, size_t d) :
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

	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::ConstColIterator {
	public:
		ConstColIterator (typename Rep::const_iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename Rep::const_iterator> (p, (ptrdiff_t)stride),
			      Subiterator<typename Rep::const_iterator> (p + (ptrdiff_t)(len * stride), (ptrdiff_t)stride)),
			_stride (stride)
		{}

		ConstColIterator (const ConstCol& col, size_t stride) :
			_col (col),
			_stride (stride)
		{}

		ConstColIterator () {}

		ConstColIterator (const ConstColIterator& rowp) :
			_col (rowp._col)
		{}

		ConstColIterator& operator= (const ConstColIterator& rowp)
		{
			_col = rowp._col;
			_stride = rowp._stride;
			return *this;
		}

		ConstColIterator& operator++ ()
		{
			_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + 1, (ptrdiff_t)_stride),
					 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + 1, (ptrdiff_t)_stride));
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
			_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride),
					 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride));
			return *this;
		}

		ConstCol operator[] (int i) const
		{
			return ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + i, _stride),
					 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + i, _stride));
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

	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::ColIterator {
	public:
		ColIterator (typename Rep::iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename Rep::iterator> (p, (long)stride),
			      Subiterator<typename Rep::iterator> (p + (ptrdiff_t)(len * stride),(long) stride)), _stride (stride)
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
			_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + 1, (const long)_stride),
				    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + 1,(const long) _stride));
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
			_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + i, _stride),
				    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + i, _stride));
			return *this;
		}

		Col operator[] (int i) const
		{
			return Col (Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).begin ().operator-> () + i, _stride),
				    Subiterator<typename Rep::iterator> (const_cast<Col&> (_col).end ().operator-> () + i, _stride));
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
	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::IndexedIterator {
		size_t _r_index;
		size_t _c_index;
		size_t _dim;
		typename Rep::iterator _begin;
		typedef typename _Field::Element value_type;

	public:
		IndexedIterator (const size_t  &dim,
				 const size_t  &r_index,
				 const size_t  &c_index,
				 const typename Rep::iterator &begin) :
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

	template < class _Field, class _Rep >
	class BlasMatrix< _Field, _Rep >::ConstIndexedIterator {
		size_t _r_index;
		size_t _c_index;
		size_t _dim;
		typedef typename _Field::Element value_type;
		typename Rep::const_iterator _begin;

	public:
		ConstIndexedIterator (const size_t  &my_dim,
				      const size_t  &r_index,
				      const size_t  &c_index,
				      const typename Rep::const_iterator &begin) :
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
	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::Iterator BlasMatrix< _Field, _Rep >::Begin ()
	{
		return _rep.begin ();
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::Iterator BlasMatrix< _Field, _Rep >::End ()
	{
		return _rep.end ();
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstIterator BlasMatrix< _Field, _Rep >::Begin () const
	{
		return _rep.begin ();
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstIterator BlasMatrix< _Field, _Rep >::End () const
	{
		return _rep.end ();
	}

	/*   Indexed  */

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::IndexedIterator BlasMatrix< _Field, _Rep >::IndexedBegin ()
	{
		return IndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::IndexedIterator BlasMatrix< _Field, _Rep >::IndexedEnd ()
	{
		return IndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstIndexedIterator BlasMatrix< _Field, _Rep >::IndexedBegin () const
	{
		return ConstIndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstIndexedIterator BlasMatrix< _Field, _Rep >::IndexedEnd () const
	{
		return ConstIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	/*  Row  */

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::RowIterator BlasMatrix< _Field, _Rep >::rowBegin ()
	{
		return RowIterator (_rep.begin (), _col, _col);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::RowIterator BlasMatrix< _Field, _Rep >::rowEnd ()
	{
		return RowIterator (_rep.end (), _col, _col);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstRowIterator BlasMatrix< _Field, _Rep >::rowBegin () const
	{
		return ConstRowIterator (_rep.begin (), _col, _col);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstRowIterator BlasMatrix< _Field, _Rep >::rowEnd () const
	{
		return ConstRowIterator (_rep.end (), _col, _col);
	}

	/*  Col */

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ColIterator BlasMatrix< _Field, _Rep >::colBegin ()
	{
		return  typename BlasMatrix< _Field, _Rep >::ColIterator (_rep.begin (), _col, _row);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ColIterator BlasMatrix< _Field, _Rep >::colEnd ()
	{
		return  typename BlasMatrix< _Field, _Rep >::ColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstColIterator BlasMatrix< _Field, _Rep >::colBegin () const
	{
		return  typename BlasMatrix< _Field, _Rep >::ConstColIterator (_rep.begin (), _col, _row);
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstColIterator BlasMatrix< _Field, _Rep >::colEnd () const
	{
		return  typename BlasMatrix< _Field, _Rep >::ConstColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
	}

	/*  operators */
	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::Row BlasMatrix< _Field, _Rep >::operator[] (size_t i)
	{
		return Row (_rep.begin () +(ptrdiff_t)( i * _col), _rep.begin () + (ptrdiff_t)(i * _col +_col));
	}

	template < class _Field, class _Rep >
	typename BlasMatrix< _Field, _Rep >::ConstRow BlasMatrix< _Field, _Rep >::operator[] (size_t i) const
	{
		return Row (_rep.begin () +(ptrdiff_t) (i * _col), _rep.begin () + (ptrdiff_t)( i * _col + _col));
	}

} // LinBox

//////////////////
//     MISC     //
//////////////////

namespace LinBox
{
	template < class _Field, class _Rep >
	template <class Vector>
	Vector& BlasMatrix< _Field, _Rep >::columnDensity (Vector &v) const
	{
		std::fill (v.begin (), v.end (), _row);
		return v;
	}

} // LinBox

///////////////////
//   BLACK BOX   //
///////////////////

namespace LinBox
{
#if 1
	template < class _Field, class _Rep >
	template <class Vector1, class Vector2>
	Vector1&  BlasMatrix< _Field, _Rep >::apply (Vector1& y, const Vector2& x) const
	{
		// std::cout << "prepare apply1 Matrix" << std::endl;
		BlasSubmatrix<constSelf_t> A(*this);
		// std::cout << "...................." << std::endl;
		A.apply(y,x);
		// std::cout << ".......done........." << std::endl;
		return y;
	}
#endif

#if 1
	template < class _Field, class _Rep >
	template<class _Vrep>
	BlasVector<_Field,_Vrep>&  BlasMatrix< _Field, _Rep >::apply (BlasVector<_Field,_Vrep>& y, const BlasVector<_Field,_Vrep>& x) const
	{
		// std::cout << "prepare apply2 Matrix" << std::endl;
		BlasSubmatrix<constSelf_t> A(*this);
		// std::cout << "...................." << std::endl;
		A.apply(y,x);
		// std::cout << ".......done........." << std::endl;

		return y;
	}
#endif

	template < class _Field, class _Rep >
	template <class Vector1, class Vector2>
	Vector1&  BlasMatrix< _Field, _Rep >::applyTranspose (Vector1& y, const Vector2& x) const
	{
		// std::cout << "prepare applyT Matrix" << std::endl;
		BlasSubmatrix<constSelf_t> A(*this);
		// std::cout << "...................." << std::endl;
		A.applyTranspose(y,x);
		// std::cout << ".......done........." << std::endl;

		return y;
	}

	template < class _Field, class _Rep >
	const _Field& BlasMatrix< _Field, _Rep >::field() const
	{
		return *_field;
	}

#if 0 /* why not ? */
	template < class _Field, class _Rep >
	_Field& BlasMatrix< _Field, _Rep >::field()
	{
		return const_cast<_Field&>(_field );
	}
#endif
}


#endif // __LINBOX_densematrix_blas_matrix_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
