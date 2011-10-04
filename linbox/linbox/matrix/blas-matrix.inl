/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/matrix/blas-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *               Clément Pernet <clement.pernet@imag.fr>
 *               Brice Boyer    <bboyer@imag.fr>
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

/*!@internal
 * @file matrix/blas-matrix.inl
 * @ingroup matrix
 * A \c BlasMatrix<\c _element > represents a matrix as an array of
 * <code>_element</code>s.
 */

#ifndef __LINBOX_blas_matrix_INL
#define __LINBOX_blas_matrix_INL

/////////////////
//   PRIVATE   //
/////////////////

namespace LinBox
{
	template<class _Element>
	void BlasMatrix<_Element>::createBlasMatrix (const BlasMatrix<_Element>& A)
	{
		createBlasMatrix(A,MatrixContainerCategory::BlasContainer());
	}

	template<class _Element>
	void BlasMatrix<_Element>::createBlasMatrix (const _Element * v)
	{
		const_pointer iter_v = v ;
		const_pointer v_end = v+(_col*_row) ;
		Iterator  iter_addr = this->Begin();
		for (; v != v_end ; ++v, ++iter_addr)
			*iter_addr = *v;
	}

	template<class _Element>
	void BlasMatrix<_Element>::createBlasMatrix (const std::vector< _Element> & v)
	{
		typename std::vector< _Element>::const_iterator iter_value = v.begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != v.end(); ++iter_value,++iter_addr)
			*iter_addr = *iter_value;
	}


	template<class _Element>
	template <class _Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const _Matrix& A,
						     MatrixContainerCategory::BlasContainer)
	{
		typename _Matrix::ConstIterator         iter_value = A.Begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != A.End(); ++iter_value,++iter_addr)
			*iter_addr = *iter_value;
	}

	template<class _Element>
	template <class Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const Matrix& A,
						     MatrixContainerCategory::Container)
	{
		//!@bug With both iterators, it is Segfaulting !!!!
		typename Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();
		for (;iter_index != A.IndexedEnd(); ++iter_index)
			setEntry( iter_index.rowIndex(),
				  iter_index.colIndex(),
				  A.getEntry(iter_index.rowIndex(),iter_index.colIndex())
				);
	}

	template<class _Element>
	template <class Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const Matrix& A,
						     MatrixContainerCategory::Blackbox)
	{
		typedef typename Matrix::Field Field;
		typename Field::Element one, zero;
		Field F = A.field();
		F. init(one, 1);
		F. init(zero, 0);

		std::vector<typename Field::Element> e(A.coldim(), zero), tmp(A.rowdim());
		ColIterator col_p;

		typename BlasMatrix< _Element>::Col::iterator elt_p;

		typename std::vector<typename Field::Element>::iterator e_p, tmp_p;


		for (col_p = colBegin(), e_p = e.begin();
		     e_p != e.end(); ++ col_p, ++ e_p)
		{

			F.assign(*e_p, one);

			A.apply (tmp, e);

			for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
			     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)

				F.assign(*elt_p, *tmp_p);

			F.assign(*e_p, zero);
		}
	}

	template<class _Element>
	template <class _Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const _Matrix& A,
						     const size_t i0,const size_t j0,
						     const size_t m, const size_t n,
						     MatrixContainerCategory::Container)
	{

		typename _Matrix::ConstIterator         iter_value = A.Begin();
		typename _Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

		for (;iter_value != A.End(); ++iter_value,++iter_index){
			size_t i,j;
			i=iter_index.rowIndex();
			j=iter_index.colIndex();
			if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
				setEntry(i-i0, j-j0, *iter_value);
		}
	}

	template<class _Element>
	template <class Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const Matrix& A,
						     const size_t i0,const size_t j0,
						     const size_t m, const size_t n,
						     MatrixContainerCategory::BlasContainer)
	{

		typename Matrix::ConstIterator         iter_value = A.Begin();
		typename Matrix::ConstIndexedIterator  iter_index = A.IndexedBegin();

		for (;iter_value != A.End(); ++iter_value,++iter_index){
			size_t i,j;
			i=iter_index.rowIndex();
			j=iter_index.colIndex();
			if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
				setEntry(i-i0, j-j0, *iter_value);
		}
	}

	template<class _Element>
	template <class Matrix>
	void BlasMatrix<_Element>::createBlasMatrix (const Matrix& A,
						     const size_t i0,const size_t j0,
						     const size_t m, const size_t n,
						     MatrixContainerCategory::Blackbox)
	{
		throw(NotImplementedYet(__func__,__FILE__,__LINE__,
					"need to be implemented by succesive apply"));
	}

} // LinBox

//////////////////
// CONSTRUCTORS //
//////////////////

namespace LinBox
{


	template <class _Element>
	BlasMatrix< _Element>::BlasMatrix () :
		_row(0),_col(0),_rep(0),_ptr(NULL)
	{ }



	template <class _Element>
	template<class T>
	BlasMatrix< _Element>::BlasMatrix ( unsigned int m, T n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		// linbox_check(n>=0);
		// makePointer();
	}

	template <class _Element>
	template<class T>
	BlasMatrix< _Element>::BlasMatrix (long m, T n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		linbox_check(n>=0);
		linbox_check(m>=0);
		// makePointer();
	}

	template <class _Element>
	template<class T>
	BlasMatrix< _Element>::BlasMatrix (unsigned long m, T  n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		//!@todo
		// linbox_check_non_neg(n);
		// linbox_check(n>=0);
		// makePointer();
	}

	template <class _Element>
	template<class T>
	BlasMatrix< _Element>::BlasMatrix (int m, T n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		linbox_check(n>=0);
		linbox_check(m>=0);
		// makePointer();
	}

	template <class _Element>
	template<class T>
	BlasMatrix< _Element>::BlasMatrix ( Integer & m, T n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		//!@todo check m,n not too big ?
		linbox_check(n>=0);
		linbox_check(m>=0);
		// makePointer();
	}



	template <class _Element>
	template< class Field >
	BlasMatrix< _Element>::BlasMatrix(MatrixStream<Field>& ms) :
		_row(0),_col(0),_rep(0)
	{
		if( !ms.getArray(_rep) || !ms.getRows(_row) || !ms.getColumns(_col) )
			throw ms.reportError(__FUNCTION__,__LINE__);
		_ptr = &_rep[0];
	}

	template <class _Element>
	template <class Matrix>
	BlasMatrix< _Element>::BlasMatrix (const Matrix &A) :
		_row(A.rowdim()),_col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0])
	{
		// makePointer();
		createBlasMatrix(A, typename MatrixContainerTrait<Matrix>::Type());
	}

	template <class _Element>
	template <class Matrix>
	BlasMatrix< _Element>::BlasMatrix (const Matrix& A,
					 const size_t i0, const size_t j0,
					 const size_t m, const size_t n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		// makePointer();
		createBlasMatrix(A, i0, j0, m, n,
				 typename MatrixContainerTrait<Matrix>::Type());
	}

	template <class _Element>
	template<class _Matrix, class _Field>
	BlasMatrix< _Element>::BlasMatrix (const _Matrix &A,  const _Field &F) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0])
	{
		// makePointer();
		typename _Matrix::template rebind<_Field>()(*this,A,F);
	}

	template <class _Element>
	BlasMatrix< _Element>::BlasMatrix (const BlasMatrix< _Element>& A) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0])
	{
		// makePointer();
		createBlasMatrix(A);
	}

	template <class _Element>
	BlasMatrix< _Element>::BlasMatrix (const std::vector< _Element>& v, size_t m, size_t n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		linbox_check(v.size() == m*n);
		// makePointer();
		createBlasMatrix(v);
	}

	template <class _Element>
	BlasMatrix< _Element>::BlasMatrix (const _Element * v, size_t m, size_t n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0])
	{
		// makePointer();
		createBlasMatrix(v);
	}

	template <class _Element>
	BlasMatrix< _Element>::~BlasMatrix ()
	{
		// if (_ptr)
		// free(_ptr);
	}

} // LinBox

///////////////////
//      I/O      //
///////////////////

namespace LinBox
{

	template <class _Element>
	template <class Field>
	std::istream &BlasMatrix< _Element>::read (std::istream &file, const Field &F)
	{
#if 0
		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;

		if (m*n < _row*_col)
			cerr<<"NOT ENOUGH ELEMENT TO READ\n";
		else {
			for (p = Begin (); p != End (); ++p) {
				integer tmp;
				file>>tmp;cout<<tmp<<endl;
				//file.ignore(1);
				F.read (file, *p);
			}
		}
#endif


		Iterator p;
		int m,n;
		char c;
		file>>m>>n>>c;
		// std::cout << m << 'x' << n << ':' << c << std::endl;
		_row = m; _col = n;

	 _Element zero;
		F.init(zero,0UL);
		// _ptr.resize(_row * _col, zero);
		resize(_row,_col);

		if ((c != 'M') && (c != 'm')) {
			for (p = Begin (); p != End (); ++p) {
				//file.ignore(1);
				F.read (file, *p);
			}

		}
		else { // sparse file format - needs fixing
			int i, j;
			while (true)
			{
				file >> i >> j;
				//file.ignore(1);
				//if (! file) break;
				if (i+j <= 0) break;
				// std::cout << i << ',' << j << ':' ;
				F.read (file, _ptr[_col*(i-1) + j-1]);
				// std::cout << _ptr[_col*(i-1) + j-1] << std::endl;
			}
		}

		return file;
	}

	template <class _Element>
	template <class Field>
	std::ostream& BlasMatrix< _Element>::write (std::ostream &os, const Field &F,
						  bool mapleFormat) const
	{

		ConstRowIterator p;

		if (!mapleFormat) {
			integer c;
			int wid;




			F.cardinality (c);

			if (c >0)
				wid = (int) ceil (log ((double) c) / M_LN10);
			else {
				integer tmp;
				size_t max=0;
				ConstIterator it = Begin();
				for (; it != End(); ++it){
					F.convert(tmp,*it);
					if (tmp.bitsize() > max)
						max= tmp.bitsize();
				}
				wid= (int) ceil ((double)max / M_LN10)+1;
			}

			for (p = rowBegin (); p != rowEnd ();++p) {
				typename ConstRow::const_iterator pe;

				os << "  [ ";

				for (pe = p->begin (); pe != p->end (); ++pe) {
					os.width (wid);
					/*!  @warning
					 * matrix base does not provide this field(), maybe should?
					 * _M.field ().write (os, *pe);
					 * os << *pe;
					 * fixed by using extra field
					 */

					F.write (os, *pe);
					os << " ";
				}

				os << "]" << std::endl;
			}
		}
		else {

			os << "Matrix( " << rowdim() << ',' << coldim() << ",[" ;
			for (p = rowBegin (); p != rowEnd (); ) {
				typename ConstRow::const_iterator pe;

				os << " [ ";

				for (pe = p->begin (); pe != p->end (); ) {
					F.write (os, *pe);
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
		return os;
	}

	template <class _Element>
	std::ostream& BlasMatrix< _Element>::write (std::ostream &os,
						  bool mapleFormat) const
	{

		ConstRowIterator p;
		if (!mapleFormat) {
			integer c;
			int wid;




			integer tmp;
			size_t max=0;
			ConstIterator it = Begin();
			for (; it != End(); ++it){
				tmp = (integer) *it;
				if (tmp.bitsize() > max)
					max= tmp.bitsize();
			}
			wid= (int) ceil ((double)max / M_LN10)+1;

			for (p = rowBegin (); p != rowEnd (); ++p) {
				typename ConstRow::const_iterator pe;

				os << "  [ ";

				for (pe = p->begin (); pe != p->end (); ++pe) {
					os.width (wid);
					os << *pe;
					os << " ";
				}

				os << "]" << std::endl;
			}
		}
		else {

			os << "Matrix( " << rowdim() << ',' << coldim() << ",[" ;
			for (p = rowBegin (); p != rowEnd (); ) {
				typename ConstRow::const_iterator pe;

				os << " [ ";

				for (pe = p->begin (); pe != p->end (); ) {
					os << *pe;
					++pe ;
					if (pe != p->end())
						os << ", ";
				}

				os << "]" ;
				++p ;
				if (p != rowEnd() )
					os << ',' << std::endl;
			}
			os << "])" ;
		}

		return os;
	}


	template <class _Element>
	BlasMatrix< _Element>& BlasMatrix< _Element>::operator= (const BlasMatrix< _Element>& A)
	{

		_col = A.coldim();
		_row = A.rowdim();
		_rep = Rep(_row*_col);
		_ptr = &_rep[0] ;
		// makePointer();
		createBlasMatrix(A);

		return *this;
	}

	template <class _Element>
	template<typename _Tp1>
	struct BlasMatrix< _Element>::rebind {
		typedef BlasMatrix<typename _Tp1::Element> other;

		void operator() (other & Ap, const Self_t& A, const _Tp1& F)
		{
#if 0
			// typedef typename BlasMatrix< _Element>::ConstIndexedIterator ConstIndexedIterator ;
			// typedef typename BlasMatrix< _Element>::ConstIterator ConstIterator ;
			ConstIterator         iter_value = A.Begin();
			ConstIndexedIterator  iter_index = A.IndexedBegin();
			typename _Tp1::Element tmp;
			for (;iter_value != A.End(); ++iter_value,++iter_index){
				F.init(  tmp, *iter_value );
				Ap.setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
			}
#endif
#if 1
			linbox_check( Ap.rowdim() == A.rowdim() );
			linbox_check( Ap.coldim() == A.coldim() ) ;
			for (size_t i = 0 ; i < A.rowdim() ; ++i)
				for (size_t j = 0 ; j < A.coldim() ; ++j) {
					F.init(Ap.refEntry(i,j),(Element)A.getEntry(i,j));
				}
#endif
			return ;
		}
	};

} // LinBox

//////////////////
//  DIMENSIONS  //
//////////////////

namespace LinBox
{

	template <class _Element>
	size_t BlasMatrix< _Element>::rowdim() const
	{
		return _row ;
	}

	template <class _Element>
	size_t BlasMatrix< _Element>::coldim() const
	{
		return _col ;
	}

	template <class _Element>
	size_t  BlasMatrix< _Element>::getStride() const
	{
		return _col;
	}

	template <class _Element>
	size_t&  BlasMatrix< _Element>::getWriteStride()
	{
		return _col;
	}

	template <class _Element>
	void BlasMatrix< _Element>::resize (size_t m, size_t n, const _Element& val )
	{
#ifndef NDEBUG
		if (_col != n)
			std::cerr << " ***Warning*** you are resizing a matrix, possibly loosing data. " << std::endl;
#endif
		_row = m;
		_col = n;
		_rep.resize (m * n, val);
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


	template <class _Element>
	typename BlasMatrix< _Element>::pointer
	BlasMatrix< _Element>::getPointer() const
	{
		return _ptr;
	}

	template <class _Element>
	typename BlasMatrix< _Element>::const_pointer &
	BlasMatrix< _Element>::getConstPointer() const
	{
		return (const_pointer)_ptr;
	}


	template <class _Element>
	typename BlasMatrix< _Element>::pointer&
	BlasMatrix< _Element>::getWritePointer()
	{
		return _ptr;
	}

	template <class _Element>
	void BlasMatrix< _Element>::setEntry (size_t i, size_t j, const _Element &a_ij)
	{
		_ptr[i * _col + j] = a_ij;
	}

	template <class _Element>
 _Element & BlasMatrix< _Element>::refEntry (size_t i, size_t j)
	{
		return _ptr[i * _col + j];
	}

	template <class _Element>
	const _Element & BlasMatrix< _Element>::getEntry (size_t i, size_t j) const
	{
		return _ptr[i * _col + j];
	}

	template <class _Element>
 _Element & BlasMatrix< _Element>::getEntry (Element &x, size_t i, size_t j) const
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
	template <class _Element>
	BlasMatrix< _Element> BlasMatrix< _Element>::transpose(BlasMatrix< _Element> & tM) const
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

	template <class _Element>
	void BlasMatrix< _Element>::transpose()
	{
		size_t r = this->rowdim() ;
		size_t c = this->coldim() ;
		if ( r == c) {
			for (size_t i = 0 ; i < r ; ++i)
				for (size_t j = i+1 ; j < c ; ++j)
					std::swap(this->refEntry(i,j),this->refEntry(j,i));
		}
		else {
			// maybe this should be possible on a DenseMatrix sharing its data with
			// a BlasMatrix but with rowdim/coldim integer members.
			// Or this has an integer rowdim/coldim member, mutable.
			throw LinBoxError("you cannot transpose a BlasMatrix in place where m != n...");
		}
	}

	template <class _Element>
	void BlasMatrix< _Element>::reverseRows()
	{
		size_t r = this->rowdim()/2 ;
		typedef UnparametricField< _Element> Field ;
		Field F ;
		VectorDomain<Field> Vd(F);
		for (size_t i = 0 ; i <  r ; ++i) {
			Vd.swap( this->rowBegin()+i,
				 this->rowBegin()+(r-1-i) );
		}

	}

	template <class _Element>
	void BlasMatrix< _Element>::reverseCols()
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

	template <class _Element>
	void BlasMatrix< _Element>::reverse()
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

	template <class _Element>
	class BlasMatrix< _Element>::ConstRowIterator {
	public:
		ConstRowIterator (const typename Rep::const_iterator& p, size_t len, size_t d) :
			_row (p, p + len), _dis (d)
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
			_row = ConstRow (_row.begin () - _dis, _row.end () - _dis);
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
			_row = ConstRow (_row.begin () + _dis, _row.end () + _dis);
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
			return ConstRowIterator (_row.begin () + _dis * i, _row.size (), _dis);
		}

		ConstRowIterator& operator += (int i)
		{
			_row = ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i);
			return *this;
		}

		ConstRow operator[] (int i) const
		{
			return ConstRow (_row.begin () + _dis * i, _row.end () + _dis * i);
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

	template <class _Element>
	class BlasMatrix< _Element>::RowIterator {
	public:
		RowIterator (const typename Rep::iterator& p, size_t len, size_t d) :
			_row (p, p + len), _dis (d)
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
			_row = Row (_row.begin () + _dis, _row.end () + _dis);
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
			_row = Row (_row.begin () - _dis, _row.end () - _dis);
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
			return RowIterator (_row.begin () + _dis * i, _row.size (), _dis);
		}

		RowIterator& operator += (int i)
		{
			_row = Row (_row.begin () + _dis * i, _row.end () + _dis * i);
			return *this;
		}

		Row operator[] (int i) const
		{
			return Row (const_cast<Row&> (_row).begin () + _dis * i,
				    const_cast<Row&> (_row).end () + _dis * i);
		}

		Row* operator-> ()
		{
			return &_row;
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

	template <class _Element>
	class BlasMatrix< _Element>::ConstColIterator {
	public:
		ConstColIterator (typename Rep::const_iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename Rep::const_iterator> (p, stride),
			      Subiterator<typename Rep::const_iterator> (p + len * stride, stride)),
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
			_col = ConstCol (Subiterator<typename Rep::const_iterator> (_col.begin ().operator-> () + 1, _stride),
					 Subiterator<typename Rep::const_iterator> (_col.end ().operator-> () + 1, _stride));
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

	template <class _Element>
	class BlasMatrix< _Element>::ColIterator {
	public:
		ColIterator (typename Rep::iterator p, size_t stride, size_t len) :
			_col (Subiterator<typename Rep::iterator> (p, stride),
			      Subiterator<typename Rep::iterator> (p + len * stride, stride)), _stride (stride)
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
			_col = Col (Subiterator<typename Rep::iterator> (_col.begin ().operator-> () + 1, _stride),
				    Subiterator<typename Rep::iterator> (_col.end ().operator-> () + 1, _stride));
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
			return &_col;
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
	template <class _Element>
	class BlasMatrix< _Element>::IndexedIterator {
		size_t _r_index;
		size_t _c_index;
		size_t _dim;
		typename Rep::iterator _begin;

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

	 _Element &operator * () const
		{
			return *(_begin + (_r_index * _dim + _c_index));
		}


	 _Element* operator -> () const
		{
			return _begin + (_r_index * _dim + _c_index);
		}


		size_t rowIndex () const
		{
			return _r_index;
		}

		size_t colIndex () const
		{
			return _c_index;
		}

		const _Element &value () const
		{
			return *(_begin + (_r_index * _dim + _c_index));
		}


	};

	template <class _Element>
	class BlasMatrix< _Element>::ConstIndexedIterator {
		size_t _r_index;
		size_t _c_index;
		size_t _dim;
		typedef _Element value_type;
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

		const _Element &operator * () const
		{
			return *(_begin + (_r_index * _dim + _c_index));
		}

		const _Element *operator -> () const
		{
			return _begin + (_r_index * _dim + _c_index);
		}

		size_t rowIndex () const
		{
			return _r_index;
		}

		size_t colIndex () const
		{
			return _c_index;
		}

		const _Element &value() const
		{
			return *(_begin + (_r_index * _dim + _c_index));
		}
	};

	/*   */

	// Entry access  view.  Size m*n vector in C (row major) order.
	template <class _Element>
	typename BlasMatrix< _Element>::Iterator BlasMatrix< _Element>::Begin ()
	{
		return _rep.begin ();
	}

	template <class _Element>
	typename BlasMatrix< _Element>::Iterator BlasMatrix< _Element>::End ()
	{
		return _rep.end ();
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstIterator BlasMatrix< _Element>::Begin () const
	{
		return _rep.begin ();
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstIterator BlasMatrix< _Element>::End () const
	{
		return _rep.end ();
	}

	/*   Indexed  */

	template <class _Element>
	typename BlasMatrix< _Element>::IndexedIterator BlasMatrix< _Element>::IndexedBegin ()
	{
		return IndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template <class _Element>
	typename BlasMatrix< _Element>::IndexedIterator BlasMatrix< _Element>::IndexedEnd ()
	{
		return IndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstIndexedIterator BlasMatrix< _Element>::IndexedBegin () const
	{
		return ConstIndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstIndexedIterator BlasMatrix< _Element>::IndexedEnd () const
	{
		return ConstIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	/*  Row  */

	template <class _Element>
	typename BlasMatrix< _Element>::RowIterator BlasMatrix< _Element>::rowBegin ()
	{
		return RowIterator (_rep.begin (), _col, _col);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::RowIterator BlasMatrix< _Element>::rowEnd ()
	{
		return RowIterator (_rep.end (), _col, _col);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstRowIterator BlasMatrix< _Element>::rowBegin () const
	{
		return ConstRowIterator (_rep.begin (), _col, _col);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstRowIterator BlasMatrix< _Element>::rowEnd () const
	{
		return ConstRowIterator (_rep.end (), _col, _col);
	}

	/*  Col */

	template <class _Element>
	typename BlasMatrix< _Element>::ColIterator BlasMatrix< _Element>::colBegin ()
	{
		return  typename BlasMatrix< _Element>::ColIterator (_rep.begin (), _col, _row);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ColIterator BlasMatrix< _Element>::colEnd ()
	{
		return  typename BlasMatrix< _Element>::ColIterator (_rep.begin ()+_col, _col, _row);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstColIterator BlasMatrix< _Element>::colBegin () const
	{
		return  typename BlasMatrix< _Element>::ConstColIterator (_rep.begin (), _col, _row);
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstColIterator BlasMatrix< _Element>::colEnd () const
	{
		return  typename BlasMatrix< _Element>::ConstColIterator (_rep.begin ()+_col, _col, _row);
	}

	/*  operators */
	template <class _Element>
	typename BlasMatrix< _Element>::Row BlasMatrix< _Element>::operator[] (size_t i)
	{
		return Row (_rep.begin () + i * _col, _rep.begin () + (i * _col +_col));
	}

	template <class _Element>
	typename BlasMatrix< _Element>::ConstRow BlasMatrix< _Element>::operator[] (size_t i) const
	{
		return Row (_rep.begin () + i * _col, _rep.begin () + ( i * _col + _col));
	}

} // LinBox

//////////////////
//     MISC     //
//////////////////

namespace LinBox
{
	template <class _Element>
	template <class Vector>
	Vector &BlasMatrix< _Element>::columnDensity (Vector &v) const
	{
		std::fill (v.begin (), v.end (), _row);
		return v;
	}

} // LinBox
#endif // __LINBOX_blas_matrix_INL
