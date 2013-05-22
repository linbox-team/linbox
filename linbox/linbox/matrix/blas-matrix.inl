/* linbox/matrix/blas-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
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
 * @file matrix/blas-matrix.inl
 * @ingroup matrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field</code>s.
 */

#ifndef __LINBOX_blas_matrix_INL
#define __LINBOX_blas_matrix_INL

/////////////////
//   PRIVATE   //
/////////////////

namespace LinBox
{
	template<class _Field>
	void BlasMatrix<_Field>::createBlasMatrix (const BlasMatrix<_Field>& A)
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

	template<class _Field>
	void BlasMatrix<_Field>::createBlasMatrix (const Element * v)
	{
		const_pointer iter_v = v ;
		const_pointer v_end = v+(_col*_row) ;
		Iterator  iter_addr = this->Begin();
		for (; v != v_end ; ++v, ++iter_addr)
			field().init(*iter_addr,*v);
	}

	template<class _Field>
	void BlasMatrix<_Field>::createBlasMatrix (const std::vector<Element> & v)
	{
		typename std::vector< Element>::const_iterator iter_value = v.begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != v.end(); ++iter_value,++iter_addr)
			field().init(*iter_addr,*iter_value);
	}


	template<class _Field>
	template <class _Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const _Matrix& A,
						   MatrixContainerCategory::BlasContainer)
	{
		linbox_check( areFieldEqual(A.field(),field()));
		typename _Matrix::ConstIterator         iter_value = A.Begin();
		Iterator  iter_addr = this->Begin();
		for (;iter_value != A.End(); ++iter_value,++iter_addr)
			field().init(*iter_addr, *iter_value);
	}

	template<class _Field>
	template <class Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const Matrix& A,
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

	template<class _Field>
	template <class Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const Matrix& A,
						   MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field()) );

		std::vector<Element> e(A.coldim(), field().zero), tmp(A.rowdim());
		ColIterator col_p;

		typename BlasMatrix< _Field>::Col::iterator elt_p;
		typename std::vector<Element>::iterator e_p, tmp_p;


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

	template<class _Field>
	template <class _Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const _Matrix& A,
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

	template<class _Field>
	template <class Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const Matrix& A,
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

	template<class _Field>
	template <class Matrix>
	void BlasMatrix<_Field>::createBlasMatrix (const Matrix& A,
						   const size_t i0,const size_t j0,
						   const size_t m, const size_t n,
						   MatrixContainerCategory::Blackbox)
	{
		linbox_check( areFieldEqual(A.field(),field() ) );


		std::vector<Element> e(A.coldim(), field().zero), tmp(A.rowdim());
		ColIterator col_p;

		typename BlasMatrix< _Field>::Col::iterator elt_p;
		typename std::vector<Element>::iterator e_p, tmp_p;


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


	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F) :
		_row(0),_col(0),_rep(0),_ptr(NULL),
		_field(&F),_MD(F),_VD(F),_use_fflas(false)
	{ }

	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix () //:
			//_row(0),_col(0),_rep(0),_ptr(NULL),
			//_field(Field()),_MD(_field ),_VD(_field )
		{}

	template <class _Field>
	void BlasMatrix< _Field>::init(const _Field &F, size_t r, size_t c)
	{
		_field = &F; _row = r; _col = c;
		 _rep.resize(r*c, F.zero);
		_ptr = &_rep[0];
		_VD.init(F); _MD.init(F);
	}

	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix ( const _Field &F, const uint32_t& m, const T& n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}

	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F, const int64_t& m, const T& n) :
		_row((size_t)m),_col((size_t)n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		linbox_check(n>=0);
		linbox_check(m>=0);
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}

#ifdef __GNUC__
#ifndef __x86_64__
#if (__GNUC__ == 4 && __GNUC_MINOR__ ==4 && __GNUC_PATCHLEVEL__==5)
	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F, const long & m, const T& n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		linbox_check(n>=0);
		linbox_check(m>=0);
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}
#endif
#endif
#endif

#if defined(__APPLE__) || (defined(__s390__) && !defined(__s390x__))
	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F, const unsigned long & m, const T& n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		linbox_check(n>=0);
		linbox_check(m>=0);
		_use_fflas = Protected::checkBlasApply(field(),_col);
	}
#endif


	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F, const uint64_t &m, const T & n) :
		_row(m),_col((size_t)n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		//!@todo
		// linbox_check_non_neg(n);
		// linbox_check(n>=0);
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}

	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F, const int32_t & m, const T &n) :
		_row((size_t) m),_col((size_t)n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		linbox_check(isPositive<T>(n));
		linbox_check(m>=0);
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}

	template <class _Field>
	template<class T>
	BlasMatrix< _Field>::BlasMatrix ( const _Field &F, const Integer & m, const T &n) :
		_row(m),_col(n),_rep(_row*_col, F.zero),_ptr(&_rep[0]),
		_field(&F),_MD(F),_VD(F)
	{
		//!@todo check m,n not too big ?
		linbox_check(n>=0);
		linbox_check(m>=0);
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}



	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix(MatrixStream<_Field>& ms) :
		_row(0),_col(0),_rep(0),
		_field(&(ms.getField())),_MD(field() ),_VD(field() )
	{
		if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
			throw ms.reportError(__FUNCTION__,__LINE__);
		_ptr = &_rep[0];
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}

	template <class _Field>
	template <class StreamVector>
	BlasMatrix< _Field>::BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) :
		_row(stream.size ()), _col(stream.dim ()), _rep(_row*_col), _ptr(&_rep[0]),
		_field (&F), _MD (F), _VD(F)
	{
		StreamVector tmp(F);
		typename BlasMatrix<Field>::RowIterator p;

		VectorWrapper::ensureDim (tmp, stream.dim ());

		for (p = BlasMatrix<Field>::rowBegin (); p != BlasMatrix<Field>::rowEnd (); ++p) {
			stream >> tmp;
			_VD.copy (*p, tmp);
		}
		_use_fflas = Protected::checkBlasApply(field(), _col);
	}


	template <class _Field>
	template <class Matrix>
	BlasMatrix< _Field>::BlasMatrix (const Matrix &A) :
		_row(A.rowdim()),_col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )

	{
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(A, typename MatrixContainerTrait<Matrix>::Type());
	}

	template <class _Field>
	template <class Matrix>
	BlasMatrix< _Field>::BlasMatrix (const Matrix& A,
					 const size_t i0, const size_t j0,
					 const size_t m, const size_t n) :
		_row(m),_col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )
	{
		_use_fflas = Protected::checkBlasApply(field(), _col);
		// makePointer();
		createBlasMatrix(A, i0, j0, m, n,
				 typename MatrixContainerTrait<Matrix>::Type());
	}

	template <class _Field>
	template<class _Matrix>
	BlasMatrix< _Field>::BlasMatrix (const _Matrix &A,  const _Field &F) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F),_MD(field() ),_VD(field() )
	{
		_use_fflas = Protected::checkBlasApply(field(), _col);
		typename _Matrix::template rebind<_Field>()(*this,A);
	}

	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix (const BlasMatrix< _Field>& A) :
		_row(A.rowdim()), _col(A.coldim()),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&(A.field())),_MD(field() ),_VD(field() )
	{
		_use_fflas = Protected::checkBlasApply(field(), _col);
		// makePointer();
		createBlasMatrix(A);
	}

	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F,
					 const std::vector<typename _Field::Element>& v,
					 size_t m, size_t n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F),_MD(field() ),_VD(field() )
	{
		linbox_check(v.size() == m*n);
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(v);
	}

	template <class _Field>
	BlasMatrix< _Field>::BlasMatrix (const _Field &F,
					 const typename _Field::Element * v,
					 size_t m, size_t n) :
		_row(m), _col(n),_rep(_row*_col),_ptr(&_rep[0]),
		_field(&F), _MD(field() ),_VD(field() )
	{
		// makePointer();
		_use_fflas = Protected::checkBlasApply(field(), _col);
		createBlasMatrix(field(),v);
	}

	template <class _Field>
	BlasMatrix< _Field>::~BlasMatrix ()
	{
		// if (_ptr)
		// free(_ptr);
	}

	template <class _Field>
	void BlasMatrix< _Field>::random()
	{
		typename _Field::Element x; field().init(x);
		typename _Field::RandIter r(field());
		for (size_t i = 0; i < rowdim(); ++i)
			for (size_t j = 0; j < coldim(); ++j)
				setEntry(i, j, r.random(x));
	}
} // LinBox

///////////////////
//      I/O      //
///////////////////

namespace LinBox
{

	template <class _Field>
	std::istream &BlasMatrix< _Field>::read (std::istream &file)
	{
		MatrixStream<Field> ms(field(), file);
		if( !ms.getArray(_rep) || !ms.getDimensions(_row, _col) )
			throw ms.reportError(__FUNCTION__,__LINE__);
		_ptr = &_rep[0];
		_use_fflas = Protected::checkBlasApply(field(), _col);
		return file;
	}

#if 0 /* now writing is done in BlasSubmatrix. */

	template <class _Field>
	std::ostream& BlasMatrix< _Field>::write (std::ostream &os) {
		// ....
		return os;
	}

	template <class _Field>
	std::ostream& BlasMatrix< _Field>::write (std::ostream &os,
						  enum LinBoxTag::Format f) const
	{

		ConstRowIterator p;

		switch(f) {
		case (LinBoxTag::FormatPlain) : /*  raw output */
			{
				integer c;
				int wid;

				field().cardinality (c);
				if (c >0) {
					wid = (int) ceil (log ((double) c) / M_LN10);
				}
				else {
					integer tmp;
					size_t max=0;
					ConstIterator it = Begin();
					for (; it != End(); ++it){
						field().convert(tmp,*it);
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
                                                field().write (os, *pe);
						os << " ";
					}

					os << "]" << std::endl;
				}
			}
			break;
		case (LinBoxTag::FormatMaple) : /*  maple format */
			{
				os << "Matrix( " << rowdim() << ',' << coldim() << ",[" << std::endl;
				for (p = rowBegin (); p != rowEnd (); ) {
					typename ConstRow::const_iterator pe;

					os << " [ ";

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
		case (LinBoxTag::FormatHTML) : /*  html format */
			{
				os << "<table border=\"1\">" ;
				for (p = rowBegin (); p != rowEnd (); ) {
					typename ConstRow::const_iterator pe;

					os << "<tr>";

					for (pe = p->begin (); pe != p->end (); ) {
						field().write (os << "<td>", *pe)<<"</td>";
						++pe ;
					}

					os << "</tr>" << std::endl;
					++p ;
				}
				os << "</table>" ;
			}
			break;
		case (LinBoxTag::FormatLaTeX) : /*  latex format (pmatrix) */
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
#endif

	template <class _Field>
	BlasMatrix< _Field>& BlasMatrix< _Field>::operator= (const BlasMatrix< _Field>& A)
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
	template <class _Field>
	template<typename _Tp1>
	struct BlasMatrix< _Field>::rebind {
		typedef BlasMatrix<_Tp1> other;

		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
			// typedef typename BlasMatrix< _Field>::ConstIndexedIterator ConstIndexedIterator ;
			// typedef typename BlasMatrix< _Field>::ConstIterator ConstIterator ;
			ConstIterator         iter_value = A.Begin();
			ConstIndexedIterator  iter_index = A.IndexedBegin();
			typename _Tp1::Element tmp;
			for (;iter_value != A.End(); ++iter_value,++iter_index){
				Ap.field().init(  tmp, *iter_value );
				Ap.setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
			}
		}
		};
#endif

#if 1 /*  HOM */
	template <class _Field>
	template<typename _Tp1>
	struct BlasMatrix< _Field>::rebind {
		typedef BlasMatrix<_Tp1> other;

		void operator() (other & Ap, const Self_t& A) {
			typedef typename BlasMatrix<_Field>::ConstIterator ConstSelfIterator ;
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

	template <class _Field>
	size_t BlasMatrix< _Field>::rowdim() const
	{
		return _row ;
	}

	template <class _Field>
	size_t BlasMatrix< _Field>::coldim() const
	{
		return _col ;
	}

	template <class _Field>
	size_t  BlasMatrix< _Field>::getStride() const
	{
		return _col;
	}

	template <class _Field>
	size_t&  BlasMatrix< _Field>::getWriteStride()
	{
		return _col;
	}

	template <class _Field>
	void BlasMatrix< _Field>::resize (size_t m, size_t n, const Element& val )
	{
#ifndef NDEBUG
		if (_col > 0 && _col != n)
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


	template <class _Field>
	typename BlasMatrix< _Field>::pointer
	BlasMatrix< _Field>::getPointer() const
	{
		return _ptr;
	}

	template <class _Field>
	typename BlasMatrix< _Field>::const_pointer &
	BlasMatrix< _Field>::getConstPointer() const
	{
		return (const_pointer)_ptr;
	}


	template <class _Field>
	typename BlasMatrix< _Field>::pointer&
	BlasMatrix< _Field>::getWritePointer()
	{
		return _ptr;
	}

	template <class _Field>
	void BlasMatrix< _Field>::setEntry (size_t i, size_t j, const Element &a_ij)
	{
		_ptr[i * _col + j] = a_ij;
	}

	template <class _Field>
	typename _Field::Element & BlasMatrix< _Field>::refEntry (size_t i, size_t j)
	{
		return _ptr[i * _col + j];
	}

	template <class _Field>
	const typename _Field::Element & BlasMatrix< _Field>::getEntry (size_t i, size_t j) const
	{
		return _ptr[i * _col + j];
	}

	template <class _Field>
	typename _Field::Element & BlasMatrix< _Field>::getEntry (Element &x, size_t i, size_t j) const
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
	template <class _Field>
	BlasMatrix< _Field> BlasMatrix< _Field>::transpose(BlasMatrix< _Field> & tM) const
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
			size_t start, next, i;
			T tmp;

			for (start = 0; start <= w * h - 1; ++start) {
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

	template <class _Field>
	template<bool _IP>
	void BlasMatrix< _Field>::transpose()
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
				BlasMatrix<_Field> MM(*this);
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

	template <class _Field>
	void BlasMatrix< _Field>::transpose()
	{
		this->transpose<false>();
	}

	template <class _Field>
	void BlasMatrix< _Field>::reverseRows()
	{
		size_t r = this->rowdim()/2 ;
		for (size_t i = 0 ; i <  r ; ++i) {
			_VD.swap( this->rowBegin()+i,
				 this->rowBegin()+(r-1-i) );
		}

	}

	template <class _Field>
	void BlasMatrix< _Field>::reverseCols()
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

	template <class _Field>
	void BlasMatrix< _Field>::reverse()
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

	template <class _Field>
	class BlasMatrix< _Field>::ConstRowIterator {
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

	template <class _Field>
	class BlasMatrix< _Field>::RowIterator {
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

	template <class _Field>
	class BlasMatrix< _Field>::ConstColIterator {
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

	template <class _Field>
	class BlasMatrix< _Field>::ColIterator {
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
	template <class _Field>
	class BlasMatrix< _Field>::IndexedIterator {
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

	template <class _Field>
	class BlasMatrix< _Field>::ConstIndexedIterator {
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
	template <class _Field>
	typename BlasMatrix< _Field>::Iterator BlasMatrix< _Field>::Begin ()
	{
		return _rep.begin ();
	}

	template <class _Field>
	typename BlasMatrix< _Field>::Iterator BlasMatrix< _Field>::End ()
	{
		return _rep.end ();
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstIterator BlasMatrix< _Field>::Begin () const
	{
		return _rep.begin ();
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstIterator BlasMatrix< _Field>::End () const
	{
		return _rep.end ();
	}

	/*   Indexed  */

	template <class _Field>
	typename BlasMatrix< _Field>::IndexedIterator BlasMatrix< _Field>::IndexedBegin ()
	{
		return IndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template <class _Field>
	typename BlasMatrix< _Field>::IndexedIterator BlasMatrix< _Field>::IndexedEnd ()
	{
		return IndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstIndexedIterator BlasMatrix< _Field>::IndexedBegin () const
	{
		return ConstIndexedIterator (coldim (), 0, 0, _rep.begin ());
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstIndexedIterator BlasMatrix< _Field>::IndexedEnd () const
	{
		return ConstIndexedIterator (coldim (), rowdim (), 0, _rep.begin ());
	}

	/*  Row  */

	template <class _Field>
	typename BlasMatrix< _Field>::RowIterator BlasMatrix< _Field>::rowBegin ()
	{
		return RowIterator (_rep.begin (), _col, _col);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::RowIterator BlasMatrix< _Field>::rowEnd ()
	{
		return RowIterator (_rep.end (), _col, _col);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstRowIterator BlasMatrix< _Field>::rowBegin () const
	{
		return ConstRowIterator (_rep.begin (), _col, _col);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstRowIterator BlasMatrix< _Field>::rowEnd () const
	{
		return ConstRowIterator (_rep.end (), _col, _col);
	}

	/*  Col */

	template <class _Field>
	typename BlasMatrix< _Field>::ColIterator BlasMatrix< _Field>::colBegin ()
	{
		return  typename BlasMatrix< _Field>::ColIterator (_rep.begin (), _col, _row);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ColIterator BlasMatrix< _Field>::colEnd ()
	{
		return  typename BlasMatrix< _Field>::ColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstColIterator BlasMatrix< _Field>::colBegin () const
	{
		return  typename BlasMatrix< _Field>::ConstColIterator (_rep.begin (), _col, _row);
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstColIterator BlasMatrix< _Field>::colEnd () const
	{
		return  typename BlasMatrix< _Field>::ConstColIterator (_rep.begin ()+(ptrdiff_t)_col, _col, _row);
	}

	/*  operators */
	template <class _Field>
	typename BlasMatrix< _Field>::Row BlasMatrix< _Field>::operator[] (size_t i)
	{
		return Row (_rep.begin () +(ptrdiff_t)( i * _col), _rep.begin () + (ptrdiff_t)(i * _col +_col));
	}

	template <class _Field>
	typename BlasMatrix< _Field>::ConstRow BlasMatrix< _Field>::operator[] (size_t i) const
	{
		return Row (_rep.begin () +(ptrdiff_t) (i * _col), _rep.begin () + (ptrdiff_t)( i * _col + _col));
	}

} // LinBox

//////////////////
//     MISC     //
//////////////////

namespace LinBox
{
	template <class _Field>
	template <class Vector>
	Vector& BlasMatrix< _Field>::columnDensity (Vector &v) const
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
	template <class _Field>
	template <class Vector1, class Vector2>
	Vector1&  BlasMatrix< _Field>::apply (Vector1& y, const Vector2& x) const
	{
		//_stride ?
		if (_use_fflas){
			//!@bug this supposes &x[0]++ == &x[1]
                        // PG: try to discover stride of x and y (not use it works on every platform)
                        size_t ldx,ldy;
                        ldx=(size_t)(&x[1] - &x[0]);
                        ldy=(size_t)(&y[1] - &y[0]);

			FFLAS::fgemv((typename Field::Father_t) field(), FFLAS::FflasNoTrans,
				      _row, _col,
				      field().one,
				      _ptr, getStride(),
				      &x[0],ldx,
				      field().zero,
				      &y[0],ldy);
		}
		else {
			_MD. vectorMul (y, *this, x);
#if 0
			typename BlasMatrix<_Field>::ConstRowIterator i = this->rowBegin ();
			typename Vector1::iterator j = y.begin ();

			for (; j != y.end (); ++j, ++i)
				_VD.dot (*j, *i, x);
#endif
		}
		return y;
	}
	template <class _Field>
	BlasVector<_Field>&  BlasMatrix< _Field>::apply (BlasVector<_Field>& y, const BlasVector<_Field>& x) const
	{
		//_stride ?
		if (_use_fflas){
			//!@bug this supposes &x[0]++ == &x[1]
                        // PG: try to discover stride of x and y (not use it works on every platform)
                        size_t ldx,ldy;
                        ldx=(size_t)x.getStride();
                        ldy=(size_t)y.getStride();

			FFLAS::fgemv((typename Field::Father_t) field(), FFLAS::FflasNoTrans,
				      _row, _col,
				      field().one,
				      _ptr, getStride(),
				      x.getPointer(),x.getStride(),
				      field().zero,
				      y.getWritePointer(),y.getStride());
		}
		else {
			_MD. vectorMul (y, *this, x);
		}
		return y;
	}

	template <class _Field>
	template <class Vector1, class Vector2>
	Vector1&  BlasMatrix< _Field>::applyTranspose (Vector1& y, const Vector2& x) const
	{

		//_stride ?
		if (_use_fflas) {
                        // PG: try to discover stride of x and y (not use it works on every platform)
                        size_t ldx,ldy;
                        ldx=(size_t)(&x[1] - &x[0]);
                        ldy=(size_t)(&y[1] - &y[0]);

			FFLAS::fgemv((typename Field::Father_t) field(), FFLAS::FflasTrans,
				      _row, _col,
				      field().one,
				      _ptr, getStride(),
				      &x[0],ldx,
				      field().zero,
				      &y[0],ldy);
		}
		else {
			typename BlasMatrix<_Field>::ConstColIterator i = this->colBegin ();
			typename Vector1::iterator j = y.begin ();
			for (; j != y.end (); ++j, ++i)
				_VD.dot (*j, x, *i);
		}

		return y;
	}

	template <class _Field>
	const _Field& BlasMatrix< _Field>::field() const
	{
		return *_field;
	}

#if 0 /* why not ? */
	template <class _Field>
	_Field& BlasMatrix< _Field>::field()
	{
		return const_cast<_Field&>(_field );
	}
#endif
}


#endif // __LINBOX_blas_matrix_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

