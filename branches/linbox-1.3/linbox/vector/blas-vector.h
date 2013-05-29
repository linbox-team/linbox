/* linbox/matrix/blas-vector.h
 * Copyright (C) 2013 the LinBox group
 *
 * Written by :
 * BB <bbboyer@ncsu.edu>
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


/*! @file vector/blas-vector.h
 * @ingroup vector
 * A \c BlasVector<\c _Field > represents a vector as an array of
 * <code>_Field::Element</code>s.
 *
 */

#ifndef __LINBOX_vector_blas_vector_H
#define __LINBOX_vector_blas_vector_H

#include "linbox/util/debug.h"
#include "linbox/algorithms/linbox-tags.h"

namespace LinBox { /* BlasVector */

	template<class _Field>
	class BlasSubvector ;

	template<class _Field>
	class BlasVector {

	public:
		typedef _Field                                Field;
		typedef typename Field::Element             Element;    //!< Element type
		typedef typename RawVector<Element>::Dense      Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef typename Rep::pointer               pointer;    //!< pointer type to elements
		typedef const pointer                 const_pointer;    //!< const pointer type
		typedef BlasVector<_Field>                   Self_t;    //!< Self type
		typedef BlasSubvector<_Field>         subVectorType;    //!< Submatrix type
		typedef BlasVector<_Field>               vectorType;    //!< matrix type
                typedef BlasVector<_Field>                 blasType;    //!< blas type
	protected:
		size_t			       _size;
		size_t                       _stride;
		Rep			        _rep;
		pointer			        _ptr;
		const Field		    * _field;
	private:
		void createBlasVector(const BlasVector<_Field> & V)
		{
			//! use std::copy somehow ?
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,V.getEntry(i));
		}

		void createBlasVector(const BlasVector<_Field> & V
				      , const size_t i0, const size_t str)
		{
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,V.getEntry(i0+i*str));
		}

		void createBlasVector(const BlasSubvector<_Field> & V)
		{
			//! use std::copy somehow ?
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,V.getEntry(i));
		}

		void createBlasVector ( const Element * V)
		{
			//! use std::copy somehow ?
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,V[i]);
		}

		void createBlasVector ( const Element * V
					, const size_t i0, const size_t str)
		{
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,V[i0+i*str]);
		}


		void createBlasVector ( const std::vector<Element> & V)
		{
			createBlasVector(&V[0]);
		}

		void createBlasVector ( const std::vector<Element> & V
					, const size_t i0, const size_t str)
		{
			createBlasVector(&V[0],i0,str);
		}

		/** Gets a vector from data in a BlasMatrix.
		 * Vector starts at \c (i0,j0) and next element is \p stride away.
		 * For instance, row \c i correpsonds to args \c (i,0,1) and col \c j is \c (0,j,lda)
		 * It is possible to get vectors that extend over multiple rows/columns
		 * and \p stride need not be equal to \c 1 or \c lda.
		 */
		void createBlasVector ( const BlasMatrix<Field> & A, size_t i0, size_t j0, size_t str)
		{
			if ((str == 1) && (i0+_size*str<A.coldim())){
				for (size_t i = 0 ; i < _size ; ++i)
					setEntry(i,A.getEntry(i0+i,j0));
			}
			else if ((str == A.coldim()) && j0+_size*str<A.rowdim()) {
				for (size_t i = 0 ; i < _size ; ++i)
					setEntry(i,A.getEntry(i0,j0+i));
			}
			else{ /* is this always faster ? */
				Element * Aptr = A.getPointer(i0,j0);
				for (size_t i = 0 ; i < _size ; ++i)
					setEntry(i,Aptr[i*str]);
			}
		}

		void createBlasVector ( const BlasSubmatrix<Field> & A, size_t i0, size_t j0, size_t str)
		{
			Element * Aptr = A.getPointer(i0,j0);
			for (size_t i = 0 ; i < _size ; ++i)
				setEntry(i,Aptr[i*str]);
		}


	public:
		BlasVector (const _Field &F)  :
			_size(0),_stride(1),_rep(0),_ptr(NULL),
			_field(&F)
		{}

		BlasVector () {} ;

		void init(const _Field & F, size_t n = 0)
		{
			_field = &F;
			_size = n;
		       	_stride=1 ;
			_rep.resize(n, F.zero);
			_ptr = &_rep[0];
		}

#ifdef __GNUC__
#ifndef __x86_64__
#if (__GNUC__ == 4 && __GNUC_MINOR__ ==4 && __GNUC_PATCHLEVEL__==5)
		BlasVector (const _Field &F, const long &m) :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}
#endif
#endif
#endif

#if defined(__APPLE__) || (defined(__s390__) && !defined(__s390x__))
		BlasVector (const _Field &F, const unsigned long &m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}

#endif

		BlasVector (const _Field &F, const uint64_t &m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}


		BlasVector (const _Field &F, const int64_t &m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}


		BlasVector (const _Field &F, const uint32_t &m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}

		BlasVector (const _Field &F, const int32_t &m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}

		BlasVector (const _Field &F, const Integer & m)  :
			_size((size_t)m),_stride(1),_rep((size_t)_size, F.zero),_ptr(&_rep[0]),_field(&F)
		{}

		BlasVector (const BlasVector<_Field> &V)  :
			_size(V.size()),_stride(1),_rep(V.size(), V.field().zero),_ptr(&_rep[0]),_field(&(V.field()))
		{
			createBlasVector(V);
		}

		BlasVector (const BlasSubvector<_Field> &V)  :
			_size(V.size()),_stride(1),_rep(V.size(), V.field().zero),_ptr(&_rep[0]),_field(&(V.field()))
		{
			createBlasVector(V);
		}

		BlasVector (const BlasMatrix<_Field> &A, size_t k, enum LinBoxTag::Direction f )  :
			_size((f==LinBoxTag::Row)?(A.rowdim()):(A.coldim())),_stride(1),_rep(_size, A.field().zero),_ptr(&_rep[0]),_field(&(A.field()))
			{
				if (f==LinBoxTag::Row)
					createBlasVector(A,k,0,1);
				else // LinBoxTag::Col
					createBlasVector(A,0,k,A.coldim());

			}

		BlasVector (const BlasSubmatrix<_Field> &A, size_t k, enum LinBoxTag::Direction f )  :
			_size((f==LinBoxTag::Row)?(A.rowdim()):(A.coldim())),_stride(1),_rep(_size, A.field().zero),_ptr(&_rep[0]),_field(&(A.field()))
			{
				if (f==LinBoxTag::Row)
					createBlasVector(A,k,0,1);
				else // LinBoxTag::Col
					createBlasVector(A,0,k,A.stride());
			}

		BlasVector (const BlasMatrix<_Field> &A, size_t n, size_t i0, size_t j0, size_t str )  :
			_size(n),_stride(1),_rep(_size, A.field().zero),_ptr(&_rep[0]),_field(&(A.field()))
		{
			createBlasVector(A,i0,j0,str);
		}

		// ~BlasVector () ;

		BlasVector<_Field>& operator= (const BlasVector<_Field>& V)
		{
			if ( &V == this)
				return *this;

			_size = V.size();
			_stride = 1;
			_rep = Rep(_size);
			_ptr = &_rep[0] ;
			createBlasVector(V);

			return *this;
		}

		// template<typename _Tp1>
		// struct rebind ;

		// write

		// read


		size_t size() const
		{
			return _size;
		}

		//!should never ever be used
		size_t getStride() const
		{
			return _stride;
		}

		size_t stride() const { return getStride() ;}

		//!should never ever be used
		size_t& getWriteStride()
		{
			return _stride;
		}

		void resize (size_t n, const Element& val = Element())
		{
#ifndef NDEBUG
			if (_size < n)
				std::cerr << " ***Warning*** you are resizing a matrix, possibly loosing data. " << std::endl;
#endif
			_size = n;
			_rep.resize (n, val);
		}

		Rep & refRep() { return _rep ; }

		pointer getPointer() const
		{
			return _ptr;
		}

		const_pointer &getConstPointer() const
		{
			return (const_pointer)_ptr;
		}

		pointer& getWritePointer()  //! @bug should be called refPointer()
		{
			return _ptr;
		}

		void setEntry (size_t i, const Element &a_i)
		{
			_ptr[i] = a_i;
		}

		Element &refEntry (size_t i)
		{
			return _ptr[i];
		}

		const Element &getEntry (size_t i) const
		{
			return _ptr[i];
		}

		Element &getEntry (Element &x, size_t i) const
		{
			x = _ptr[i];
			return x;
		}

		void random()
		{
			typename _Field::Element x; field().init(x);
			typename _Field::RandIter r(field());
			for (size_t i = 0; i < size(); ++i)
				setEntry(i, r.random(x));
		}


		const _Field& field() const { return _field ;}
		_Field & field() { return _field; }


	};// BlasVector

	template <class _Field>
	class BlasSubvector {
	public :
		typedef _Field                                  Field;
		typedef typename Field::Element                 Element;      //!< Element type
		typedef BlasSubvector<_Field>                   Self_t;       //!< Self type
		typedef typename RawVector<Element>::Dense      Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef typename Rep::pointer                   pointer;    //!< pointer type to elements
		typedef const pointer                           const_pointer;    //!< const pointer type
		typedef Self_t                                  subVectorType;    //!< Subvector type
		typedef BlasVector<_Field>                      vectorType;    //!< vector type
                typedef BlasVector<_Field>                      blasType;    //!< blas type


	protected:
		Rep &_Vec;                     //!< Parent raw vector
		size_t _size;                   //!< size of Subvector
		size_t _i0;                    //!< beginning of Subvector in \p _Vec
		size_t _stride ;               //!< number of columns in \p _Vec (or stride of \p _Vec)
		Field & _field;

	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*  constructors */

		/** NULL constructor.  */
		// BlasSubvector () ;

		BlasSubvector (const BlasVector<_Field> &V,
			       size_t ibeg,
			       size_t siz,
			       size_t str
			      ) :
			_Vec (const_cast<Rep&>(V.refRep())),
			_size(siz),_i0 (ibeg),_stride(str)
			,_field(V.field())
		{}


		BlasSubvector (const BlasVector<_Field> &V) :
			_Vec (const_cast<Rep&>(V.refRep())),
			_size(V.size()),_i0 (0),_stride(1)
			,_field(V.field())
		{}



		BlasSubvector (const BlasSubvector<_Field> &SV,
			       size_t ibeg,
			       size_t siz
			       // , size_t stride //!stride in new subvector as if SV were of stride 1
			      ) :
			_Vec (SV._Vec),
			_size(siz),_i0 (SV._i0*SV._stride+ibeg)
			,_stride(SV._stride)
			// ,_stride(SV._stride*stride)
			,_field(SV.field())
		{}


		/** Copy constructor.
		 * @param SM Subvector to copy
		 */
		BlasSubvector (const BlasSubvector<_Field> &SV) :
			_Vec (SV._Vec),
			_size(SV._size),_i0 (SV._i0)
			,_stride(SV._stride)
			,_field(SV.field())
		{}

		BlasSubvector (const BlasMatrix<_Field> &M
			       , size_t beg
			       , enum LinBoxTag::Direction f ) :
			_Vec (const_cast<Rep&>(M.refRep()))
			,_size((f==LinBoxTag::Row)?(M.coldim()):(M.rowdim()))
			,_i0 ((f==LinBoxTag::Row)?(beg*M.coldim()):(beg))
			,_stride((f==LinBoxTag::Row)?(1):(M.coldim()))
			,_field(M.field())
			{}

		BlasSubvector (const BlasSubmatrix<_Field> &M
			       , size_t beg
			       , enum LinBoxTag::Direction f ) :
			_Vec (const_cast<Rep&>(M.refRep()))
			,_size((f==LinBoxTag::Row)?(M.coldim()):(M.rowdim()))
			,_i0 ((f==LinBoxTag::Row)?(M.offset()+beg*M.stride()):(M.offset()+beg))
			,_stride((f==LinBoxTag::Row)?(1):(M.stride()))
			,_field(M.field())
			{}


		//! @todo more general subvectors

		/*  Members  */

		BlasSubvector &operator = (const BlasSubvector<_Field> &SV)
		{
			if ( &SV == this)
				return *this ;
			_Vec=SV._Vec ; //!@todo use functions
			_size = SV.size();
			_i0=SV._i0;
			_stride = SV.stride();
			_field = SV.field();
			return *this ;

		}

		// template<typename _Tp1>
		// struct rebind ;

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		size_t size() const {return _size;}

		size_t getStride() const {return _stride ; }
		size_t stride() const { return getStride() ;}


		///////////////////
		//      I/O      //
		///////////////////

		// template<class Field>
		// std::istream& read (std::istream &file/*, const Field& field*/);


		// std::ostream &write (std::ostream &os) const;


		//////////////////
		//   ELEMENTS   //
		//////////////////

		pointer getPointer() const { return &(_Vec[_i0]); }

		const_pointer &getConstPointer() const { return &(_Vec[_i0]); }


		pointer& getWritePointer() { return &(_Vec[_i0]); }


		void setEntry (size_t i, const Element &a_i)
		{
			_Vec[_i0+i*_stride] = a_i;
		}

		Element &refEntry (size_t i)
		{
			return _Vec[_i0+i*_stride] ;
		}

		const Element &getEntry (size_t i) const
		{
			return _Vec[_i0+i*_stride] ;
		}

		Element &getEntry (Element &x, size_t i) const
		{
			x = _Vec[_i0+i*_stride] ;
			return x;
		}


		///////////////////
		//   ITERATORS   //
		///////////////////

		// class Iterator  ;
		// class ConstIterator ;

		// class IndexedIterator ;
		// class ConstIndexedIterator ;


		const _Field& field() const { return _field ;}
		_Field & field() { return _field; }
	};

} // LinBox

#include "blas-vector.inl"

#endif // __LINBOX_vector_blas_vector_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
