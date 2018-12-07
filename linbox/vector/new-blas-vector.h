/* linbox/matrix/blas-vector.h
 * Copyright (C) 2013 the LinBox group
 *               2018 revamped by Pascal Giorgi 
 *
 * Written by :
 * Brice Boyer (briceboyer) <boyer.brice@gmail.com>
 *               Pascal Giorgi  pascal.giorgi@lirmm.fr
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

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-tags.h"
#include "linbox/field/hom.h"
#include "linbox/vector/vector.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/subvector.h"
#include "linbox/vector/subiterator.h"

#include "fflas-ffpack/fflas/fflas.h"
namespace LinBox { /* BlasVector */


	template<class _Field, class _Storage>
	class BlasVector : public Subvector<Subiterator<typename _Storage::iterator>> {

	public:
		typedef _Field                                   Field;
		typedef typename Field::Element                Element;    //!< Element type
        typedef typename Field::Element_ptr        Element_ptr;    //!< Element type
		typedef _Storage                               Storage;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef BlasVector<_Field,_Storage>             Self_t;    //!< Self type
		typedef Self_t                              vectorType;    //!< vector type
		typedef BlasSubvector<Self_t>            subVectorType;    //!< SubVector type
        typedef BlasSubvector<const Self_t> constSubVectorType;    //!< const SubVector type



        /* iterators */
        /* all iterator method are inherited from Subvector*/
        /* Danger: all inherited iterators must be updated after _rep is modified. This is achieved with the method updateIterators() */
		typedef Subvector<Subiterator<typename _Storage::iterator>> Father_t;
		typedef typename Father_t::iterator                         iterator;
		typedef typename Father_t::const_iterator             const_iterator;


	protected:
		size_t			       _size;
		Storage	    	        _rep;
		Element_ptr			    _ptr;
		const Field		     &_field;

	private:
        // adjust the iterators in Father_t
		void updateStorage()
		{
            _size= _rep.size();
            _ptr = _rep.data();
			Father_t::_begin = iterator (_rep.begin() , 1);
			Father_t::_end   = iterator (_rep.begin()+(ptrdiff_t)_size , 1);
		}

	public:

        /** Resize the vector to the given dimension.
		 * @param n vector dimension
         * used for allocating every BlasVector
         */
		void resize (size_t n) {
            _rep.resize(n);
            update_storage();
            for (auto& it:_rep)
                field().init(it);
        }


        
        //////////////////
		// CONSTRUCTORS //
		//////////////////

        // PG: not sure we need it somewhere
		//BlasVector () {} ;

        /*! Copy Constructor of a vector (copying data).
		 * @param V vector to be copied.
		 */
		BlasVector (const Self_t &V) : Father_t(), _field(V.field())
        {
            resize(V.size());
            for (auto it=_rep.begin(), jt=V.begin();it!=_rep.end();it++,jt++)
                field().assign(*it,*jt);
        }


		/*! Allocates a new \f$ 0 \f$ vector (shaped and ready).*/
		BlasVector (const _Field &F): Father_t(), _field(V.field())
        {
            resize(0);
        }

        // This constructor is templated because if SizeType was a size_t, then calling it with a signed const litteral 
        // (e.g. v(F,3)) would fail to launch this constructor, but rather go in the templated one (_Field, VectorBase) 
        // on OSX where long can not be cast to size_t).
        /*! Allocates a new \f$ m \f$ zero vector (shaped and ready).
		 * @param F
		 * @param m vector dimension
		 */
        template<class SizeType, typename std::enable_if<std::is_arithmetic<SizeType>::value, int>::type=0>
		BlasVector (const _Field &F, const SizeType &m): Father_t(), _field(F)
        {
            resize(m);
        }

        /*! Create a BlasMatrix from an iterator of elements
		 * @param F Field of the created vector
		 * @param it iterator to be copied from
		 * @param m vector dimension
		 * @param inc increment value for iterating
		 */
        template<typename ConstIterator>
		BlasVector(const _Field & F, const ConstIterator& jt, const size_t m) : Father_t(), _field(F)
        {
            resize(m);
            for (auto it=_rep.begin();it!=_rep.end();it++,jt++)
                field().assign(*it,*jt);        
        }

        
        /*! Create a BlasVector from another vector defined over a different field (use homomorphism if it exists)
		 * @param F Field of the created vector
		 * @param V Vector to be copied
		 */
		template<class OtherVector>//, typename std::enable_if<!std::is_arithmetic<VectorBase>::value, int>::type=0>
		BlasVector (const OtherVector & V, const _Field & F) : Father_t(), _field(F)
        {
            resize(V.size());
            typename OtherVector::template rebind<_Field>()(*this,V);        
        }
        
        /// Destructor.
		~BlasVector () {}

        //! operator = (copying data)
		Self_t& operator= (const Self_t& A) {
            if (this!=&A){
                _field=A.field();
                resize(A.size());
                FFLAS::fassign(field(),_size,_ptr,1,A._ptr,1);                               
            }
            return *this;
        }
        
        //! operator = (copying data from different vector type)
        template<class _Vector>
        Self_t& operator= (const _Vector& A){
            if (this!=&A){
                _field=A.field();
                resize(A.size());
                for (auto it=_rep.begin(), jt=A.begin();it!=_rep.end();it++,jt++)
                    field().assign(*it,*jt);                        
            }
            return *this;
        }
        
		//! Rebind operator
        template<typename _Tp1>
		struct rebind {
			typedef BlasVector<_Tp1> other;
            
			void operator() (other & Ap, const Self_t& A) {
				typedef typename Self_t::ConstIterator ConstSelfIterator ;
				typedef typename other::Iterator OtherIterator ;
				OtherIterator    Ap_i = Ap.Begin();
				ConstSelfIterator A_i = A.Begin();
				Hom<Field, _Tp1> hom(A. field(), Ap. field()) ;
				for ( ; A_i != A. End(); ++ A_i, ++ Ap_i)
					hom.image (*Ap_i, *A_i);
			}
		};


		/////////////////
		//  ACCESSORS  //
		/////////////////

        const _Field& field() const { return _field;}
        
        // dimension of the vector
        size_t size() const{ return _size; }

        /*!_@internal
		 * Get read-only access to the vector  data.
		 */
		Element_ptr getPointer() const { return _ptr; }
        Element_ptr& getWritePointer() { return _ptr; }        
        
        /** add an element at the end of the vector (possibly invalidate existing iterators)
		 * @param e element to add
		 */
		void push_back(const Element & e)
        {
            _rep.push_back(e);
            update_storage();            
        }

        // clear memory used by the vector
		void clear(void)
        {
            _rep.clear();
            update_storage();            
        }
        
        // reserve space for the vector
		void reserve(const size_t &m)
        {
            _rep.reserve(m);
            update_storage();            
        }

		void setEntry (size_t i, const Element &a_i){ field().assign(_rep[i],a_i); }
		
		Element &refEntry (size_t i){ return _ptr[i]; }

		const Element &getEntry (size_t i) const { return _ptr[i]; }
		
		Element &getEntry (Element &x, size_t i) const{	return field().assign(x,_ptr[i]); }

		// write
		std::ostream &write ( std::ostream &os, LINBOX_enum(Tag::FileFormat) fmt = Tag::FileFormat::Pretty ) const;
        //read
		std::istream &read ( std::istream &os, LINBOX_enum(Tag::FileFormat) fmt = Tag::FileFormat::Pretty );
        

	};// BlasVector

	template <class Field, class _Rep>
	struct VectorTraits< BlasVector<Field,_Rep> > {
		typedef typename VectorCategories::DenseVectorTag VectorCategory;
		typedef BlasVector<Field,_Rep>                          VectorType;
	};


    template <typename Field, typename Rep>
    bool operator==(const BlasVector<Field,Rep>& v1, const BlasVector<Field,Rep>& v2)
    {
        if (v1.size() != v2.size()) return false;
        auto i1=v1.begin();
        auto i2=v2.begin();
        for( ; i1 != v1.end(); ++i1, ++i2)
            if (*i1 != *i2) return false;
        return true;
    }


	template<>
	Integer BlasVector<Givaro::ZRing<Integer> >::magnitude() const
	{
		Integer max_elt(0);
		for (size_t i = 0 ; i < size() ; ++i)
			if (max_elt < Givaro::abs(_ptr[i]))
				max_elt = Givaro::abs(_ptr[i]) ;
		return max_elt ;
	}


	template<class T>
	std::ostream& operator<< (std::ostream & o, const BlasVector<T> & Mat)
	{
		return Mat.write(o);
	}

} // LinBox

namespace LinBox { /*  BlasSubvector */




	template <class _Vector > // inherit that from owner ?
	class BlasSubvector :
		public Subvector<Subiterator<typename _Vector::Rep::iterator > >
    /* public  BlasVector<typename _Vector::Field, typename _Vector::Rep> */ {
	public :

        typedef _Field                                   Field;
		typedef typename Field::Element                Element;    //!< Element type
        typedef typename Field::Element_ptr        Element_ptr;    //!< Element type
		typedef BlasSubvector<_Vector>                  Self_t;       //!< Self type
		typedef typename _Vector::Storage              Storage;    //!< Actually a <code>std::vector<Element></code> (or alike.)        

        typedef BlasVector<Field,Storage>           vectorType;    //!< vector type
		typedef Self_t                              vectorType;    //!< vector type
		typedef BlasSubvector<Self_t>            subVectorType;    //!< SubVector type
        typedef BlasSubvector<const Self_t> constSubVectorType;    //!< const SubVector type

		typedef typename _Vector::Field                   Field;
		typedef typename Field::Element                 Element;      //!< Element type


		typedef typename Rep::pointer                   pointer;    //!< pointer type to elements
		typedef const pointer                     const_pointer;    //!< const pointer type
		typedef Self_t                            subVectorType;    //!< Subvector type

		typedef BlasVector<Field,Rep>                  blasType;    //!< blas type


	protected:
		Rep &_Vec;                     //!< Parent raw vector
		size_t _size;                   //!< size of Subvector
		size_t _i0;                    //!< beginning of Subvector in \p _Vec
		size_t _inc ;               //!< number of columns in \p _Vec (or stride of \p _Vec)
		const Field & _field;

	public:
		typedef Subvector<Subiterator<typename Rep::iterator > > Father_t;
		typedef typename Father_t::iterator             iterator;
		typedef typename Father_t::const_iterator const_iterator;


		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*  constructors */

		/** NULL constructor.  */
		// BlasSubvector () ;

		BlasSubvector (BlasVector<Field,Rep> &V,
                       size_t ibeg,
                       size_t Stride,
                       size_t Size
                       ) :
			Father_t(),
			_Vec ((V.refRep())),
			_size(Size),_i0 (ibeg),_inc(Stride)
			,_field(V.field())
		{
			setIterators();
		}


		BlasSubvector (const Field & F, std::vector<Element> &V) :
			Father_t(),
			_Vec (V),
			_size(V.size()),_i0 (0),_inc(1)
			,_field(F)
		{
			// not tested
			setIterators();
		}

#if 0 /* impossible */
		BlasSubvector (const Field & F, Subvector<typename Rep::iterator, typename Rep::const_iterator> &V) :
			Father_t(),
			_Vec (V),
			_size(V.size()),_i0 (0),_inc(1)
			,_field(F)
		{
			// not tested
			setIterators();
		}
#endif

		BlasSubvector (const Field & F, const std::vector<Element> &V) :
			Father_t(),
			_Vec (V),
			_size(V.size()),_i0 (0),_inc(1)
			,_field(F)
		{
			std::cout << "oops, copy (?)" << std::endl;
			// not tested
			setIterators();
		}


		BlasSubvector (const Field & F, const std::vector<Element> &V
                       , const size_t ibeg, const size_t Stride, const size_t Size) :
			Father_t(),
			_Vec (V),
			_size(Size),_i0 (ibeg),_inc(Stride)
			// could have _i0 = 0 and start at V+ibeg
			,_field(F)
		{
			// not tested
			setIterators();
		}

		//! @todo subvector of ptr


		BlasSubvector (const BlasSubvector<_Vector> &SV,
                       size_t ibeg,
                       size_t Stride,
                       size_t Size
                       ) :
			Father_t(),
			_Vec (SV._Vec),
			_size(Size),_i0 (SV._i0+SV._inc*ibeg)
			,_inc(SV._inc*Stride)
			,_field(SV.field())
		{
			// not tested
			setIterators();
		}

		/** Copy constructor.
		 * @param SM Subvector to copy
		 */
		BlasSubvector (const BlasSubvector<_Vector> &SV) :
			Father_t(),
			_Vec (SV._Vec),
			_size(SV._size),_i0 (SV._i0)
			,_inc(SV._inc)
			,_field(SV.field())
		{
			// not tested
			setIterators();
		}

#if 0 /*  from BlasMatrix (should be a Row/Col in BlasMatrix, not here... */
		BlasSubvector (const BlasMatrix<Field,Rep> &M
                       , size_t ibeg
                       , LINBOX_enum (Tag::Direction) f ) :
			Father_t(),
			_Vec (const_cast<Rep&>(M.refRep()))
			,_size((f==Tag::Direction::Row)?(M.coldim()):(M.rowdim()))
			,_i0 ((f==Tag::Direction::Row)?(ibeg*M.coldim()):(ibeg))
			,_inc((f==Tag::Direction::Row)?(1):(M.coldim()))
			,_field(M.field())
		{
			setIterators();
		}

		template<class _Matrix>
		BlasSubvector (const BlasSubmatrix<_Matrix> &M
                       , size_t ibeg
                       , LINBOX_enum (Tag::Direction) f ) :
			Father_t(),
			_Vec (const_cast<Rep&>(M.refRep()))
			,_size((f==Tag::Direction::Row)?(M.coldim()):(M.rowdim()))
			,_i0 ((f==Tag::Direction::Row)?(M.offset()+ibeg*M.stride()):(M.offset()+ibeg))
			,_inc((f==Tag::Direction::Row)?(1):(M.stride()))
			,_field(M.field())
		{
			setIterators();
		}
#endif


		//! @todo more general subvectors

		/*  Members  */

		BlasSubvector &operator = (const BlasSubvector<_Vector> &SV)
		{
			if ( &SV == this)
				return *this ;
			_Vec=SV._Vec ; //!@todo use functions, not =
			_size = SV.size();
			_i0=SV._i0;
			_inc = SV.stride();
			// _field = SV.field();
			linbox_check(field().characteristic() == SV.field().characteristic());
			setIterators();
			return *this ;

		}

		// template<typename _Tp1>
		// struct rebind ;

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		size_t size() const {return _size;}

		size_t getStride() const {return _inc ; }
		size_t stride() const { return getStride() ;}


		///////////////////
		//      I/O      //
		///////////////////

		// std::ostream &write (std::ostream &os) const;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		pointer getPointer() const { return &(_Vec[_i0]); }

		const_pointer &getConstPointer() const { return &(_Vec[_i0]); }


		pointer& getWritePointer() { return &(_Vec[_i0]); }


		const Element& setEntry (size_t i, const Element &a_i)
		{
			return _Vec[_i0+i*_inc] = a_i;
		}

		Element &refEntry (size_t i)
		{
			return _Vec[_i0+i*_inc] ;
		}

		const Element &getEntry (size_t i) const
		{
			return _Vec[_i0+i*_inc] ;
		}

		Element &getEntry (Element &x, size_t i)
		{
			x = _Vec[_i0+i*_inc] ;
			return x;
		}


#if 0 /* using Father_t */
		///////////////////
		//   ITERATORS   //
		///////////////////

		class Iterator  ;
		class ConstIterator ;

		class IndexedIterator ;
		class ConstIndexedIterator ;
#endif


		const Field& field() const { return _field ;}
		// Field & field() { return _field; }

		// write (same as BlasVector)
		std::ostream &write ( std::ostream &os, LINBOX_enum(Tag::FileFormat) fmt = Tag::FileFormat::Pretty ) const
		{
			switch(fmt) {
			case (Tag::FileFormat::Pretty) :
				{
					os << '[' ;
					for(const_iterator it= this->begin();it != this->end(); ++it)
						field().write(os, *it) << " ";
					return	os << ']' ;
				}
			case (Tag::FileFormat::Maple) :
				{
					os << '<' ;
					for(const_iterator it=this->begin();it != this->end(); ) {
						field().write(os, *it);
						++it ;
						if (it != this->end())
							os << ',' ;
					}
					return	os << '>' ;
				}

			default :
				return os << "not implemented" ;
			}
		}

	private:
		void setIterators()
		{
			Father_t::_begin = iterator (_Vec.begin()+(ptrdiff_t)_i0 , (ptrdiff_t)_inc);
			Father_t::_end   = iterator (_Vec.begin()+(ptrdiff_t)(_i0 + _size*_inc) , (ptrdiff_t)_inc);
		}

	};

	template<class T>
	std::ostream& operator<< (std::ostream & o, const BlasSubvector<T> & Mat)
	{
		return Mat.write(o);
	}

} // LinBox

namespace LinBox { /*  traits */

	// this could also be a member of BlasVector
	template<class _Field, class _Rep>
	struct ContainerTraits<BlasVector<_Field,_Rep> > {
		typedef ContainerCategories::Vector ContainerCategory ;
	};

	//! @todo remove vectors
	template<class _Rep>
	struct ContainerTraits<std::vector<_Rep> > {
		typedef ContainerCategories::Vector ContainerCategory ;
	};

	template<class _Vector>
	struct ContainerTraits<BlasSubvector<_Vector> > {
		typedef ContainerCategories::Vector ContainerCategory ;
	};

}

#include "blas-vector.inl"

#endif // __LINBOX_vector_blas_vector_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
