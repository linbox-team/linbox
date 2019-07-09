/* linbox/matrix/blas-vector.h
 * Copyright (C) 2013 the LinBox group
 *               2019 Pascal Giorgi 
 *
 * Written by :
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


/*! @file vector/dense-vector.h
 * @ingroup vector
 * A \c DenseVector<\c _Field > represents a vector as an array of
 * <code>_Field::Element</code>s.
 *
 */

#ifndef __LINBOX_vector_dense_subvector_H
#define __LINBOX_vector_dense_subvector_H

#include "linbox/vector/subiterator.h"

namespace LinBox {



    template <typename _Vector>
    class VectorEltPointer {
    public:
        typedef typename _Vector::Field::Element_ptr     pointer;
    };
    template <typename _Vector>
    class VectorEltPointer <const _Vector> {
    public:
        typedef typename _Vector::Field::ConstElement_ptr pointer;
    };
    
    template<class _Vector>
    class BlasSubvector {

    public:
        typedef typename _Vector::_Field                   Field;
        typedef typename Field::Element                  Element;    //!< Element type
        typedef typename _Vector::Storage                Storage;    //!< Actually a <code>std::vector<Element></code> (or alike.)
        typedef BlasSubvector<_Vector>                    Self_t;    //!< Self type
        typedef typename VectorEltPointer<_Vector>::pointer                    pointer;    //!< pointer type to elements
        typedef typename VectorEltPointer<const _Vector>::pointer        const_pointer;    //!< const pointer type to elements
        typedef _Vector                               vectorType;    //!< vector type
        typedef Self_t                             subVectorType;    //!< SubVector type
        typedef BlasSubvector<const _Vector> constSubVectorType;    //!< const SubVector type


        /* Forward declaration of iterators */
        typedef Subiterator<pointer>             iterator;
        typedef Subiterator<const_pointer> const_iterator;


    protected:
		pointer	     		    _ptr;
        size_t			       _size;
        size_t                  _inc;
		const Field		     &_field;

    public:


        //////////////////
		// CONSTRUCTORS //
		//////////////////

        /** Constructor from an existing @ref BlasVector and dimensions.
         * \param V Pointer to @ref BlasVector of which to construct subvector
         * \param beg Starting idx
         * \param dim dimension
         * \param inc distance between two element
         */
        BlasSubvector (vectorType &V,
                        size_t beg,
                        size_t dim,
                        size_t inc);

        /** Constructor from an existing @ref BlasSubvector and dimensions.
         * \param V Pointer to @ref DenseSubector of which to construct subvector
         * \param beg Starting idx
         * \param dim dimension
         * \param inc distance between two element
         */
        BlasSubvector (Self_t &V,
                        size_t beg,
                        size_t dim,
                        size_t inc);

        
        /** Constructor from an existing @ref BlasMatrix
         * \param V Pointer to @ref BlasVector of which to construct submatrix
         */
        BlasSubvector (vectorType &V);

        /** Constructor from a raw pointer
         * \param ptr Pointer to @ref first element of the vector of which to construct submatrix
         */
        BlasSubvector (const Field& F,
                       pointer ptr,
                       size_t dim,
                       size_t inc);
        

        
		//! Rebind operator
        template <typename _TP1, typename _Rep2 = typename Rebind<Storage, _TP1>::other >
        struct rebind;

        

		/////////////////
		//  ACCESSORS  //
		/////////////////

        const Field& field() const { return _field;}
        
        // dimension of the vector
        size_t size() const{ return _size; }

        /** Get access to the vector  data.
         * @return a pointer to the first element of the vector : constness is according to _Vector constness
		 */
		const_pointer getPointer() const { return _ptr; }
        pointer       getPointer()       { return _ptr; }
        const_pointer getConstPointer() const { return _ptr;}

        /** Get the increment in the vector
         * @return the inc value of the subvector
         */
        size_t getInc() const {return _inc;}
        

		void setEntry (size_t i, const Element &a_i){ field().assign(_ptr[i],a_i); }
		
		Element &refEntry (size_t i){ return _ptr[i]; }

		const Element &getEntry (size_t i) const { return _ptr[i]; }
		
		Element &getEntry (Element &x, size_t i) const{	return field().assign(x,_ptr[i]); }

		// write
		std::ostream &write ( std::ostream &os, Tag::FileFormat fmt = Tag::FileFormat::Pretty ) const {
			switch(fmt) {
			case (Tag::FileFormat::Pretty) :
				{
					os << '[' ;
					for(size_t i=0; i<_size; i++)
						field().write(os, *(_ptr+_inc*i)) << " ";
                    return	os << ']' ;
				}
			case (Tag::FileFormat::Maple) :
				{
					os << '<' ;
                    for(size_t i=0; i<_size; i++){ 
						field().write(os, *(_ptr+_inc*i));
						if (i != _size-1)
							os << ',' ;
					}
					return	os << '>' ;
				}
			default :
				return os << "not implemented" ;
			}
        }
        
        //read
		std::istream &read ( std::istream &is, Tag::FileFormat fmt = Tag::FileFormat::Pretty ) {
            return is;
        }

        ///////////////////
		//   ITERATORS   //
		///////////////////
		iterator               begin  (void)       { return iterator(_ptr,_inc); }
		const_iterator         begin  (void) const { return iterator(_ptr,_inc); }
		iterator               end    (void)       { return iterator(_ptr+_inc*_size,_inc); }
		const_iterator         end    (void) const { return iterator(_ptr+_inc*size,_inc); }

		// Element access
		Element&       operator[] (size_t n)       { return _ptr[n*_inc]; }
		const Element& operator[] (size_t n) const { return _ptr[n*_inc]; }

		Element& at(size_t n) {
            if (n<_size)
                return _ptr[n*_inc];
            else throw std::out_of_range("out of range"); //out of range error message.
        }

		const Element&  at(size_t n) const {
            if (n<_size)
                return _ptr[n*_inc];
            else throw std::out_of_range("out of range"); //out of range error message.
        }

		Element&       front (void)       { return *_ptr;}
		const Element& front (void) const { return *_ptr;}
		Element&       back  (void)       { return *(_ptr+(_size-1)*_inc);}
		const Element& back  (void) const { return *(_ptr+(_size-1)*_inc);}


    };
    
    template <class Vector>
    std::ostream& operator<< (std::ostream & os, const BlasSubvector<Vector> & V) {
		return V.write(os);
	}



    
} // LinBox
#endif
// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
