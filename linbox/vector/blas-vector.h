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


/*! @file vector/blas-vector.h
 * @ingroup vector
 * A \c BlasVector<\c _Field > represents a vector as an array of
 * <code>_Field::Element</code>s.
 *
 */

#ifndef __LINBOX_vector_blas_vector_H
#define __LINBOX_vector_blas_vector_H

#include <vector>
#include <iterator>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-tags.h"
#include "linbox/field/hom.h"
#include "linbox/field/field-traits.h"
#include "linbox/field/rebind.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/blas-subvector.h"

#include "fflas-ffpack/fflas/fflas.h"

namespace LinBox {

    template<class _Field, class _Storage>
    class BlasVector {

    public:
        typedef _Field                                   Field;
        typedef typename Field::Element                Element;    //!< Element type
        typedef Element                             value_type;
        typedef typename Field::Element_ptr        Element_ptr;    //!< Element type
        typedef Element_ptr                            pointer;
        typedef typename Field::ConstElement_ptr   ConstElement_ptr;    //!< Element type
        typedef _Storage                               Storage;    //!< Actually a <code>std::vector<Element></code> (or alike.)
        typedef typename Storage::reference          reference;
        typedef typename Storage::const_reference  const_reference;
        typedef BlasVector<_Field,_Storage>             Self_t;    //!< Self type
        typedef Self_t                              vectorType;    //!< vector type
        typedef BlasSubvector<Self_t>            subVectorType;    //!< SubVector type
        typedef BlasSubvector<const Self_t> constSubVectorType;    //!< const SubVector type


        /* iterators */
        typedef typename Storage::iterator                           iterator;
        typedef typename Storage::const_iterator               const_iterator;
        typedef std::reverse_iterator<iterator>	             reverse_iterator;
		typedef std::reverse_iterator<const_iterator>  const_reverse_iterator;

    protected:
        size_t                 _size=0;
        Storage                 _rep;
        Element_ptr             _ptr;
        const Field          &_field;

    public:

        void resize (size_t n){
// std::cout<<"BlasVector resize: "<<_ptr<<" ("<<_size<<") to ";
            _rep.resize(n);
            _ptr = _rep.data();
// std::cout<<_ptr<<" ("<<n<<")"<<std::endl;
            for (size_t i=_size;i<n;i++)
                field().init(_rep[i]);
            _size = n;
// std::cout<<"BlasVector resize end."<<std::endl;
        }

        void resize (size_t n, const Element& val){
            _rep.resize(n);
            _ptr = _rep.data();
            for (size_t i=_size;i<n;i++){
                field().init(_rep[i]);
                field().assign(_rep[i],val);
            }
            _size = n;
        }

         //////////////////
        // CONSTRUCTORS //
        //////////////////


        /*! Allocates a new \f$ 0 \f$ vector (shaped and ready).*/
        BlasVector (const _Field &F) : _field(F) {
            resize(0);
        }

        /*! Copy Constructor of a vector (copying data).
         * @param V vector to be copied.
         */
        BlasVector (const Self_t &V) : _field(V.field()) {
            resize(V.size());
            auto jt(V.begin());
            for (auto it=_rep.begin();it!=_rep.end();++it,++jt)
                field().assign(*it,*jt);
        }


        // This constructor is templated because if SizeType was a size_t, then calling it with a signed const litteral
        // (e.g. v(F,3)) would fail to launch this constructor, but rather go in the templated one (_Field, VectorBase)
        // on OSX where long can not be cast to size_t).
        /*! Allocates a new \f$ m \f$ zero vector (shaped and ready).
         * @param F
         * @param m vector dimension
         */
        template<class SizeType, typename std::enable_if<std::is_arithmetic<SizeType>::value, int>::type=0>
        BlasVector (const _Field &F, const SizeType &m) : _field(F) {
            resize(m);
        }

        template<class SizeType, typename std::enable_if<std::is_arithmetic<SizeType>::value, int>::type=0>
        BlasVector (const _Field &F, const SizeType &m, const Element& e) : _field(F) {
            resize(m,e);
        }

        /*! Create a BlasVector from an iterator of elements
         * @param F Field of the created vector
         * @param it iterator to be copied from
         * @param m vector dimension
         * @param inc increment value for iterating
         */
        template<typename ConstIterator, typename std::enable_if<!std::is_arithmetic<ConstIterator>::value, int>::type=0>
        BlasVector(const _Field & F, const ConstIterator& jt, const size_t m) :  _field(F) {
            resize(m);
            ConstIterator jtt(jt);
            for (auto it=_rep.begin();it!=_rep.end();++it,++jtt)
                field().assign(*it,*jtt);
        }

         /*! Create a BlasVector from a container of elements
         * @param F Field of the created vector
         * @param V iterator to be copied from
         */
        template<class VectorBase, typename std::enable_if<!std::is_arithmetic<VectorBase>::value, int>::type=0>
        BlasVector (const _Field & F, const VectorBase & V)  :  _field(F) {
            resize(V.size());
            typename VectorBase::const_iterator jt=V.begin();
            for (auto it=_rep.begin();it!=_rep.end();++it,++jt)
                field().assign(*it,*jt);
        }


        /*! Create a BlasVector from another vector defined over a different field (use homomorphism if it exists)
         * @param F Field of the created vector
         * @param V Vector to be copied
         */
        template<class OtherVector>//, typename std::enable_if<!std::is_arithmetic<VectorBase>::value, int>::type=0>
        BlasVector (const OtherVector & V, const _Field & F) :  _field(F) {
            resize(V.size());
            typename OtherVector::template rebind<_Field>()(*this,V);
        }

        //! operator = (copying data)
        Self_t& operator= (const Self_t& A) {
            if (this!=&A){
                linbox_check (areFieldEqual(_field, A.field()));
                resize(A.size());
                FFLAS::fassign(field(),_size,A._ptr,1,_ptr,1);
            }
            return *this;
        }

        //! operator = (copying data from different vector type)
        template<class _Vector>
        Self_t& operator= (const _Vector& A){
            if (this!=&A){
                linbox_check (_field !=A.field());
                resize(A.size());
                for (auto it=_rep.begin(), jt=A.begin();it!=_rep.end();it++,jt++)
                    field().assign(*it,*jt);
            }
            return *this;
        }


        template<class _Vector>
        Self_t& copy(const _Vector& A){
            return *this=A;
        }

        //! Rebind operator
        template<typename _Tp1, typename _Rep2 = typename Rebind<Storage, _Tp1>::other>
        struct rebind {
            typedef BlasVector<_Tp1, _Rep2> other;

            void operator() (other & Ap, const Self_t& A) {
                typedef typename Self_t::const_iterator ConstSelfIterator ;
                typedef typename other::iterator OtherIterator ;
                OtherIterator    Ap_i = Ap.begin();
                ConstSelfIterator A_i = A.begin();
                Hom<Field, _Tp1> hom(A. field(), Ap. field()) ;
                for ( ; A_i != A. end(); ++ A_i, ++ Ap_i)
                    hom.image (*Ap_i, *A_i);
            }
        };


        /////////////////
        //  ACCESSORS  //
        /////////////////

        const _Field& field() const { return _field;}

        // dimension of the vector
        size_t size() const{ return _size; }
        size_t max_size() const{ return _size; }

        /*!_@internal
         * Get access to the vector  data.
         */
        ConstElement_ptr getPointer() const { return _ptr; }
        Element_ptr      getPointer()       { return _ptr; }
        ConstElement_ptr getConstPointer() const { return _ptr;}

        const Storage &getRep() const { return _rep ; }

        /** Get the increment in the vector
         * @return the inc value of the subvector
         */
        size_t getInc() const {return 1;}


        /** add an element at the end of the vector (possibly invalidate existing iterators)
         * @param e element to add
         */
        void push_back(const Element & e) {
            _rep.push_back(e);
            _ptr= _rep.data();
            _size++;
        }

        // clear memory used by the vector
        void clear(void) {
            _rep.clear();
            _ptr=_rep.data();
            _size=0;
        }

        // reserve space for the vector
        void reserve(const size_t &m) {
            _rep.reserve(m);
            _ptr=_rep.data();
        }

        void setEntry (size_t i, const Element &a_i){
// std::cout<<"BV: "<<" "<<&(*_ptr)<<" "<<i<<" "<<a_i<<std::endl;
//             field().assign(_rep.at(i),a_i);
            field().assign(_rep[i],a_i);
        }

        reference refEntry (size_t i){ return _ptr[i]; }

        const_reference getEntry (size_t i) const { return _ptr[i]; }

        Element& getEntry (Element &x, size_t i) const{ return field().assign(x,_ptr[i]); }

        // write
        std::ostream &write ( std::ostream &os, Tag::FileFormat fmt = Tag::FileFormat::Pretty ) const {
            switch(fmt) {
            case (Tag::FileFormat::Pretty) :
                {
                    os << '[' ;
                    for(auto it= _rep.begin();it != _rep.end(); ++it)
                        field().write(os, *it) << " ";
                    return	os << ']' ;
                }
            case (Tag::FileFormat::Maple) :
                {
                    os << '<' ;
                    for(auto it= _rep.begin();it != _rep.end(); ++it) {
                        field().write(os, *it);
                        if (it != this->end()-1)
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

        // Iterators
        iterator               begin  (void)       { return _rep.begin(); }
        const_iterator         begin  (void) const { return _rep.begin(); }
        iterator               end    (void)       { return _rep.end(); }
        const_iterator         end    (void) const { return _rep.end(); }
        reverse_iterator       rbegin (void)       { return reverse_iterator (end()); }
		const_reverse_iterator rbegin (void) const { return reverse_iterator (end()); }
		reverse_iterator       rend   (void)       { return reverse_iterator (begin()); }
		const_reverse_iterator rend   (void) const { return reverse_iterator (begin()); }
        iterator               Begin  (void)       { return this->begin(); }
        const_iterator         Begin  (void) const { return this->begin(); }
        iterator               End    (void)       { return this->end(); }
        const_iterator         End    (void) const { return this->end(); }

        // Element access
        reference       operator[] (size_t n)       { return _rep[n]; }
        const_reference operator[] (size_t n) const { return _rep[n]; }

        reference       at(size_t n) { return _rep.at(n);}
        const_reference at(size_t n) const {return _rep.at(n);}

        reference       front (void)       { return _rep.front(); }
        const_reference front (void) const { return _rep.front(); }
        reference       back  (void)       { return _rep.back(); }
        const_reference back  (void) const { return _rep.back; }

        bool empty() const {return _rep.empty();}


        // Miscelleanous
        void random(){
            typename _Field::Element x; field().init(x);
            typename _Field::RandIter r(field());
            for (size_t i = 0; i < size(); ++i)
                setEntry(i, r.random(x));
        }

        template<class RandIter>
        void random( RandIter r){
            typename _Field::Element x; field().init(x);
            for (size_t i = 0; i < size(); ++i){
                setEntry(i, r.random(x));
            }
        }

    };


    template <class Field, class Storage>
    bool operator==(const BlasVector<Field,Storage>& v1, const BlasVector<Field,Storage>& v2) {
                if (v1.size() != v2.size()) return false;
                auto i1=v1.begin();
                auto i2=v2.begin();
                for( ; i1 != v1.end(); ++i1, ++i2)
                        if (*i1 != *i2) return false;
                return true;
    }

    template <class Field, class Storage>
    std::ostream& operator<< (std::ostream & os, const BlasVector<Field,Storage> & V) {
        return V.write(os);
    }

    template <class Field, class _Rep>
    struct VectorTraits< BlasVector<Field,_Rep> > {
        typedef typename VectorCategories::DenseVectorTag VectorCategory;
        typedef BlasVector<Field,_Rep>                        VectorType;
    };


    template<class _Vector>
    struct VectorTraits< BlasSubvector<_Vector> > {
        typedef typename VectorCategories::DenseVectorTag VectorCategory;
        typedef BlasSubvector<_Vector>                        VectorType;
    };


    // this could also be a member of BlasVector
    template<class _Field, class _Storage>
    struct ContainerTraits<BlasVector<_Field,_Storage> > {
        typedef ContainerCategories::Vector ContainerCategory ;
    };

    template<class _Vector>
    struct ContainerTraits<BlasSubvector<_Vector> > {
        typedef ContainerCategories::Vector ContainerCategory ;
    };




} // LinBox
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
