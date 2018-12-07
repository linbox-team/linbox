/* linbox/matrix/blas-vector.inl
 * Copyright (C) 2013 the LinBox group
 *               2018 revamped by Pascal Giorgi 
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

/*! @internal
 * @file vector/blas-vector.inl
 * @ingroup vector
 *
 */


#ifndef __LINBOX_vector_blas_vector_INL
#define __LINBOX_vector_blas_vector_INL

namespace LinBox {
   
    template<class _Field, class _Storage>
    template<typename Vector>
    void  BlasVector<_Field,_Storage>::init(){
        for (auto &it: _rep)
            _field.init(it) ;
    }
    
    template<class _Field, class _Storage>
    template<typename Vector>
    void  BlasVector<_Field,_Storage>::initBlasVector(const Vector & V){
        iterator it = _rep.begin();
        typename Vector::const_iterator jt = V.begin();
        for ( ; it != _rep.end(); ++it, jt++)
            _field.init(*it, *jt) ;
    }
    
    template<class _Field, class _Storage>
    template<typename Vector>
    void 	BlasVector<_Field,_Storage>::assignBlasVector(const Vector & V){
        iterator it = _rep.begin();
        typename Vector::const_iterator jt = V.begin();
        for ( ; it != _rep.end(); ++it, jt++){
            _field.init(*it);
            _field.assign(*it,*jt);
        }
    }

    template<class _Field, class _Storage>
    BlasVector<_Field,_Storage>::BlasVector (const Self_t &V) : Father_t(),
                                   _size(V.size()),_rep(V.size()),_ptr(_rep.data()), _field(V.Field())
    {
        updateIterators();
        assignBlasVector(V);
    }


            
    template<class _Field, class _Storage>
    BlasVector<_Field,_Storage>::BlasVector (const _Field &F) : Father_t(),
                                   _size(0),_rep(0),_ptr(_rep.data()), _field(F)
    {
        updateIterators();
    }

    template<class _Field, class _Storage>
    template<class SizeType, typename std::enable_if<std::is_arithmetic<SizeType>::value, int>::type=0>
    BlasVector<_Field,_Storage>::BlasVector (const _Field &F, const SizeType &m,const Element& val)
        : Father_t(),_size(m),_rep(m),_ptr(_rep.data()), _field(F)
    {
        updateIterators();        
        for (auto& it:_rep){
            field().init(it);
            field().assign(it,val);
        }
    }

    template<class _Field, class _Storage>
    template<class SizeType, typename std::enable_if<std::is_arithmetic<SizeType>::value, int>::type=0>
    BlasVector<_Field,_Storage>::BlasVector (const _Field &F, const SizeType &m)
        : Father_t(),_size(m),_rep(m),_ptr(_rep.data()), _field(F)
    {
        updateIterators();
        init();
    }

    template<class _Field, class _Storage>
    template<typename ConstIterator>
    BlasVector<_Field,_Storage>::BlasVector(const _Field & F, const ConstIterator& it, const size_t m)
        : Father_t(),_size(m),_rep(m),_ptr(_rep.data()), _field(F)
    {
        updateIterators();
        init();
    }


    template<class _Field, class _Storage>
    template<class VectorBase>//, typename std::enable_if<!std::is_arithmetic<VectorBase>::value, int>::type=0>
    BlasVector<_Field,_Storage>::BlasVector (const VectorBase & V, const _Field & F);
        

    template<class _Field, class _Storage>
    Self_t& BlasVector<_Field,_Storage>::operator= (const Self_t& A) ;
    
    template<class _Field, class _Storage>
    template<class _Matrix>
    Self_t& BlasVector<_Field,_Storage>::operator= (const _Matrix& A);
    
    


    


    
} // namespace LinBox 

#endif // __LINBOX_vector_blas_vector_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
