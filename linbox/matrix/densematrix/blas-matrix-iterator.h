/*
 * Copyright (C) 2018 Pascal Giorgi 
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file matrix/densematrix/blas-matrix-iterator.h
 * @ingroup densematrix
 *
 */

#ifndef __LINBOX_matrix_densematrix_blas_matrix_iterator_H
#define __LINBOX_matrix_densematrix_blas_matrix_iterator_H

#include "linbox/vector/blas-vector.h"

namespace LinBox {
    
    template < class _Field, class _Storage, class _subVectorType>
    class BlasMatrixIterator {
    public:

         BlasMatrixIterator(){}
  
        BlasMatrixIterator (const _Field & F, typename _subVectorType::pointer ptr,  size_t size, size_t inc, size_t dist) :
            _data(F, ptr, inc, size), _dist(dist) {
            // std::cout<<"BLASMATRIX ITERATOR CST: "<<_data<<std::endl;
            //std::cout<<"BLASMATRIX ITERATOR CST: "<<ptr<<" "<<size<<" "<<inc<<" "<<dist<<" "<<&F<<std::endl;
        }
        
        BlasMatrixIterator(const BlasMatrixIterator &it) : _data(it._data), _dist(it._dist) {}

        BlasMatrixIterator& operator=(const BlasMatrixIterator &it) {
            if (&it != this){
                _data = it._data;
                _dist = it._dist;
            }
            return *this;
        }

        BlasMatrixIterator& operator --() {
            _data= _subVectorType (_data.field(), _data.getPointer() - _dist, _data.getInc(), _data.size());
            return *this;
        }

        
        BlasMatrixIterator  operator-- (int) {
            BlasMatrixIterator tmp (*this);
            --*this;
            return tmp;
        }

        BlasMatrixIterator& operator++ (){
            // std::cout<<"BLASMATRIX ITERATOR ++: "<<_data.getPointer() + _dist<<"->";
            // std::cout<<_data<<std::endl;
            _data= _subVectorType (_data.field(), _data.getPointer() + _dist, _data.getInc(), _data.size());

            return *this;
        }

        BlasMatrixIterator  operator++ (int){
            BlasMatrixIterator tmp (*this);
            ++*this;
            return tmp;
        }

        BlasMatrixIterator operator+ (int i) {
            return BlasMatrixIterator (_data.field(), _data.getPointer() + i*_dist, _data.size (), _data.getInc(), _dist);
        }

        BlasMatrixIterator& operator += (int i) {
            _data= _subVectorType(_data.field(), _data.getPointer() + i*_dist, _data.getInc(), _data.size ());
            return *this;
        }

        _subVectorType operator[] (int i) const {
            return _subVectorType (_data.field(), _data.getPointer() + i*_dist, _data.getInc(), _data.size ());
        }

        _subVectorType* operator-> () { return &_data;}


        _subVectorType& operator* () {
            //std::cout<<"BLASMATRIX ITERATOR *: "<<_data<<std::endl;
            return _data;
        }


        bool operator!= (const BlasMatrixIterator& c) const
        {

            //std::cout<<"Checking BMIter :"<<_data.getPointer()<<" <> "<< c._data.getPointer()<<std::endl;
            return  _data.getPointer()!= c._data.getPointer();
        }
        
        operator BlasMatrixIterator<_Field,_Storage,typename _subVectorType::constSubVectorType> ()
        {
            return BlasMatrixIterator<_Field,_Storage, typename _subVectorType::constSubVectorType>(_data.field(), _data.getPointer(), _data.size(), _data.getInc(),_dist);
        }
        
    private:
        _subVectorType  _data;
        size_t          _dist;
    };




    template < class _Field, class _Pointer, class _Element>
    class BlasMatrixIndexedIterator {
        size_t _r_index;
        size_t _c_index;
        size_t _c_dim;
        size_t _stride;
        _Pointer _begin;
        typedef _Element value_type;

    public:
        BlasMatrixIndexedIterator (const size_t  &c_dim,
                                   const size_t  &stride, 
                                   const size_t  &r_index,
                                   const size_t  &c_index,
                                   const _Pointer &begin) :
            _r_index (r_index), _c_index (c_index), _c_dim (c_dim), _stride(stride), _begin (begin) {}

        BlasMatrixIndexedIterator () :
            _r_index (0), _c_index (0), _c_dim (1), _begin (0) {}

       
		BlasMatrixIndexedIterator (const BlasMatrixIndexedIterator& r) :
			_r_index (r._r_index), _c_index (r._c_index), _c_dim (r._c_dim), _stride(r._stride), _begin (r._begin) {}
        
		BlasMatrixIndexedIterator& operator = (const BlasMatrixIndexedIterator &iter) {
			_r_index = iter._r_index;
			_c_index = iter._c_index;
			_c_dim   = iter._c_dim;
            _stride  = iter._stride,
			_begin   = iter._begin;
			return *this;
		}

		bool operator == (const BlasMatrixIndexedIterator &iter) const {
			return (_r_index == iter._r_index) &&
                (_c_index == iter._c_index) &&
                (_c_dim == iter._c_dim) &&
                (_stride == iter._stride) &&
                (_begin==iter._begin);
		}

		bool operator != (const BlasMatrixIndexedIterator& iter) const {
			return (_r_index != iter._r_index) ||
                (_c_index != iter._c_index) ||
                (_c_dim != iter._c_dim) ||
                (_stride != iter._stride) ||
                (_begin!=iter._begin);
		}

		BlasMatrixIndexedIterator &operator ++ () {
			++_c_index;

			if (_c_index == _c_dim) {
				_c_index = 0;
				++_r_index;
			}

			return *this;
		}


		BlasMatrixIndexedIterator operator ++ (int) {
			BlasMatrixIndexedIterator tmp = *this;
			++(*this);
			return tmp;
		}

		BlasMatrixIndexedIterator &operator -- () {
			if (_c_index)
				--_c_index;
			else {
				--_r_index;
				_c_index = _c_dim - 1;
			}

			return *this;
		}

		BlasMatrixIndexedIterator operator -- (int) {
			BlasMatrixIndexedIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * () const {
			return *(_begin + _r_index * _stride + _c_index);
		}

		value_type * operator -> () const {
			return _begin + _r_index * _stride + _c_index;
		}

		size_t rowIndex () const {
			return _r_index;
		}

		size_t colIndex () const {
			return _c_index;
		}

		const value_type &value () const {
			return *(_begin + _r_index * _stride + _c_index);
		}


	};


    
} // end namespace LinBox

#endif // __LINBOX_matrix_densematrix_blas_matrix_iterator_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


