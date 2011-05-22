/* Copyright (C) 2007 LinBox
 * Written by JG Dumas
 *
 *
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_lazy_product_H
#define __LINBOX_lazy_product_H


namespace LinBox 
{

        // Lazy computation of the product of the moduli
    struct LazyProduct : public std::vector< Integer > {
        typedef std::vector< Integer > Father_t;
    protected:
        bool                _tobecomputed;
    public:

        LazyProduct() : Father_t(), _tobecomputed(true) {}

        void initialize(const Integer& i) {
            _tobecomputed = false;
            this->resize(0);
            this->push_back(i);
        }
            
        bool mulin(const Integer& i) {
            if (this->size()) {
                if (i != this->back()) {
                    this->push_back( i );
                    return _tobecomputed = true;
                } else {
                    return _tobecomputed;
                }
            
            } else {
                this->push_back( i );
                return _tobecomputed = false;
            }
        }
      
        bool mulin(const LazyProduct& i) {
            this->insert(this->end(), i.begin(), i.end());
            return _tobecomputed = (this->size()>1);
        }
      
        Integer & operator() () {
            if (_tobecomputed) {
                Father_t::const_iterator iter = this->begin();
                Father_t::iterator       prod = this->begin();
                for(++iter; iter != this->end(); ++iter)
                    *prod *= *iter;
                this->resize(1);
                _tobecomputed = false;
            }
            return this->back();
        }

        bool noncoprime(const Integer& i) const {
            Integer g;
            for(Father_t::const_iterator iter = this->begin(); iter != this->end(); ++iter)
                if ( gcd(g,i,*iter) > 1) return true;
            return false;
        }   
       
        friend std::ostream& operator<< (std::ostream& o, const LazyProduct& C) {
            o << "{";
            for(Father_t::const_iterator refs = C.begin();
                refs != C.end() ;
                ++refs )
                o << (*refs) << " " ;
            return o << "}";
        }
        
    };
    
}


#endif //__LINBOX_lazy_product_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
