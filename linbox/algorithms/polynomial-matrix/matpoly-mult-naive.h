/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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


#ifndef __LINBOX_MATPOLY_MUL_NAIVE_H
#define __LINBOX_MATPOLY_MUL_NAIVE_H

#include <typeinfo>

#include "linbox/util/error.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/polynomial-matrix.h"
#include <algorithm>
namespace LinBox
{
        
	template <class _Field>
	class PolynomialMatrixNaiveMulDomain {
	private:
		const _Field            *_field;
		BlasMatrixDomain<_Field>   _BMD;

	public:
        typedef _Field Field; 

        inline const Field& field() const { return *_field; }

        PolynomialMatrixNaiveMulDomain(const Field &F) :
			_field(&F), _BMD(F){}
                

        // c must be allocated with the right size
        template<typename Matrix1,typename Matrix2,typename Matrix3>
        void mul(Matrix1& c, const Matrix2&  a, const Matrix3&  b) const
        {
            //std::clog<<"Mul matpoly: "<<a.rowdim()<<"x"<<a.coldim()<<" by "<<b.rowdim()<<"x"<<b.coldim()<<std::endl;
            // std::clog<<"A="<<a<<std::endl;
            // std::clog<<"B="<<b<<std::endl;
            // std::clog<<"C="<<c<<std::endl;
            
            for (size_t k=0;k<a.size()+b.size()-1;k++){
                auto c_tmp=c[k];
                size_t idx_min= (k+1<b.size()?0:k+1-b.size());
                size_t idx_max=std::min(k,a.size()-1);
                // a[idx_min].write(std::clog, Tag::FileFormat::Plain)<<std::endl;
                // b[k-idx_min].write(std::clog, Tag::FileFormat::Plain)<<std::endl;;

                _BMD.mul(c_tmp,a[idx_min],b[k-idx_min]);               
                for (size_t i=idx_min+1;i<=idx_max;i++){
                    // std::clog<<"MUL["<<k<<"]=";
                    // a[i].write(std::clog, Tag::FileFormat::Plain)<<std::endl;
                    // b[k-i].write(std::clog, Tag::FileFormat::Plain)<<std::endl;;
                    _BMD.axpyin(c_tmp,a[i],b[k-i]);
                }
                // std::clog<<"c_tmp["<<k<<"]=";c_tmp.write(std::clog, Tag::FileFormat::Plain)<<std::endl;
                c.setMatrix(c_tmp,k); 
            }
            // std::clog<<"C="<<c<<std::endl;
            // std::clog<<"-----------"<<std::endl;
        }                          
        
        template<typename Matrix1,typename Matrix2,typename Matrix3>
        void midproduct(Matrix1& c, const Matrix2&  a, const Matrix3&  b,
                        bool smallLeft=true, size_t n0=0,size_t n1=0) const
        {
            //cout<<"naive midprod "<<a.size()<<"x"<<b.size()<<"->"<<c.size()<<endl;
            size_t hdeg = ((n0==0)?c.size():n0);
			size_t deg  = ((n1==0)?2*hdeg-1:n1);
            //cout<<"("<<hdeg-1<<","<<deg-1<<")"<<endl;
            if (smallLeft){
                for (size_t k=hdeg-1;k<std::min(a.size()+b.size()-1,deg);k++){
                    size_t idx_b=std::min(k,b.size()-1);
                    size_t idx_a=k-idx_b;
                    //cout<<k<<" : "<<idx_a<<"---"<<idx_b<<endl;
                    auto c_tmp=c[k-hdeg+1];
                    _BMD.mul(c_tmp,a[idx_a],b[idx_b]);
                    for (size_t j=idx_a+1;j<=std::min(k,a.size()-1);++j){
                        //cout<<k<<" : "<<j<<"---"<<k-j<<endl;
                        _BMD.axpyin(c_tmp,a[j],b[k-j]);                        
                    }
                    c.setMatrix(c_tmp,k-hdeg+1);
                }
            }else {
                for (size_t k=hdeg-1;k<std::min(a.size()+b.size()-1,deg);k++){
                    size_t idx_a=std::min(k,a.size()-1);
                    size_t idx_b=k-idx_a;
                    auto c_tmp=c[k-hdeg+1];
                    _BMD.mul(c_tmp,a[idx_a],b[idx_b]);
                    for (size_t j=idx_b+1;j<=std::min(k,b.size()-1);++j){
                        _BMD.axpyin(c_tmp,a[k-j],b[j]);
                    }
                    c.setMatrix(c_tmp,k-hdeg+1);
                }
            }
        }               
    };
	
} // end of namespace LinBox

#endif //__LINBOX_MATPOLY_MUL_NAIVE_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
