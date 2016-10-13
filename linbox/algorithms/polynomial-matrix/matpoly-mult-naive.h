/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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

#include "linbox/util/error.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/polynomial-matrix.h"
#include <algorithm>
namespace LinBox
{
        
	template <class Field>
	class PolynomialMatrixNaiveMulDomain {
	private:
		const Field            *_field;
		BlasMatrixDomain<Field>   _BMD;

		// Polynomial matrix stored as a matrix of polynomial
		typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
		// Polynomial matrix stored as a polynomial of matrix
		typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;
                template <PMStorage S> using _PMatrix = PolynomialMatrix<PMType::matfirst,S,Field>;
	public:

                inline const Field& field() const { return *_field; }

                PolynomialMatrixNaiveMulDomain(const Field &F) :
			_field(&F), _BMD(F){}
                

                void mul(PMatrix&  c, const PMatrix&  a, const PMatrix&  b)
                {
                        for (size_t k=0;k<a.size()+b.size()-1;k++){
                                size_t idx_min= (k+1<b.size()?0:k+1-b.size());
                                size_t idx_max=std::min(k,a.size()-1);                                
                                _BMD.mul(c[k],a[idx_min],b[k-idx_min]);
                                for (size_t i=idx_min+1;i<=idx_max;i++){ 
                                        _BMD.axpyin(c[k],a[i],b[k-i]);
                                }
                        }
                }                

                // c must be allocated with the right size
                template<typename Matrix1,typename Matrix2,typename Matrix3>
                void mul(Matrix1& c, const Matrix2&  a, const Matrix3&  b)
                {
                        PMatrix a2(field(),a.rowdim(),a.coldim(),a.size());
			PMatrix b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			PMatrix c2(field(),c.rowdim(),c.coldim(),a.size()+b.size()-1);
			mul(c2, a2, b2);
			c.copy(c2,0,c2.size()-1);
                }                          


                
                // c must be allocated with the right size 
                // a and b can have a size smaller than required
                template<typename Matrix1,typename Matrix2,typename Matrix3>
                void midproduct(Matrix1& c, const Matrix2&  a, const Matrix3&  b,
                                bool smallLeft=true, size_t n0=0,size_t n1=0) 
                {
                        PMatrix a2(field(),a.rowdim(),a.coldim(),a.size());
			PMatrix b2(field(),b.rowdim(),b.coldim(),b.size());
			a2.copy(a,0,a.size()-1);
			b2.copy(b,0,b.size()-1);
			PMatrix c2(field(),c.rowdim(),c.coldim(),c.size());
			midproduct(c2, a2, b2, smallLeft, n0, n1);
			c.copy(c2,0,c2.size()-1);
                }                          
                void midproduct(PMatrix&  c, const PMatrix&  a, const PMatrix&  b,
                                bool smallLeft=true, size_t n0=0,size_t n1=0) {
                        //cout<<"naive midprod "<<a.size()<<"x"<<b.size()<<"->"<<c.size()<<endl;
                        size_t hdeg = ((n0==0)?c.size():n0);
			size_t deg  = ((n1==0)?2*hdeg-1:n1);
                        //cout<<"("<<hdeg-1<<","<<deg-1<<")"<<endl;
                        if (smallLeft){
                                // for (size_t k=hdeg-1;k<std::min(a.size()+b.size()-1,deg);k++){
                                //         _BMD.mul(c[k-hdeg+1],a[0],b[k]);
                                //         for (size_t j=1;j<std::min(hdeg,a.size());++j)
                                //                 _BMD.axpyin(c[k-hdeg+1],a[j],b[k-j]);
                                //}
                                for (size_t k=hdeg-1;k<std::min(a.size()+b.size()-1,deg);k++){
                                        size_t idx_b=std::min(k,b.size()-1);
                                        size_t idx_a=k-idx_b;
                                        //cout<<k<<" : "<<idx_a<<"---"<<idx_b<<endl;
                                        _BMD.mul(c[k-hdeg+1],a[idx_a],b[idx_b]);
                                        for (size_t j=idx_a+1;j<=std::min(k,a.size()-1);++j){
                                                //cout<<k<<" : "<<j<<"---"<<k-j<<endl;
                                                _BMD.axpyin(c[k-hdeg+1],a[j],b[k-j]);
                                        }
                                }
                        }else {
                                for (size_t k=hdeg-1;k<std::min(a.size()+b.size()-1,deg);k++){
                                        size_t idx_a=std::min(k,a.size()-1);
                                        size_t idx_b=k-idx_a;
                                        _BMD.mul(c[k-hdeg+1],a[idx_a],b[idx_b]);
                                        for (size_t j=idx_b+1;j<=std::min(k,b.size()-1);++j){
                                                _BMD.axpyin(c[k-hdeg+1],a[k-j],b[j]);
                                        }
                                }
                        }
                }               
        };
	
} // end of namespace LinBox

#endif //__LINBOX_MATPOLY_MUL_NAIVE_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

