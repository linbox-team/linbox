/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *                     Romain Lebreton
 *
 * Written by Pascal Giorgi   <pascal.giorgi@lirmm.fr>
 *            Romain Lebreton <lebreton@lirmm.fr>
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
#ifndef __LINBOX_matpoly_mult_ftt_wordsize_INL
#define __LINBOX_matpoly_mult_ftt_wordsize_INL

#include "linbox/field/modular.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/polynomial-matrix.h"


namespace LinBox {

	/**************************************************************
         * Polynomial Matrix Multiplication over Zp[x] with p <2^32 ***
         **************************************************************/
        template <>
        class PolynomialMatrixFFTMulDomain<Modular<int32_t> > {
        public:
                typedef Modular<int32_t>              Field;
                typedef Modular<integer>         LargeField;
                typedef typename Field::Element     Element;
                typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
                typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,LargeField> MatrixP_L;

        private:
                const Field            *_field;  // Read only
                uint32_t                    _p;
        public:
                inline const Field & field() const { return *_field; }

                PolynomialMatrixFFTMulDomain (const Field& F) : _field(&F), _p(F.cardinality()) {}

                template<typename Matrix1, typename Matrix2, typename Matrix3>
                void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b) {
                        uint32_t pts= 1<<(integer(a.size()+b.size()-1).bitsize());
                        if (_p< 536870912  &&  ((_p-1) % pts)==0){
                                PolynomialMatrixFFTPrimeMulDomain MulDom(field());
                                MulDom.mul(c,a,b);
                        }
                        else {  // use computation with Modular<integer>
                                // -> could be optimized in some cases (e.g. output entries less than 2^64)
                                LargeField Fp(_p);
                                PolynomialMatrixFFTMulDomain<LargeField> MulDom(Fp);
                                MatrixP_L a2(Fp,a.rowdim(),a.coldim(),a.size());
                                MatrixP_L b2(Fp,b.rowdim(),b.coldim(),b.size());
                                MatrixP_L c2(Fp,c.rowdim(),c.coldim(),c.size());
                                a2.copy(a,0,a.size()-1);
                                b2.copy(b,0,b.size()-1);
                                MulDom.mul(c2,a2,b2);
                                c.copy(c2,0,c.size()-1);
                        }
                }

                template<typename Matrix1, typename Matrix2, typename Matrix3>
                void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
                                 bool smallLeft=true, size_t n0=0,size_t n1=0) {
                        uint32_t pts= 1<<(integer(a.size()+b.size()-1).bitsize());
                        if (_p< 536870912  &&  ((_p-1) % pts)==0){
                                PolynomialMatrixFFTPrimeMulDomain MulDom(field());
                                MulDom.midproduct(c,a,b,smallLeft,n0,n1);
                        }
                        else {  // use computation with Modular<integer>
                                //-> could be optimized in some cases (e.g. output entries less than 2^64)
                                LargeField Fp(_p);
                                PolynomialMatrixFFTMulDomain<LargeField> MulDom(Fp);
                                MatrixP_L a2(Fp,a.rowdim(),a.coldim(),a.size());
                                MatrixP_L b2(Fp,b.rowdim(),b.coldim(),b.size());
                                MatrixP_L c2(Fp,c.rowdim(),c.coldim(),c.size());
                                a2.copy(a,0,a.size()-1);
                                b2.copy(b,0,b.size()-1);
                                MulDom.midproduct(c2,a2,b2,smallLeft,n0,n1);
                                c.copy(c2,0,c.size()-1);
                        }
                }


        };


}//end of namespace LinBox

#endif // __LINBOX_matpoly_mult_ftt_wordsize_INL
