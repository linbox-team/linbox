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

#include "linbox/ring/modular.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-three-primes.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-fast.inl"
namespace LinBox {

	/**************************************************************
     * Polynomial Matrix Multiplication over Zp[x] with p <2^32 ***
     **************************************************************/
    template <class T1, class T2>
    class PolynomialMatrixFFTMulDomain<Givaro::Modular<T1,T2> > {
    public:
        typedef Givaro::Modular<T1,T2>                Field;
        typedef Givaro::Modular<integer>         LargeField;
        typedef typename Field::Element             Element;
        typedef PolynomialMatrix<Field,PMType::polfirst>         MatrixP;
        typedef PolynomialMatrix<LargeField,PMType::polfirst>  MatrixP_L;

    private:
        const Field            *_field;  // Read only
        uint64_t                    _p;
    public:
        inline const Field & field() const { return *_field; }

        PolynomialMatrixFFTMulDomain (const Field& F) : _field(&F), _p(F.cardinality()) {}

        template<typename Matrix1, typename Matrix2, typename Matrix3>
        void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, size_t max_rowdeg=0) const {

            size_t deg  = (max_rowdeg?max_rowdeg:a.size()+b.size()-2); //size_t deg  = a.size()+b.size()-1;
            c.resize(deg+1);
            size_t pts  = 1; while (pts <= deg) { pts<<=1; }
            if ( _p< 536870912ULL  &&  ((_p-1) % pts)==0){
                PolynomialMatrixFFTPrimeMulDomain<Field> MulDom(field());
                MulDom.mul(c,a,b, max_rowdeg);
            }
            else {
                if (_p< 536870912ULL){
                    PolynomialMatrixThreePrimesFFTMulDomain<Field> MulDom(field());
                    MulDom.mul(c,a,b, max_rowdeg);
                }
				else {
					// use computation with Givaro::Modular<integer>
					// -> could be optimized in some cases (e.g. output entries less than 2^64)
					FFT_PROFILE_START(2);
					LargeField Fp(_p);
					PolynomialMatrixFFTMulDomain<LargeField> MulDom(Fp);
					MatrixP_L a2(Fp,a.rowdim(),a.coldim(),a.size());
					MatrixP_L b2(Fp,b.rowdim(),b.coldim(),b.size());
					MatrixP_L c2(Fp,c.rowdim(),c.coldim(),c.size());
					a2.copy(a,0,a.degree());
					b2.copy(b,0,b.degree());
					FFT_PROFILING(2,"converting rep of polynomial matrix input");
					MulDom.mul(c2,a2,b2, max_rowdeg);
					c.copy(c2,0,c.degree());
					FFT_PROFILING(2,"converting rep of polynomial matrix output");
				}
            }
        }

        template<typename Matrix1, typename Matrix2, typename Matrix3>
        void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
                         bool smallLeft=true, size_t n0=0,size_t n1=0) const {
            uint64_t pts= 1<<(integer((uint64_t)a.size()+b.size()-1).bitsize());
            if (_p< 536870912ULL  &&  ((_p-1) % pts)==0){
				//std::cout<<"MIDP: Staying with FFT Prime Field"<<std::endl;
                PolynomialMatrixFFTPrimeMulDomain<Field> MulDom(field());
                MulDom.midproduct(c,a,b,smallLeft,n0,n1);
            }
			else {
				if (_p< 536870912ULL){
					PolynomialMatrixThreePrimesFFTMulDomain<Field> MulDom(field());
					MulDom.midproduct(c,a,b,smallLeft,n0,n1);
				}
				else {  // use computation with Givaro::Modular<integer>
					//-> could be optimized in some cases (e.g. output entries less than 2^64)
					FFT_PROFILE_START(2);
					//std::cout<<"MIDP: Switching to Large Field"<<std::endl;
					LargeField Fp(_p);
					PolynomialMatrixFFTMulDomain<LargeField> MulDom(Fp);
					MatrixP_L a2(Fp,a.rowdim(),a.coldim(),a.size());
					MatrixP_L b2(Fp,b.rowdim(),b.coldim(),b.size());
					MatrixP_L c2(Fp,c.rowdim(),c.coldim(),c.size());
					a2.copy(a,0,a.size()-1);
					b2.copy(b,0,b.size()-1);
					FFT_PROFILING(2,"converting rep of polynomial matrix input");
					MulDom.midproduct(c2,a2,b2,smallLeft,n0,n1);
					c.copy(c2,0,c.size()-1);
					FFT_PROFILING(2,"converting rep of polynomial matrix output");
				}
			}
		}


    };


}//end of namespace LinBox

#endif // __LINBOX_matpoly_mult_ftt_wordsize_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
