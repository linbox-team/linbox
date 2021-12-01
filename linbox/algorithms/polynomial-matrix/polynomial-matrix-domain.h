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


#ifndef __LINBOX_POLYNOMIAL_MATRIX_DOMAIN_H
#define __LINBOX_POLYNOMIAL_MATRIX_DOMAIN_H

#define KARA_DEG_THRESHOLD  2 
#define FFT_DEG_THRESHOLD   2

#include "linbox/algorithms/polynomial-matrix/matpoly-mult-naive.h"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-kara.h"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft.h"
#include <algorithm>


        
namespace LinBox
{     
	template <class _Field> 
	class PolynomialMatrixMulDomain {
	public:
		PolynomialMatrixKaraDomain<_Field>       _kara;
		PolynomialMatrixFFTMulDomain<_Field>      _fft;
		PolynomialMatrixNaiveMulDomain<_Field>  _naive;
		const _Field*                           _field;
	public:
        typedef _Field Field; 
		PolynomialMatrixMulDomain (const Field &F) :
			_kara(F), _fft(F), _naive(F), _field(&F) {}

		inline const Field& field() const {return *_field;}

		template< class PMatrix1,class PMatrix2,class PMatrix3>
		void mul(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b, size_t max_rowdeg=0) const
		{
			size_t d = a.size()+b.size();
            if (d > FFT_DEG_THRESHOLD){
                    //std::cout<<"PolMul FFT"<<std::endl;
				_fft.mul(c,a,b);
            }
			else
				if ( d > KARA_DEG_THRESHOLD){
                        //std::cout<<"PolMul Kara"<<std::endl;
					_kara.mul(c,a,b);
                }
				else {
                        //std::cout<<"PolMul Naive"<<std::endl;
					_naive.mul(c,a,b);
                }
#ifdef CHECK_MATPOL_MUL                       
            check_mul(c,a,b,c.size());
#endif
        }

		template< class PMatrix1,class PMatrix2,class PMatrix3>
		void midproduct (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b) const
		{
			size_t d = b.size();
			if (d > FFT_DEG_THRESHOLD)
				_fft.midproduct(c,a,b);
			else
				if ( d > KARA_DEG_THRESHOLD)
					_kara.midproduct(c,a,b);
				else
					_naive.midproduct(c,a,b);
#ifdef CHECK_MATPOL_MIDP                       
                        check_midproduct(c,a,b);
#endif

		}

		template< class PMatrix1,class PMatrix2,class PMatrix3>
		void midproductgen (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b, bool smallLeft=true, size_t n0=0, size_t n1=0) const
		{
			if ( c.size() <= 4)
				_naive.midproduct(c,a,b,smallLeft,n0,n1);
			else
				_fft.midproduct(c,a,b,smallLeft,n0,n1);
#ifdef CHECK_MATPOL_MIDP                       
                        check_midproduct(c,a,b,smallLeft,n0,n1);
#endif               
		}

	};

	template<class Field>
	class PolynomialMatrixDomain : public PolynomialMatrixMulDomain<Field>, public PolynomialMatrixAddDomain<Field> {
	public:
		PolynomialMatrixDomain (const Field& F) : PolynomialMatrixMulDomain<Field>(F), PolynomialMatrixAddDomain<Field>(F) {}
	};




	template <typename Field, typename PMatrix1,typename PMatrix2,typename PMatrix3>
	class HalflineMPDomain {
	private:
		PMatrix1                        *_c;
		const PMatrix2                  *_a;
		const PMatrix3                  *_b;
		typename PMatrix3::const_view   _b1;
		typename PMatrix1::plain       _tmp;
		size_t                           _i;
		size_t                           _d;
		PolynomialMatrixDomain<Field>  _PMD;
		BlasMatrixDomain<Field>        _BMD;
	public:
		double mid;
		const Field& field(){return _PMD.field();}

		HalflineMPDomain(const Field& F):  _i(1), _d(0), _PMD(F), _BMD(F) {}

		// compute c = x^k ab mod x^2k
		// a is of sike k+1
		// b is of size 2k
		// c is of size k
		HalflineMPDomain(const Field& F, PMatrix1& c, const PMatrix2& a, const PMatrix3 &b,
				 size_t k)
                        : _c(&c), _a(&a), _b(&b), _b1(b,k,std::min(b.size()-1,2*k-1)), _tmp(F,c.rowdim(),c.coldim(),2*k-1),
                          _i(1), _d(std::min(b.size()-1,2*k-1)-k+1), _PMD(F), _BMD(F), mid(0)
		{
			linbox_check(c.size()==k);
			linbox_check(a.size()<=k+1);
			//linbox_check(b.size()==2*k);
			Timer chrono;
			chrono.start();
			// compute the high product of c= x^(k-1) a[1,k]b[0,k-1] mod x^2k
			typename PMatrix2::const_view a0= a.at(1,a.size()-1);
			typename PMatrix3::const_view b0= b.at(0,k-1);
			_PMD.mul(_tmp,a0,b0);
			c.copy(_tmp,k-1,2*k-2);
			chrono.stop();
			mid+=chrono.usertime();

		}

		void update(size_t s=1){
			for (size_t i=0;i<s;i++){
				++(*this);
			}
		}

		void operator++() {
			if (!terminated()){ // compute product at step _i
				//cout<<"read coeff <="<<_d+_i-1<<" of "<< 2*_d<<" ..."<<endl;;
				size_t m = twoValuation(_i);
				size_t step=  1ULL<<m;
				//typename PMatrix3::const_view bb = _b->at(0, 2*_d-1);

				if (step<_d) {
					// cout<<"("<<    0    <<","<< min(2*step-2,_a->size()-1)<<") x";
					// cout<<"("<< _i-step <<","<<      _i-1 <<") ->";
					// cout<<"("<< _i-1    <<","<< _i+step-2 <<") "<<endl;;

					typename PMatrix1::view     tmp = _tmp.at(0  , step-1);
					typename PMatrix1::view       c = _c->at(_i-1, _i+step-2);
					typename PMatrix2::const_view a = _a->at(0, std::min(2*step-2,_a->size()-1));
					typename PMatrix3::const_view b = _b1.at(_i-step, _i-1);
					Timer chrono;
					chrono.start();
					_PMD.midproductgen(tmp,a,b,false);//,step,2*step);
					chrono.stop();
					_PMD.addin(c,tmp);
					mid+=chrono.usertime();
				}
				else {
					// compute the last diagonal element
					for(size_t i=0;i<std::min(_d,_a->size());i++)
						_BMD.axpyin((*_c)[_d-1],(*_a)[i],_b1[_d-i-1]);
					/*
					   cout<<"checking result of online MP ..."<<_d<<"x"<<2*_d<<"...";
					   PMatrix1 R(field(),_c->rowdim(),_c->coldim(),_c->size());
					   typename PMatrix3::const_view b = _b->at(0, 2*_d-1);
					   _PMD.midproductgen(R,*_a, b, true, _d+1,2*_d);
					   if (R==(*_c))
					   cout<<"done"<<endl;
					   else
					   cout<<"error"<<endl;
					   */
				}
				_i++;
			}
		}

		bool terminated() const {return _i>_d;}

		inline size_t twoValuation(size_t x){
			size_t i=0;
			while (x!=0 && !(x&0x1)){
				i++;
				x>>=1;
			}
			return i;
		}


	};

} // end of namespace LinBox

#endif //__LINBOX_matpoly_mult_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
