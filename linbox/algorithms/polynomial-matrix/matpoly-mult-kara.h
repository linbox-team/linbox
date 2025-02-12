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
#ifndef __LINBOX_matpoly_mult_kara_H
#define __LINBOX_matpoly_mult_kara_H

#include "linbox/algorithms/polynomial-matrix/matpoly-add-domain.h"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-naive.h"
#include "linbox/matrix/polynomial-matrix.h"
#include <algorithm>

#ifndef KARA_DEG_THRESHOLD 
#define KARA_DEG_THRESHOLD 1
#endif

namespace LinBox
{
	template<class _Field>
	class PolynomialMatrixKaraDomain {
	private:
        const _Field                            *_field;
        PolynomialMatrixAddDomain<_Field>          _PMD;
        PolynomialMatrixNaiveMulDomain<_Field>     _NMD;
        double _timeMul, _timeAdd;
	public:
        typedef _Field Field;
		typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP;
        typedef PolynomialMatrix<Field, PMType::matfirst> PMatrix;

        
        const Field & field() const { return *_field; }

        PolynomialMatrixKaraDomain(const Field &F) :
            _field(&F), _PMD(F), _NMD(F){}

        // c must be allocated with the right size
        template<typename Matrix1,typename Matrix2,typename Matrix3>
        void mul(Matrix1 &c, const Matrix2 &a, const Matrix3 &b) const
        {
        //     PMatrix a2(field(),a.rowdim(),a.coldim(),a.size());
		// 	PMatrix b2(field(),b.rowdim(),b.coldim(),b.size());
		// 	a2.copy(a,0,a.size()-1);
		// 	b2.copy(b,0,b.size()-1);
		// 	PMatrix c2(field(),c.rowdim(),c.coldim(),a.size()+b.size()-1);
		// 	mul(c2, a2, b2);
		// 	c.copy(c2,0,c2.size()-1);
        // }
		
        // void mul(PMatrix &c, const PMatrix &a, const PMatrix &b) const {
			linbox_check(c.size() >= (a.size()+b.size()-1));

            if (b.coldim()!=a.rowdim() || a.rowdim()!= a.coldim()) // ony square matrix for memory purpose
                throw ::LinBox::PreconditionFailed (__func__, __FILE__, __LINE__, "Polynomial Matrix Error -> multiplication with Karatsuba only works for square matrices\n");
            //PMatrix t(*_field,c.rowdim(),c.coldim(),std::max(a.size(),b.size()));

            PMatrix t(*_field,c.rowdim(),c.coldim(),std::max(a.size(),b.size()));
            //Matrix1 t(*_field,c.rowdim(),c.coldim(),std::max(a.size(),b.size()));
            //std::cout<<"KAra Mul TMP="<<t<<std::endl;
#ifdef KARA_TIMING
            _timeMul=_timeAdd=0.;
            Givaro::Timer chrono;
            chrono.start();
#endif
            Karatsuba_mul(c,a,b,t);
#ifdef KARA_TIMING
            cout<<" Karatsuba Multiplication:    add  : "<<_timeAdd<<" s"<<endl;
            cout<<" Karatsuba Multiplication:    mul  : "<<_timeMul<<" s"<<endl;
            cout<<" Karatsuba Multiplication:  total  : "<<chrono.userElapsedTime()<<" s"<<endl<<endl;
#endif
		}

        // c must be of size n
		// b must be of size 2n-1
		// a is of size <=n
		// compute  c= (a*b x^(-n)) mod x^n

        template<typename Matrix1,typename Matrix2,typename Matrix3>
         void midproduct(Matrix1 &c, const Matrix2 &a, const Matrix3 &b) const
         {
        //     PMatrix a2(field(),a.rowdim(),a.coldim(),a.size());
		// 	PMatrix b2(field(),b.rowdim(),b.coldim(),b.size());
		// 	a2.copy(a,0,a.size()-1);
		// 	b2.copy(b,0,b.size()-1);
		// 	PMatrix c2(field(),c.rowdim(),c.coldim(),c.size());
		// 	midproduct(c2, a2, b2);
		// 	c.copy(c2,0,c2.size()-1);
        // }
                
        // void midproduct(PMatrix &c, const PMatrix &a, const PMatrix &b) const {
			linbox_check(2*c.size()-1==b.size());
            PMatrix t(field(),c.rowdim(),c.coldim(),3*c.size());
            //Matrix1 t(field(),c.rowdim(),c.coldim(),3*c.size());
            
            // Rk: if n is a power of two, only 2n-1 extra storage is needed
            //     if n is odd, 5(n-1)/2 extra storage is needed
#ifdef KARA_TIMING
            _timeMul=_timeAdd=0.;
            Givaro::Timer chrono;
            chrono.start();
#endif
            //PMatrix a2(field(),a.rowdim(), a.coldim(), c.size());
            //a2.copy(a,0,a.size()-1);
            Karatsuba_midproduct(c,a,b,t);
#ifdef KARA_TIMING
            cout<<" Karatsuba Multiplication:    add  : "<<_timeAdd<<" s"<<endl;
            cout<<" Karatsuba Multiplication:    mul  : "<<_timeMul<<" s"<<endl;
            cout<<" Karatsuba Multiplication:  total  : "<<chrono.userElapsedTime()<<" s"<<endl<<endl;
#endif
		}


	protected:

        // PG -> does not work if A and B have different size NEED TO BE MODIFIED
		template<typename PMatrix1,typename PMatrix2,typename PMatrix3, typename PMatrix4>
        void Karatsuba_mul(PMatrix1 &C, const PMatrix2 &A, const PMatrix3& B, PMatrix4 &TMP) const {		
            //cout<<"Kar mul: "<<A.size()<<"x"<<B.size()<<"->"<<C.size()<<" ("<<TMP.size()<<")"<<endl;
            //cout<<"Kar:"<< A.rowdim()<<"x"<<A.coldim()<< "by "<<B.rowdim()<<"x"<<B.coldim()<<endl;            
            if ((A.size() <=KARA_DEG_THRESHOLD) || (B.size()<=KARA_DEG_THRESHOLD)) {
#ifdef KARA_TIMING
                Givaro::Timer chrono;
                chrono.start();
#endif
                _NMD.mul(C,A,B);
#ifdef KARA_TIMING
                _timeMul+=chrono.userElapsedTime();
#endif
            }
			else {
#ifdef KARA_TIMING
                chrono.start();
#endif
				size_t dA, dB, p, p2, q,qA, qB, n2;
                dA = A.size();
                dB = B.size();
                n2=dA+dB;
                p=std::max((dA>>1)+(dA&1),(dB>>1) +(dB&1));
                p2=p<<1;
                qA= (dA>p?dA - p:0);
                qB= (dB>p?dB - p:0);
                q=std::max(qA,qB);
                // typename PMatrix::const_view Ah,Al;
                // typename PMatrix::const_view Bh,Bl;
                // typename PMatrix::view C0,C1,C2,TMP0;
                auto Al=A.at(0,p-1); auto Ah=A.at(p,dA-1);
                auto Bl=B.at(0,p-1); auto Bh=B.at(p,dB-1);

				// add low and high terms of A in C[0,p-1]
                auto C0=C.at(0,p-1);
                _PMD.add(C0,Al,Ah);                
                
                // add low and high terms of B in C[p,2p-1]
                auto C1=C.at(p,p2-1);
                _PMD.add(C1,Bl,Bh);
#ifdef KARA_TIMING
                _timeAdd+=chrono.userElapsedTime();
#endif
                // multiply the sums in TMP[0,p2-2] using C[p2,n2-2] as temporary
                auto C2=C.at(p2,n2-2);
                auto TMP0=TMP.at(0,p2-2);
                Karatsuba_mul(TMP0,C0,C1,C2);

				// multiply high terms in C[p2,n2-2] using C[0,q-1] as temporary
                C0=C.at(0,q-1);
                Karatsuba_mul(C2, Ah,Bh, C0);
#ifdef KARA_TIMING
                chrono.start();
#endif
                // subtract the high product from the product of the sums
                _PMD.subin(TMP0,C2);

                // add TMP to the result at the right position (add when only necessary)
                for (size_t i=p;i<=p2-1;i++) C.setMatrix(TMP[i-p],i);//C[i]=TMP[i-p];
                if (p>=2){
                    C2=C.at(p2,p2+p-2);
                    TMP0=TMP.at(p,p2-2);
                    _PMD.addin(C2,TMP0);
                }
#ifdef KARA_TIMING
                _timeAdd+=chrono.userElapsedTime();
#endif

                // multiply low terms in TMP[0,2p-2] using C[2p,n2-2]] as temporary
                TMP0=TMP.at(0,p2-2);
                C0=C.at(0,p-1);
                Karatsuba_mul(TMP0,Al,Bl,C0);

#ifdef KARA_TIMING
                chrono.start();
#endif
                // subtract the low product from result at the right place
                C1=C.at(p,p+p2-2);
                _PMD.subin(C1,TMP0);

                // add TMP to the result at the right position (add when only necessary)
                for (size_t i=0;i<=p-1;i++)  C.setMatrix(TMP[i],i);//C[i]=TMP[i];
                if (p>=2){
                    C1=C.at(p,p2-2);
                    TMP0=TMP.at(p,p2-2);
                    _PMD.addin(C1,TMP0);
                }
#ifdef KARA_TIMING
                _timeAdd+=chrono.userElapsedTime();
#endif
            }
        }

        // a is of size n
		// b is of size 2n-1
		// c is of size n
		// compute  c= (a*b x^(-n)) mod x^n
        // TMP must be of size 2n
        template<typename PMatrix1,typename PMatrix2,typename PMatrix3, typename PMatrix4>
        void Karatsuba_midproduct(PMatrix1 &C, const PMatrix2 &A, const PMatrix3& B, PMatrix4 &TMP) const {

            //cout<<A.size()<<"x"<<B.size()<<"->"<<C.size()<<" ("<<TMP.size()<<")\n";

            if ((A.size() <=KARA_DEG_THRESHOLD) || (B.size()<=KARA_DEG_THRESHOLD)) {
#ifdef KARA_TIMING
                Givaro::Timer chrono;
                chrono.start();
#endif
                //cout<<"-------------------------------------------------"<<endl;
                _NMD.midproduct(C,A,B);
#ifdef KARA_TIMING
                _timeMul+=chrono.userElapsedTime();
#endif
            }
			else {
#ifdef KARA_TIMING
                chrono.start();
#endif
                //cout<<"* BEGIN ***********************************************"<<endl;
                size_t n,n0,n1,s0,s1;
                n=C.size();
                n0=n>>1; // n0 <= n1
                n1=n-n0;
                s0=2*n0-1;
                s1=2*n1-1;
                // cout<<"n="<<n<<endl;
                // cout<<"n0="<<n0<<endl;
                // cout<<"n1="<<n1<<endl;
                // cout<<"s0="<<s0<<endl;
                // cout<<"s1="<<s1<<endl;
                // cout<<"TMP SIZE:"<<TMP.size()<<endl;
                // typename PMatrix::const_view Ah,Al;
                // typename PMatrix::const_view Bh,Bm0,Bm1,Bl;
                // typename PMatrix::view C0,C1,TMP0,TMP1,TMP2;
                auto Al=A.at(0,n0-1);        // size: n0
                auto Ah=A.at(n0,n-1);        // size: n1
                auto Bl=B.at(0,s1-1);        // size: 2n1-1
                auto Bh=B.at(2*n1,2*n1+s0-1);// size: 2n0-1
                auto Bm0=B.at(n1,n1+s0-1);   // size: 2n0-1
                auto Bm1=B.at(n1,n1+s1-1);   // size: 2n1-1
                //cout<<"@@@@@@@@@@@@@@@@"<<endl;
                //cout.flush();

				// add Bl and Bm1 in C[0,s1-1]

                auto C0=C.at(0,s1-1);
                _PMD.add(C0,Bl,Bm1);
#ifdef KARA_TIMING
                _timeAdd+=chrono.userElapsedTime();
#endif
                // 1) TMP[0,n1-1]= midproduct(Ah,Bl+Bm1) using TMP[n1,n1+s1] as temporary
                auto TMP0=TMP.at(0,n1-1);
                auto TMP1=TMP.at(n1,4*n1-1);
                Karatsuba_midproduct(TMP0,Ah,C0,TMP1);
                //cout<<"alpha=:"<<TMP0<<endl;
#ifdef KARA_TIMING
                chrono.start();
#endif
                // copy TMP[0,n1-1] to C[0,n1-1]
                C0=C.at(0,n1-1);
                C0.copy(TMP0);

                if(n0!=n1) {
                    //TMP[2*n-1]=C[n1-1];
                    {
                        auto varTMP= TMP[2*n-1];
                        varTMP.copy(C[n1-1]);
                    }
                    // compute Ah-x*Al in  C[n0,n-1]
                    auto C1=C.at(n0,n-1);
                    C1.copy(Ah,0,n1-1);
                    C1=C.at(n1,n-1);
                    _PMD.subin(C1,Al);
                    C1=C.at(n0,n-1);

                    // 2) T[0,n1-1]= midproduct(Ah-Al,Bm1) using TMP[n1,n1+s1] as temporary
                    Karatsuba_midproduct(TMP0,C1,Bm1,TMP1);
                    //cout<<"beta=:"<<TMP0<<endl;
                    {
                        auto varTMP=C[n1-1];
                        varTMP.copy(TMP[2*n-1]);
                    }
                    _PMD.subin(C0,TMP0);

                    C1=C.at(n1,n-1);
                    TMP1=TMP.at(0,n0-1);
                    C1.copy(TMP1);

                    TMP0=TMP.at(0,s0-1);
                    TMP1=TMP.at(s0,s0+n0-1);
                    auto TMP2=TMP.at(s0+n0,s0+4*n0-1);

                    _PMD.add(TMP0,Bm0,Bh);

                    // 3) T[s0,s0+n0-1]= midproduct(Al,Bm0+Bh) using TMP[s0+n0,2so+n0] as temporary
                    Karatsuba_midproduct(TMP1,Al,TMP0,TMP2);
                    //cout<<"gamma=:"<<TMP1<<endl;

                    C1=C.at(n1,n-1);
                    _PMD.addin(C1,TMP1);
                }
                else {
                    // add Bh and Bm in TMP[0,s0-1]
                    TMP0=TMP.at(0,s0-1);
                    _PMD.add(TMP0,Bh,Bm0);

                    // 2) C[n1,n-1]= midproduct(Al,Bh+Bm0) using TMP[s0,2s0] as temporary
                    auto C1=C.at(n1,n-1);
                    TMP1=TMP.at(s0,s0+3*n0-1);
                    Karatsuba_midproduct(C1,Al,TMP0,TMP1);
                    //cout<<"gamma=:"<<TMP0<<endl;

                    // compute Ah-Al in  T[0,n1-1]
                    TMP0=TMP.at(0,n1-1);
                    TMP1=TMP.at(n1,2*n1-1);
                    auto TMP2=TMP.at(2*n1,2*n1+3*n0-1);

                    _PMD.sub(TMP0,Ah,Al);

                    // 3) T[n1,2n1-1]= midproduct(Ah-Al,Bm0) using TMP[2n1,2n1+s0] as temporary
                    Karatsuba_midproduct(TMP1,TMP0,Bm0,TMP2);
                    //cout<<"beta=:"<<TMP1<<endl;

                    C0=C.at(0,n1-1);
                    _PMD.subin(C0,TMP1);
                    _PMD.addin(C1,TMP1);
                }

                //cout<<"* END ***********************************************"<<endl;
            }
        }

	}; // end of class KaratsubaMulDomain<Field, Matrix>

} // end of namespace LinBox
#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
