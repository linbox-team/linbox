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

#ifndef __LINBOX_matpoly_mult_ftt_H
#define __LINBOX_matpoly_mult_ftt_H

#include "linbox/util/error.h"
#include "linbox/util/debug.h"
#include "linbox/util/timer.h"
#include "linbox/util/commentator.h"
#include <linbox/randiter/random-fftprime.h>
#include <linbox/randiter/random-prime.h>
#include "linbox/integer.h"
#include <givaro/zring.h>
#include "linbox/ring/modular.h"
#include "givaro/givtimer.h"
#include <sstream>
#include <iostream>

#ifdef FFT_PROFILER
#ifndef FFT_PROF_LEVEL
int  FFT_PROF_LEVEL=1;
#endif
Givaro::Timer mychrono[3];
#define FFT_PROF_MSG_SIZE 35
#define FFT_PROFILE_START(lvl)  mychrono[lvl].clear();mychrono[lvl].start();

#define FFT_PROFILING(lvl,msg)                                          \
    if (lvl>=FFT_PROF_LEVEL) {                                          \
        mychrono[lvl].stop();std::cout<<"FFT("<<lvl<<"):";              \
        std::cout.width(FFT_PROF_MSG_SIZE);std::cout<<std::left<<msg<<" : "; \
        std::cout.precision(6);std::cout<<mychrono[lvl]<<std::endl;		\
        mychrono[lvl].clear();mychrono[lvl].start();                    \
    }
  
#ifdef HAVE_OPENMP								
#define FFT_PROFILE_GET(lvl,x)                                          \
    // mychrono.stop();(x)+=mychrono.realtime();mychrono.clear();mychrono.start();
#else
#define FFT_PROFILE_GET(lvl,x)                                          \
    mychrono[lvl].stop();(x)+=mychrono[lvl].usertime();mychrono[lvl].clear();mychrono[lvl].start();
#endif
#define FFT_PROFILE(lvl,msg,x)                                          \
    if ((lvl)>=FFT_PROF_LEVEL) {                                        \
        std::cout<<"FFT: ";                                             \
        std::cout.width(FFT_PROF_MSG_SIZE);std::cout<<std::left<<msg<<" : "; \
        std::cout.precision(6);std::cout<<x<<" s"<<std::endl;           \
    }
#else
#define FFT_PROFILE_START(lvl)
#define FFT_PROFILING(lvl,msg)
#define FFT_PROFILE_GET(lv,x)
#define FFT_PROFILE(lvl,msg,x)
#endif // FFT_PROFILER


#ifndef FFT_DEG_THRESHOLD   
#define FFT_DEG_THRESHOLD   4
#endif

namespace LinBox
{
    template<typename Field>
    bool check_mul (const PolynomialMatrix<Field, PMType::matfirst> &c,
                    const PolynomialMatrix<Field, PMType::matfirst> &a,
                    const PolynomialMatrix<Field, PMType::matfirst> &b,size_t deg) {

        typedef typename PolynomialMatrix<Field, PMType::matfirst>::Matrix Matrix;
        Matrix C1(c.field(),c.rowdim(),c.coldim()),C2(c.field(),c.rowdim(),c.coldim());
        Matrix A1(c.field(),a.rowdim(),a.coldim());
        Matrix B1(c.field(),b.rowdim(),b.coldim());
        BlasMatrixDomain< typename Matrix::Field>  BMD(c.field());
        for (size_t k=0;k<deg;k++)
            BMD.addin(C1,c[k]);
        for (size_t k=0;k<a.size();k++)
            BMD.addin(A1,a[k]);
        for (size_t k=0;k<b.size();k++)
            BMD.addin(B1,b[k]);
    
        BMD.mul(C2,A1,B1);
        bool correct=BMD.areEqual(C1,C2);
        std::ostream& report = LinBox::commentator().report();
        report<<"Checking polynomial matrix mul "
              <<a.rowdim()<<"x"<<a.coldim()<<"["<<a.size()<<"]"
              <<b.rowdim()<<"x"<<b.coldim()<<"["<<b.size()<<"]"
              <<" ... "<<(correct?"done":"error")<<std::endl;
        if (!correct){
            report<<"error with field : ";
            c.field().write(report);
            report<<std::endl;
            // std::cerr<<"A:="<<a<<";"<<std::endl;
            // std::cerr<<"B:="<<b<<";"<<std::endl;
            // std::cerr<<"C:="<<c<<";"<<std::endl;
        }
        return correct;
    }
 
    template<typename Field>
    bool check_mul (const PolynomialMatrix<Field, PMType::polfirst> &c,
                    const PolynomialMatrix<Field, PMType::polfirst> &a,
                    const PolynomialMatrix<Field, PMType::polfirst> &b,size_t deg) {

        typedef typename PolynomialMatrix<Field, PMType::polfirst>::Matrix Matrix;
        Matrix C1(c.field(),c.rowdim(),c.coldim()),C2(c.field(),c.rowdim(),c.coldim());
        Matrix A1(c.field(),a.rowdim(),a.coldim());
        Matrix B1(c.field(),b.rowdim(),b.coldim());
        BlasMatrixDomain< typename Matrix::Field>  BMD(c.field());

        // C1 = c(1)
        for (size_t i=0;i<c.rowdim()*c.coldim();i++)
            for (size_t k=0;k<deg;k++)
                c.field().addin(C1.getPointer()[i], c.get(i,k));

        // A1=a(1)
        for (size_t i=0;i<a.rowdim()*a.coldim();i++)
            for (size_t k=0;k<a.size();k++)
                a.field().addin(A1.getPointer()[i], a.get(i,k));

        // B1=b(1)
        for (size_t i=0;i<b.rowdim()*b.coldim();i++)
            for (size_t k=0;k<b.size();k++)
                b.field().addin(B1.getPointer()[i], b.get(i,k));
    
        BMD.mul(C2,A1,B1);
        bool correct=BMD.areEqual(C1,C2);
        std::ostream& report = LinBox::commentator().report();
        report<<"Checking polynomial matrix mul "
              <<a.rowdim()<<"x"<<a.coldim()<<"["<<a.size()<<"]"
              <<b.rowdim()<<"x"<<b.coldim()<<"["<<b.size()<<"]"
              <<" ... "<<(correct?"done":"error")<<std::endl;
        if (!correct){
            report<<"error with field : ";
            c.field().write(report);
            report<<std::endl;
            report<<"A:="<<a<<";"<<std::endl;
            report<<"B:="<<b<<";"<<std::endl;
            report<<"C:="<<c<<";"<<std::endl;
        }
        return correct;
    }
  


    template<typename MatrixP_F>
    bool check_midproduct (const MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b, bool smallLeft=true, size_t n0=0,size_t n1=0, size_t deg=0) {
        typename MatrixP_F::Matrix C1(c.field(),c.rowdim(),c.coldim()),C2(c.field(),c.rowdim(),c.coldim());
        typename MatrixP_F::Matrix A1(c.field(),a.rowdim(),a.coldim());
        typename MatrixP_F::Matrix B1(c.field(),b.rowdim(),b.coldim());
        BlasMatrixDomain< typename MatrixP_F::Field>  BMD(c.field());

        std::ostream& report = LinBox::commentator().report();
        if (deg==0) deg=c.size();
        if (n0 == 0) n0=deg;
        if (n1 == 0) n1=2*deg-1;

        std::ostringstream myerror;
        myerror<<"A size: "<<a.size()<<std::endl;
        myerror<<"B size: "<<b.size()<<std::endl;
        myerror<<"deg"<<deg<<std::endl;
        myerror<<"n0="<<n0<<std::endl;
        myerror<<"n1="<<n1<<std::endl;
    
        for (size_t k=0;k<deg;k++)
            BMD.addin(C1,c[k]);
    

        if (smallLeft){
            size_t k=std::min(n0-1,a.size()-1);
            size_t j=n0-1-k;
            size_t t=0;
            for (;k<size_t(-1) && j<b.size() && t<deg;k--,j++,t++){
                //myerror<<"+a["<<k<<"]"<<std::endl;
                BMD.addin(A1,a[k]);
                //myerror<<"*b["<<j<<"]"<<std::endl;
                BMD.axpyin(C2,A1,b[j]);
            }
      
            for(;t<deg;t++,j++){
                //myerror<<"*b["<<j<<"]"<<std::endl;
                BMD.axpyin(C2,A1,b[j]);		
            }
	
            /* for (;j>=lastj ;j--){ */
            /* 	myerror<<"*b["<<j<<"]"<<std::endl; */
            /* 	BMD.axpyin(C2,A1,b[j]); */
            /* } */


            //myerror<<"-------------\n";
            size_t lastj=j;
            size_t lastk=k;
            j=std::min(n1-1,b.size()-1);
            k=n1-1-j;
            t=0;
            //myerror<<"lastj="<<lastj<<std::endl;
            //myerror<<"lastk="<<lastk<<std::endl;
            //myerror<<"j="<<j<<std::endl;
            //myerror<<"k="<<k<<std::endl;
            A1.zero();

            while(j>=lastj && n0==n1){
                BMD.axpyin(C2,a[k++],b[j--]);
            }
      
      
            for (;j>=lastj && k<a.size() && t<deg;k++,j--,t++){
                //myerror<<"+a["<<k<<"]"<<std::endl;
                BMD.addin(A1,a[k]);
                //myerror<<"*b["<<j<<"]"<<std::endl;
                BMD.axpyin(C2,A1,b[j]);
            }
      
            if (lastk>0 && lastk < size_t(-1) && lastj <= j){
                A1.zero();
                for(size_t t=0;t<deg;t++)
                    BMD.addin(A1,a[lastk+t]);
                BMD.axpyin(C2,A1,b[lastj]);
            }

      
        }

        else {
            report<<"Checking polynomial matrix midp  with smallLeft=false is not yet implemented...aborting";
            std::terminate();
      
        }
    
        bool correct=BMD.areEqual(C1,C2);
        report<<"Checking polynomial matrix "<<(n0==0&&n1==0?"midp ":"midp_gen ")
              <<a.rowdim()<<"x"<<a.coldim()<<"["<<a.size()<<"]"
              <<b.rowdim()<<"x"<<b.coldim()<<"["<<b.size()<<"]"
              <<" ... "<<(correct?"done":"error")<<std::endl;
        if (!correct){
            myerror<<"error with field : ";
            c.field().write(myerror);
            myerror<<std::endl;
            if (c.size()*c.rowdim() < 120) {
                myerror<<"A:="<<a<<";"<<std::endl;
                myerror<<"B:="<<b<<";"<<std::endl;
                myerror<<"C:="<<c<<";"<<std::endl;
                myerror<<"C1:=";C1.write(myerror,Tag::FileFormat::Maple)<<std::endl;
                myerror<<"C2:=";C2.write(myerror,Tag::FileFormat::Maple)<<std::endl;

            }
            report<<myerror.str()<<std::endl;
            //std::terminate();
        }
        return correct;
    }
    
    




  
    // generic handler for multiplication using FFT
    template <class Field>
    class PolynomialMatrixFFTMulDomain {
    public:
        inline const Field & field() const;

        PolynomialMatrixFFTMulDomain (const Field& F);

        template<typename Matrix1, typename Matrix2, typename Matrix3>
        void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b) const;

        template<typename Matrix1, typename Matrix2, typename Matrix3>
        void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, bool smallLeft=true, size_t n0=0,size_t n1=0) const;
    };
		
	
    //class PolynomialMatrixFFTPrimeMulDomain ;                         // Mul in Zp[x] with p <2^32, (fflas, fourier)
		
    // template <class T>
    // class PolynomialMatrixFFTMulDomain<Givaro::Modular<T> > ;        // Mul in Zp[x] with p^2 storable in type T

    // template<>
    // class PolynomialMatrixFFTMulDomain<Givaro::ZRing<integer> >;  // Mul in Z[x]

    // template <>
    // class PolynomialMatrixFFTMulDomain<Givaro::Modular<integer> > ;           // Mul in Zp[x] with p multiprecision

    // get the maximum prime for fft with modular<double> (matrix dim =k, nbr point = pts)
    uint64_t maxFFTPrimeValue(uint64_t k, uint64_t pts) {
        uint64_t prime_max=std::sqrt( (1ULL<<53) /k)+1;
        size_t c=1;
        const int fct=24;
        while (c<k && prime_max < (1UL<<26) && prime_max< pts*fct){
            prime_max=std::sqrt( (1ULL<<53) /(k/c))+1;
            c<<=1;
        }

        //std::cout<<"maxFFTPrime: pts -> "<<pts<<std::endl;
        //std::cout<<"maxFFTPrime: replacing "<<k<<" -> "<<k/c<<std::endl;
	  
        if (c>=k && c!=1){
            std::cout<<"MatPoly FFT (maxPrimeValue): impossible to find enough FFT Prime\n";
            std::terminate();
        }
	  
        return std::min(prime_max, uint64_t(Givaro::Modular<double>::maxCardinality()));
    }

    void getFFTPrime(uint64_t prime_max, size_t lpts, integer bound, std::vector<integer> &bas, size_t k, size_t d){
        size_t nbp=0;
        bool b = RandomFFTPrime::generatePrimes (bas, prime_max, bound, lpts);
        if (!b){ /* not enough FFT prime found */
            integer MM=1;
            for(std::vector<integer>::size_type i=0;i<bas.size();i++){
                MM*=bas[i];
                //std::cout<<bas[i]<<std::endl;
            }
	    
            // compute max bitsize for prime allowing three prime fft
            integer prime_max_tp=MM/uint64_t(d*k);
            while (k>1 && prime_max_tp<100) {k/=2;prime_max_tp*=2;}
            if (k<=1) {std::cout<<"getFFTPrime error: impossible to have enough primes satisfying constraints: FFLAS prime (<2^26) and FFT (2^"<<lpts<<")\n";}
	
            PrimeIterator<IteratorCategories::HeuristicTag> Rd(std::min(prime_max_tp.bitsize()/2,integer(prime_max).bitsize())-1);
#ifdef VERBOSE_FFT
            std::cout<<"MM="<<MM<<std::endl;
            std::cout<<"normal primemax: "<<prime_max_tp<<" "<<prime_max<<std::endl;
            std::cout<<"normal prime bitmax: "<<std::min(prime_max_tp.bitsize()/2,integer(prime_max).bitsize()-1)<<std::endl;
#endif
            integer tmp;
            do {
                do {tmp = *(++Rd);}
                while (MM%tmp==0 || tmp>prime_max);
                bas.push_back(tmp);
                nbp++;
                MM*=tmp;
            } while (MM<bound);	
        }
#ifdef VERBOSE_FFT      
        std::cout<<"MatPoly Multiprecision FFT : using "<<bas.size()-nbp<<" FFT primes and "<<nbp<<" normal primes "<<std::endl;
#endif
        for(auto i: bas)
            if (i>prime_max) std::cout<<"ERROR\n";
    }

	
} // end of namespace LinBox

#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-fast.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize-three-primes.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-multiprecision.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-recint.inl"
#include "linbox/algorithms/polynomial-matrix/matpoly-mult-fft-wordsize.inl"

#endif // __LINBOX_matpoly_mult_ftt_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
