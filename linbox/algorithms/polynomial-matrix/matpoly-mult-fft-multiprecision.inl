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
#ifndef __LINBOX_matpoly_mult_ftt_multiprecision_INL
#define __LINBOX_matpoly_mult_ftt_multiprecision_INL

#include <givaro/zring.h>
#include "linbox/ring/modular.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include <fflas-ffpack/field/rns-double.h>
#define MB(x) ((x)/(double)(1<<20))
#ifndef MEMINFO
#define MEMINFO ""
#endif
namespace LinBox{

  /***************************************************
   **** Polynomial Matrix Multiplication over Z[x] ***
   ***************************************************/
  template<>
  class PolynomialMatrixFFTMulDomain<Givaro::ZRing<integer> > {
  public:
    typedef Givaro::ZRing<integer>       IntField;
    //typedef Givaro::Modular<uint32_t>     ModField;
      typedef Givaro::Modular<double>                ModField;
      typedef PolynomialMatrix<ModField, PMType::polfirst> MatrixP_F; // Polynomial matrix stored as a matrix of polynomials
      typedef PolynomialMatrix<IntField, PMType::polfirst> MatrixP_I; // Polynomial matrix stored as a matrix of polynomials

  private:
    const IntField     *_field;
    integer           _maxnorm;

    template<typename PMatrix1>
    size_t logmax(const PMatrix1& A) const {
      size_t mm=A.get(0,0,0).bitsize();
      for(size_t k=0;k<A.size();k++)
	for (size_t i=0;i<A.rowdim()*A.coldim();i++){
	  size_t tmp=A.get(i,k).bitsize();
	  mm=std::max(mm,tmp);
	}
      return mm;
    }

  public: 

    inline const IntField & field() const { return *_field; }


    PolynomialMatrixFFTMulDomain (const IntField &F, const integer maxnorm=0) :
      _field(&F), _maxnorm(maxnorm) {}

    template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
    void mul (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b, size_t max_rowdeg=0) const {
      //compute a bound on the entry of the input matrix a and b
      FFT_PROFILE_START(2);
      integer maxA,maxB;
      maxA=maxB=_maxnorm;			
      if (_maxnorm==0){
	maxA=1;maxA<<=uint64_t(logmax(a));
	maxB=1;maxB<<=uint64_t(logmax(b));
      }
      integer bound=maxA*maxB*uint64_t(a.coldim())*uint64_t(std::min(a.size(),b.size()));
      if (_maxnorm==0) bound*=2; //seems to compute over Z, need to double to handle possible negative value
      FFT_PROFILING(2,"max norm computation");

      mul_crtla(c,a,b,maxA,maxB,bound,max_rowdeg);
    }

    template<typename PMatrix1, typename PMatrix2, typename PMatrix3>
    void midproduct (PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
		     bool smallLeft=true, size_t n0=0, size_t n1=0) const {
      //compute a bound on the entry of the input matrix a and b
      FFT_PROFILE_START(2);
      integer maxA,maxB;
      maxA=maxB=_maxnorm;
      if (_maxnorm==0){
          maxA=1;maxA<<=uint64_t(logmax(a));
          maxB=1;maxB<<=uint64_t(logmax(b));
      }
      //std::cout<<"MIDP RNS bound: "<<maxA<<" "<<maxB<<" "<<a.coldim()<<" "<<a.size()<<" "<<b.size()<<std::endl;
      
      //integer bound=2*maxA*maxB*integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));;
      integer bound=maxA*maxB*integer((uint64_t)a.coldim());
      if (_maxnorm==0) bound*=2; //seems to compute over Z, need to double to handle possible negative value
      
      if (smallLeft)
	bound*= (uint64_t)a.size();
      else
	bound*= (uint64_t)b.size();

      FFT_PROFILING(2,"max norm computation");

      midproduct_crtla(c,a,b,maxA,maxB,bound,smallLeft, n0,n1);
    }


    // template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
    // void mul_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
    // 		   const integer& maxA, const integer& maxB, const integer& bound) {
    //   // (convert to MatrixP representation)
    //   FFT_PROFILE_START;
    //   MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
    //   MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
    //   a2.copy(a,0,a.size()-1);
    //   b2.copy(b,0,b.size()-1);
    //   MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
    //   FFT_PROFILING(2,"converting rep of input matrices");
    //   mul_crtla(c2,a2,b2,maxA,maxB,bound);
    //   c.copy(c2,0,c.size()-1);
    //   FFT_PROFILING(2,"converting rep of output matrices");
    // }


    // void mul_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b,
    // 		   const integer& maxA, const integer& maxB, const integer& bound){

    // WARNING: Polynomial Matrix should stored as matrix of polynomial with integer coefficient 
    template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
    void mul_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
    		   const integer& maxA, const integer& maxB, const integer& bound, size_t max_rowdeg=0) const {

      FFT_PROFILE_START(2);
      linbox_check(a.coldim() == b.rowdim());
      size_t m = a.rowdim();
      size_t k = a.coldim();
      size_t n = b.coldim();
      size_t s= a.size()+b.size()-1;
      if (max_rowdeg!=0) s = max_rowdeg+1;
      c.resize(s);
      size_t lpts=0;
      size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }

      // compute bit size of feasible prime for FFLAS
      // size_t _k=k,lk=0;
      //while ( _k ) {_k>>=1; ++lk;}
      //size_t prime_bitsize= (53-lk)>>1;

      // compute max prime value for FFLAS      
      uint64_t prime_max= std::min(uint64_t(std::sqrt( (1ULL<<53) / k)+1), uint64_t(Givaro::Modular<double>::maxCardinality()));
      std::vector<integer> bas;
      getFFTPrime(prime_max,lpts,bound,bas,k,s);
      // RandomFFTPrime RdFFT(prime_bitsize);
      // if (!RdFFT.generatePrimes(lpts,bound,bas)){
      // 	std::cout<<"COULD NOT FIND ENOUGH FFT PRIME in MatPoly FFTMUL taking normal primes..."<<std::endl;
      //        exit(1);
      // }
	
      std::vector<double> basis(bas.size());
      std::copy(bas.begin(),bas.end(),basis.begin());
      FFPACK::rns_double RNS(basis);
      size_t num_primes = RNS._size;
#ifdef FFT_PROFILER
      //double tMul=0.,tCopy=0;;
      if (FFT_PROF_LEVEL<3){
	std::cout << "*** MatPoly FFT - MUL ***"<<std::endl;
	std::cout << "number of FFT primes :" << num_primes << std::endl;
	std::cout << "max prime            : "<<prime_max<<" ("<<integer(prime_max).bitsize()<<")"<<std::endl;
	std::cout << "bitsize of the output: "<<bound.bitsize()
		  <<"( "<< RNS._M.bitsize()<<" )"<<std::endl;
	std::cout <<" +++++++++++++++++++++++++++++++"<<std::endl;
      }
#endif
      FFT_PROFILING(2,"init of CRT approach");
      // reduce t_a and t_b modulo each FFT primes
      size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
      ADD_MEM(8*(n_ta+n_tb)*num_primes);
      double* t_a_mod= new double[n_ta*num_primes];
      double* t_b_mod= new double[n_tb*num_primes];
      RNS.init(1, n_ta, t_a_mod, n_ta, a.getPointer(), n_ta, maxA);
      RNS.init(1, n_tb, t_b_mod, n_tb, b.getPointer(), n_tb, maxB);
      ADD_MEM(n_ta* (maxA.bitsize()/16 + (maxA.bitsize()%16?1:0)) *8); // needed by RNS init
      DEL_MEM(n_ta* (maxA.bitsize()/16 + (maxA.bitsize()%16?1:0)) *8);
      ADD_MEM(n_tb* (maxB.bitsize()/16 + (maxB.bitsize()%16?1:0)) *8); // needed by RNS init
      DEL_MEM(n_tb* (maxB.bitsize()/16 + (maxB.bitsize()%16?1:0)) *8);

      FFT_PROFILING(2,"reduction mod pi of input matrices");

      std::vector<MatrixP_F*> c_i (num_primes);
      
      for (size_t l=0;l<num_primes;l++)
	{
	  //FFT_PROFILE_START;
	  ModField f(RNS._basis[l]);
	  MatrixP_F a_i (f, m, k, pts);
	  MatrixP_F b_i (f, k, n, pts);
	  //a_i.changeField(f);
	  //b_i.changeField(f);

	  c_i[l] = new MatrixP_F(f, m, n, pts);

	  // copy reduced data
	  for (size_t i=0;i<m*k;i++)
	    for (size_t j=0;j<a.size();j++)
	      a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
	  for (size_t i=0;i<k*n;i++)
	    for (size_t j=0;j<b.size();j++)
	      b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];	
	  //std::cout<<"a"<<l<<":="<<a_i<<";\n";
	  //std::cout<<"b"<<l<<":="<<b_i<<";\n";
	  
	  //FFT_PROFILE_GET(tCopy);
	  //PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain (f);
	  PolynomialMatrixThreePrimesFFTMulDomain<ModField> fftdomain (f);
	  integer bound=integer(RNS._basis[l]-1)*integer(RNS._basis[l]-1)
	    *integer((uint64_t) k)*integer((uint64_t)std::min(a.size(),b.size()));

	  fftdomain.mul_fft(lpts, *c_i[l], a_i, b_i, bound);
	  //std::cout<<"c"<<l<<":="<<*c_i[l]<<";\n";
	  //std::cout<<"p"<<l<<":="<<uint64_t(RNS._basis[l])<<";\n";
	  //FFT_PROFILE_GET(tMul);
	}
      //std::cout<<"MUL FFT RNS: output polmat -> allocating "<<MB(num_primes*c_i[0]->realmeminfo())<<"Mo"<<std::endl;
      //)
      FFT_PROFILING(2,"FFTprime mult+copying");
      DEL_MEM(8*(n_ta+n_tb)*num_primes);
      delete[] t_a_mod;
      delete[] t_b_mod;
      //FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
      //FFT_PROFILE(2,"FFTprime multiplication",tMul);
      
      if (num_primes < 2) {
	FFT_PROFILE_START(2);	
	c.copy(*c_i[0],0,s-1);
      } else {
	FFT_PROFILE_START(2);
	// construct contiguous storage for c_i
	size_t n_tc=m*n*s;
	//std::cout<<"MUL FFT RNS: RNS -> allocating "<<MB(n_tc*num_primes*8)<<"Mo"<<std::endl;
	ADD_MEM(8*n_tc*num_primes);
	double *t_c_mod = new double[n_tc*num_primes];
	//std::cout<<"MUL FFT RNS: output RNS -> allocating "<<MB((n_tc)*num_primes*8)<<"Mo"<<std::endl;
	for (size_t l=0;l<num_primes;l++){
	  for (size_t i=0;i<m*n;i++)
	    for (size_t j=0;j<s;j++)
	      t_c_mod[l*n_tc + (j+i*s)]= c_i[l]->get(i,j);
	  delete c_i[l];
	}
	FFT_PROFILING(2,"linearization of results mod pi");

	// reconstruct the result in C
	RNS.convert(1,n_tc,0,c.getPointer(),n_tc, t_c_mod, n_tc);
	ADD_MEM(n_tc*RNS._ldm*8);
	DEL_MEM(n_tc*RNS._ldm*8);

	//std::cout<<"MUL FFT RNS: "<<MEMINFO<<std::endl;
	//std::cout<<"----------------------------------------------"<<std::endl;
	DEL_MEM(8*n_tc*num_primes);
	delete[] t_c_mod;

      }
      FFT_PROFILING(2,"k prime reconstruction");
      // std::cout<<"CC:="<<c<<std::endl;
      // std::cout<<"<-----------------: "<<std::endl;;
    }

    // WARNING: Polynomial Matrix should stored as matrix of polynomial with integer coefficient 
    template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
    void mul_crtla2(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
		    const integer& maxA, const integer& maxB, const integer& bound) const {

      FFT_PROFILE_START(2);
      linbox_check(a.coldim() == b.rowdim());
      size_t m = a.rowdim();
      size_t k = a.coldim();
      size_t n = b.coldim();
      size_t s= a.size()+b.size()-1;
      c.resize(s);
      size_t lpts=0;
      size_t pts  = 1; while (pts < s) { pts= pts<<1; ++lpts; }

      // compute max prime value for FFLAS      
      uint64_t prime_max= std::min(uint64_t(std::sqrt( (1ULL<<53) / k)+1), uint64_t(Givaro::Modular<double>::maxCardinality()));
      std::vector<integer> bas;
      getFFTPrime(prime_max,lpts,bound,bas,k,s);
      
      std::vector<double> basis(bas.size());
      std::copy(bas.begin(),bas.end(),basis.begin());
      FFPACK::rns_double RNS(basis);
      size_t num_primes = RNS._size;


#ifdef FFT_PROFILER
      //double tMul=0.,tCopy=0;;
      if (FFT_PROF_LEVEL<3){
	std::cout << "number of FFT primes :" << num_primes << std::endl;
	std::cout << "max prime            : "<<prime_max<<" ("<<integer(prime_max).bitsize()<<")"<<std::endl;
	std::cout << "bitsize of the output: "<<bound.bitsize()
		  <<"( "<< RNS._M.bitsize()<<" )"<<std::endl;
	std::cout <<" +++++++++++++++++++++++++++++++"<<std::endl;
      }
#endif


      FFT_PROFILING(2,"init of CRT approach");
      // reduce t_a and t_b modulo each FFT primes
      size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
      std::vector<MatrixP_F*> c_i (num_primes);

 
      // loop for memory saving
      size_t CRT_NBPRIME=4;
      ADD_MEM(8*(n_ta+n_tb)*CRT_NBPRIME);
      double* t_a_mod= new double[n_ta*CRT_NBPRIME];
      double* t_b_mod= new double[n_tb*CRT_NBPRIME];
            
      for(size_t loop=0;loop<num_primes;loop+=CRT_NBPRIME){
	
	// create chunk of RNS
	size_t rns_chunk=std::min(CRT_NBPRIME,num_primes-loop); // nbr of primes in the current smallRNS basis
	std::vector<double> smallBasis(rns_chunk);
	std::copy(basis.begin()+loop,basis.begin()+loop+rns_chunk,smallBasis.begin());
	FFPACK::rns_double smallRNS(smallBasis);
	smallRNS.precompute_cst(RNS._ldm);

	smallRNS.init(1, n_ta, t_a_mod, n_ta, a.getPointer(), n_ta, maxA);
	smallRNS.init(1, n_tb, t_b_mod, n_tb, b.getPointer(), n_tb, maxB);
	FFT_PROFILING(2,"reduction mod pi of input matrices");

	for (size_t l=0;l<rns_chunk;l++)
	  {	    
	    //FFT_PROFILE_START;
	    //std::cout<<"prime: "<<(long)smallRNS._basis[l]<<std::endl;
	    ModField f(smallRNS._basis[l]);
	    MatrixP_F a_i (f, m, k, pts);
	    MatrixP_F b_i (f, k, n, pts);	
	    c_i[loop+l] = new MatrixP_F(f, m, n, pts);
	    // copy reduced data
	    for (size_t i=0;i<m*k;i++)
	      for (size_t j=0;j<a.size();j++)
		a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
	    for (size_t i=0;i<k*n;i++)
	      for (size_t j=0;j<b.size();j++)
		b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];	
	    //FFT_PROFILE_GET(tCopy);
	    //PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain (f);
	    PolynomialMatrixThreePrimesFFTMulDomain<ModField> fftdomain (f);
	    integer bound=integer(smallRNS._basis[l]-1)*integer(smallRNS._basis[l]-1)
	      *integer(uint64_t(k))*integer((uint64_t)std::min(a.size(),b.size()));
	    
	    fftdomain.mul_fft(lpts, *c_i[loop+l], a_i, b_i, bound);	
	    //FFT_PROFILE_GET(tMul);
	  }      
	FFT_PROFILING(2,"FFTprime mult+copying");
	//FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
	//FFT_PROFILE(2,"FFTprime multiplication",tMul);

      } // end of loop for memory saving
      DEL_MEM(8*(n_ta+n_tb)*CRT_NBPRIME);
      delete[] t_a_mod;
      delete[] t_b_mod;
      
      if (num_primes < 2) {
	FFT_PROFILE_START(2);	
	c.copy(*c_i[0],0,s-1);
      } else {
	FFT_PROFILE_START(2);
	// construct contiguous storage for c_i
	size_t n_tc=m*n*s;
	//std::cout<<"MUL FFT RNS: RNS -> allocating "<<MB(n_tc*num_primes*8)<<"Mo"<<std::endl;
	ADD_MEM(8*(n_tc)*num_primes);
	double *t_c_mod = new double[n_tc*num_primes];
	for (size_t l=0;l<num_primes;l++){
	  for (size_t i=0;i<m*n;i++)
	    for (size_t j=0;j<s;j++)
	      t_c_mod[l*n_tc + (j+i*s)]= c_i[l]->get(i,j);
	  delete c_i[l];
	}
	FFT_PROFILING(2,"linearization of results mod pi");

	// reconstruct the result in C
	RNS.convert(1,n_tc,0,c.getPointer(),n_tc, t_c_mod, n_tc);
	ADD_MEM(n_tc*RNS._ldm*8);
	DEL_MEM(n_tc*RNS._ldm*8);
	
	//std::cout<<"MUL FFT RNS: "<<MEMINFO<<std::endl;
	//std::cout<<"----------------------------------------------"<<std::endl;
	DEL_MEM(8*n_tc*num_primes);
	delete[] t_c_mod;
      }
      FFT_PROFILING(2,"k prime reconstruction");
      // std::cout<<"CC:="<<c<<std::endl;
      // std::cout<<"<-----------------: "<<std::endl;;
    }




    // template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
    // void midproduct_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
    // 			  const integer& maxA, const integer& maxB, const integer& bound,
    // 			  bool smallLeft=true, size_t n0=0, size_t n1=0) {
    //   // (convert to MatrixP representation)
    //   MatrixP_I a2(field(),a.rowdim(),a.coldim(),a.size());
    //   MatrixP_I b2(field(),b.rowdim(),b.coldim(),b.size());
    //   a2.copy(a,0,a.size()-1);
    //   b2.copy(b,0,b.size()-1);
    //   MatrixP_I c2(field(),c.rowdim(),c.coldim(),c.size());
    //   midproduct_crtla(c2,a2,b2,maxA,maxB,bound,smallLeft,n0,n1);
    //   c.copy(c2,0,c2.size()-1);
    // }

    // WARNING: Polynomial Matrix should stored as matrix of polynomial with integer coefficient 
    template< typename PMatrix1,typename PMatrix2, typename PMatrix3>
    void midproduct_crtla(PMatrix1 &c, const PMatrix2 &a, const PMatrix3 &b,
			  const integer& maxA, const integer& maxB, const integer& bound,
			  bool smallLeft=true, size_t n0=0, size_t n1=0) const {
      // void midproduct_crtla(MatrixP_I &c, const MatrixP_I &a, const MatrixP_I &b,
      // 			  const integer& maxA, const integer& maxB, const integer& bound,
      // 			  bool smallLeft=true, size_t n0=0, size_t n1=0) {
      FFT_PROFILE_START(2);
      linbox_check(a.coldim() == b.rowdim());
      size_t m = a.rowdim();
      size_t k = a.coldim();
      size_t n = b.coldim();
      size_t hdeg = (n0==0?c.size():n0);
      size_t deg  = (n1==0?2*hdeg-1:n1);
      linbox_check(c.size()>=deg-hdeg);
      
      if (smallLeft){
      	linbox_check(b.size()<hdeg+deg);
      }
      else
      	linbox_check(a.size()<hdeg+deg);

      //linbox_check(2*c.size()-1 == b.size());
      //size_t deg= b.size()+1;
      //size_t hdeg= deg/2;
      size_t lpts=0;
      size_t pts  = 1; while (pts < deg) { pts= pts<<1; ++lpts; }


      // compute bit size of feasible prime for FFLAS
      // size_t _k=k,lk=0;
      //while ( _k ) {_k>>=1; ++lk;}
      //size_t prime_bitsize= (53-lk)>>1;

      // compute max prime value for FFLAS
      uint64_t prime_max= std::sqrt( (1ULL<<53) / k)+1;
      std::vector<integer> bas;
      getFFTPrime(prime_max,lpts,bound,bas,k,deg);

      std::vector<double> basis(bas.size());
      std::copy(bas.begin(),bas.end(),basis.begin());
      FFPACK::rns_double RNS(basis);
      size_t num_primes = RNS._size;
#ifdef FFT_PROFILER
      double tMul=0.,tCopy=0;;
      if (FFT_PROF_LEVEL<3){
	std::cout << "*** MatPoly FFT - MIDP ***"<<std::endl;
 	std::cout << "number of FFT primes :" << num_primes << std::endl;
	std::cout << "max prime            : "<<prime_max<<" ("<<integer(prime_max).bitsize()<<")"<<std::endl;
	std::cout << "bitsize of the output: "<<bound.bitsize()
		  <<"( "<< RNS._M.bitsize()<<" )"<<std::endl;
	std::cout <<" +++++++++++++++++++++++++++++++"<<std::endl;
      }
#endif
      FFT_PROFILING(2,"init of CRT approach");
      // reduce t_a and t_b modulo each FFT primes
      size_t n_ta=m*k*a.size(), n_tb=k*n*b.size();
      ADD_MEM(8*(n_ta+n_tb)*num_primes);
      double* t_a_mod= new double[n_ta*num_primes];
      double* t_b_mod= new double[n_tb*num_primes];
      //std::cout<<"MIDP FFT RNS: input RNS -> allocating "<<MB((n_ta+n_tb)*num_primes*8)<<"Mo"<<std::endl;

      RNS.init(1, n_ta, t_a_mod, n_ta, a.getPointer(), n_ta, maxA);
      RNS.init(1, n_tb, t_b_mod, n_tb, b.getPointer(), n_tb, maxB);
      ADD_MEM(n_ta* (maxA.bitsize()/16 + (maxA.bitsize()%16?1:0)) *8); // needed by RNS init
      DEL_MEM(n_ta* (maxA.bitsize()/16 + (maxA.bitsize()%16?1:0)) *8);
      ADD_MEM(n_tb* (maxB.bitsize()/16 + (maxB.bitsize()%16?1:0)) *8); // needed by RNS init
      DEL_MEM(n_tb* (maxB.bitsize()/16 + (maxB.bitsize()%16?1:0)) *8);

      FFT_PROFILING(2,"reduction mod pi of input matrices");

      std::vector<MatrixP_F*> c_i (num_primes);

      for (size_t l=0;l<num_primes;l++){
	FFT_PROFILE_START(2);
	ModField f(RNS._basis[l]);
	MatrixP_F a_i (f, m, k, pts);
	MatrixP_F b_i (f, k, n, pts);
	c_i[l] = new MatrixP_F(f, m, n, pts);
	// copy reduced data and reversed when necessary according to midproduct algo
	//std::cout<<"hdeg-size: "<<hdeg<<" <-> "<<a.size()<<std::endl;
	for (size_t i=0;i<m*k;i++)
	  for (size_t j=0;j<a.size();j++)
	    if (smallLeft)
	      a_i.ref(i,hdeg-1-j)=t_a_mod[l*n_ta+j+i*a.size()];
	    else
	      a_i.ref(i,j)=t_a_mod[l*n_ta+j+i*a.size()];
	for (size_t i=0;i<k*n;i++)
	  for (size_t j=0;j<b.size();j++)
	    if (smallLeft)
	      b_i.ref(i,j)=t_b_mod[l*n_tb+j+i*b.size()];
	    else
	      b_i.ref(i,hdeg-1-j)=t_b_mod[l*n_tb+j+i*b.size()];
	FFT_PROFILE_GET(2,tCopy);
	//PolynomialMatrixFFTPrimeMulDomain<ModField> fftdomain (f);
	PolynomialMatrixThreePrimesFFTMulDomain<ModField> fftdomain (f);       
	integer bound2=integer(RNS._basis[l]-1)*integer(RNS._basis[l]-1)
	  *integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));
	fftdomain.midproduct_fft(lpts, *(c_i[l]), a_i, b_i, bound2, smallLeft);
				
	FFT_PROFILE_GET(2,tMul);
      }

      DEL_MEM(8*(n_ta+n_tb)*num_primes);
      delete[] t_a_mod;
      delete[] t_b_mod;

      FFT_PROFILE(2,"copying linear reduced matrix",tCopy);
      FFT_PROFILE(2,"FFTprime multiplication",tMul);

      if (num_primes < 2) {
	FFT_PROFILE_START(2);
	c.copy(*(c_i[0]),0,c.size()-1);
      } else {
	FFT_PROFILE_START(2);
	// construct contiguous storage for c_i
	double *t_c_mod;
	size_t n_tc=m*n*c.size();
	ADD_MEM(8*(n_tc)*num_primes);
	t_c_mod = new double[n_tc*num_primes];
	//std::cout<<"MIDP FFT RNS: output RNS -> allocating "<<MB((n_tc)*num_primes*8)<<"Mo"<<std::endl;
	for (size_t l=0;l<num_primes;l++){
	  for (size_t i=0;i<m*n;i++)
	    for (size_t j=0;j<c.size();j++)
	      t_c_mod[l*n_tc + (j+i*c.size())]= c_i[l]->get(i,j);
	  delete c_i[l];
	}
	FFT_PROFILING(2,"linearization of results mod pi");

	// reconstruct the result in C
	RNS.convert(1,n_tc,0,c.getPointer(),n_tc, t_c_mod, n_tc);
	ADD_MEM(n_tc*RNS._ldm*8); // needed by RNS
	DEL_MEM(n_tc*RNS._ldm*8);

	DEL_MEM(8*n_tc*num_primes);
	delete[] t_c_mod;
	
	FFT_PROFILING(2,"k prime reconstruction");
      }
    }
  };


  /***************************************************************************
   **** Polynomial Matrix Multiplication over Fp[x], with p multiprecision ***
   ***************************************************************************/
  template <>
  class PolynomialMatrixFFTMulDomain<Givaro::Modular<integer> > {
  public:
    typedef Givaro::Modular<integer>              Field;
    typedef typename Field::Element     Element;
    typedef Givaro::ZRing<integer>  IntField;
      // Polynomial matrix stored as a polynomial of matrix
      typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP_F;
      // Polynomial matrix stored as a polynomial of matrix
      typedef PolynomialMatrix<IntField, PMType::polfirst> MatrixP_I;

  private:
    const Field            *_field;  // Read only
    integer                     _p;

  public:
    inline const Field & field() const { return *_field; }

    PolynomialMatrixFFTMulDomain(const Field &F) : _field(&F) {
      field().cardinality(_p);
    }

    template<typename Matrix1, typename Matrix2, typename Matrix3>
    void mul (Matrix1 &c, const Matrix2 &a, const Matrix3 &b, size_t max_rowdeg=0) const {
      FFT_PROFILE_START(2);
      MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
      MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
      MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
      a2.copy(a,0,a.size()-1);
      b2.copy(b,0,b.size()-1);
      FFT_PROFILING(2,"converting rep of input");
      mul(c2,a2,b2, max_rowdeg);
      FFT_PROFILE_START(2);
      c.copy(c2,0,c.size()-1);
      FFT_PROFILING(2,"converting rep of output");

    }
    
    // Matrix with polynomials  
    void mul (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b, size_t max_rowdeg=0) const {

      FFT_PROFILE_START(2);
      IntField Z;      
      PolynomialMatrixFFTMulDomain<IntField> Zmul(Z,_p);
      integer bound=2*_p*_p*integer((uint64_t)a.coldim())*integer((uint64_t)std::min(a.size(),b.size()));
#ifdef TRY1
      Zmul.mul_crtla2(c,a,b,_p,_p,bound); 
#else
      Zmul.mul_crtla(c,a,b,_p,_p,bound, max_rowdeg);
#endif
      
      
      // reduce the result mod p
      FFT_PROFILE_START(2);
      for (size_t i=0;i<c.rowdim()*c.coldim();i++)
	for (size_t j=0;j<c.size();j++)
	  c.ref(i,j)%=_p;
      FFT_PROFILING(2,"reduction mod p of output");
    }

    template<typename Matrix1, typename Matrix2, typename Matrix3>
    void midproduct (Matrix1 &c, const Matrix2 &a, const Matrix3 &b,
		     bool smallLeft=true, size_t n0=0, size_t n1=0) const {

      MatrixP_F a2(field(),a.rowdim(),a.coldim(),a.size());
      MatrixP_F b2(field(),b.rowdim(),b.coldim(),b.size());
      a2.copy(a,0,a.size()-1);
      b2.copy(b,0,b.size()-1);
      MatrixP_F c2(field(),c.rowdim(),c.coldim(),c.size());
      midproduct(c2,a2,b2,smallLeft,n0,n1);
      c.copy(c2,0,c.size()-1);
    }

    void midproduct (MatrixP_F &c, const MatrixP_F &a, const MatrixP_F &b,
		     bool smallLeft=true, size_t n0=0, size_t n1=0) const {
      IntField Z;
      PolynomialMatrixFFTMulDomain<IntField> Zmul(Z,_p);
      //const MatrixP_I* a2 = reinterpret_cast<const MatrixP_I*>(&a);
      //const MatrixP_I* b2 = reinterpret_cast<const MatrixP_I*>(&b);
      //MatrixP_I* c2       = reinterpret_cast<MatrixP_I*>(&c);
      //Zmul.midproduct(*c2,*a2,*b2,smallLeft,n0,n1);
      Zmul.midproduct(c,a,b,smallLeft,n0,n1);
      // reduce the result mod p
      FFT_PROFILE_START(2);
      for (size_t i=0;i<c.rowdim()*c.coldim();i++)
	for (size_t j=0;j<c.size();j++)
	  c.ref(i,j)%=_p;
      FFT_PROFILING(2,"reduction mod p of output");
    }
  };

} // end of namespace LinBox

#endif // __LINBOX_matpoly_mult_ftt_multiprecision_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
