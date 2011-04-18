/* linbox/algorithms/lifting-container-base.h
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi pascal.giorgi@ens-lyon.fr
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

#ifndef __LINBOX_lifting_container_H
#define __LINBOX_lifting_container_H

//#define DEBUG_LC
//#define LC_CHECK_DIVISION
#undef _L

#include <vector>

#include <linbox/linbox-config.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/apply.h>
#include <linbox/algorithms/blackbox-container.h>
#include <linbox/algorithms/massey-domain.h>
#include <linbox/algorithms/blackbox-block-container.h>
#include <linbox/algorithms/block-massey-domain.h>
#include <linbox/algorithms/gauss.h>
#include <linbox/blackbox/blas-blackbox.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/compose.h>
#include <linbox/blackbox/block-hankel-inverse.h>
#include <linbox/algorithms/blas-domain.h>
#include <linbox/field/hom.h>
#include <linbox/matrix/transpose-matrix.h>
#include <linbox/blackbox/transpose.h>
//#include <linbox/algorithms/vector-hom.h>
#undef _U

namespace LinBox 
{

	/**
	 * BoundBlackbox: Sets
	 *      H_col_sqr <- H_col(A)^2,   short_col_sqr <- short_col(A)^2
	 * where H_col(A) is prod_j sqrt(sum_i a_ij^2)     ('Hadamard column bound')
	 *   short_col(A) is min_j  sqrt(sum_i a_ij^2)     ('shortest column')
	 *
	 * note: H_col is not actually a norm! but it is what we need for lifting bound computation
	 */
	template <class Ring, class ItMatrix>
	void SpecialBound(const Ring& R, typename Ring::Element& H_col_sqr, 
			  typename Ring::Element& short_col_sqr, const ItMatrix& A) {
					
		typedef typename Ring::Element Integer;
		//Integer sqsum;
		//size_t m, n, col=0;
		//n=A.coldim();
		//m=A.rowdim();
		R.init(H_col_sqr, 1);
		Integer zero;
		R.init(zero,0UL);
		
		typename ItMatrix::ConstRowIterator row= A.rowBegin();
		std::vector<Integer> tmp(A.coldim(), zero);
		for (; row != A.rowEnd(); row++){
			typename ItMatrix::ConstRow::const_iterator elm= row->begin();
			for (size_t i=0; elm != row->end(); elm++, ++i)
				R.axpyin(tmp[i], *elm, *elm);			
		}

		for (size_t i=0;i<A.coldim();++i)
			R.mulin(H_col_sqr,tmp[i]);
		short_col_sqr= *(std::min_element(tmp.begin(),tmp.end()));

		/* at this point RowIterator is better than ColIterator
		typename ItMatrix::ConstColIterator colIter;
		colIter = A.colBegin();

		for (; colIter != A.colEnd(); colIter++, col++) {
			typename ItMatrix::ConstCol::const_iterator elm;
			R.init(sqsum, 0);
			for (elm = colIter->begin(); elm != colIter->end(); elm++)
				R.axpyin(sqsum, *elm, *elm);
			R.mulin(H_col_sqr, sqsum);
			if (col == 0 || sqsum < short_col_sqr)
				short_col_sqr = sqsum;
		}
		*/
		
	}

	// in solveNonsingular, we may work with something that inherits from DenseMatrixBase
	template <class Ring>
	void BoundBlackbox(const Ring& R, typename Ring::Element& H_col_sqr, 
			   typename Ring::Element& short_col_sqr, const DenseMatrixBase<typename Ring::Element>& A) {
		SpecialBound(R, H_col_sqr, short_col_sqr, A);
	}
	template <class Ring>

	void BoundBlackbox(const Ring& R, typename Ring::Element& H_col_sqr, 
			   typename Ring::Element& short_col_sqr, const BlasBlackbox<Ring>& A) {
		SpecialBound(R, H_col_sqr, short_col_sqr, A);
	}

	// in other solvers we generally use BlasBlackbox which inherits from DenseSubmatrix
	template <class Ring>
	void BoundBlackbox(const Ring& R, typename Ring::Element& H_col_sqr, 
			   typename Ring::Element& short_col_sqr, const DenseSubmatrix<typename Ring::Element>& A) {
		SpecialBound(R, H_col_sqr, short_col_sqr, A);
	}


	template < class Ring, class Blackbox>
	void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const Blackbox& A) {

		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
	
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.apply(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter){
				sqsum += (*iter)*(*iter);
			}
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}

	template < class Ring, class Matrix1, class Matrix2>
	void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const Compose<Matrix1,Matrix2> & A) {
		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.apply(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}
	
	template < class Ring, class Matrix>
	void BoundBlackbox (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const Transpose<Matrix> & A) {
		typedef typename Ring::Element Integer;
		Integer one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		R.init(H_col_sqr, 1);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.applyTranspose(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			R.mulin(H_col_sqr, sqsum);
			if (i==0 || sqsum < short_col_sqr) 
				short_col_sqr = sqsum;
			e[i]=zero;
		}
	}


	/* 
	   This should work with blackboxes. However it is much slower if column iterators are available.
	   Furthermore the compiler always binds to this instead of the above faster version; so some 
	   trickier kind of specialization may have to be done when BoundBlackBox is to be used with true blackboxes.
	   (Or is the plural Blackboxen?)

	   template < class Ring, class IMatrix>
	   void BoundBlackbox 
	   (const Ring& R, typename Ring::Element& H_col_sqr, typename Ring::Element& short_col_sqr, const IMatrix& A) {
	   typedef typename Ring::Element Integer;
	   Integer one,zero,sqsum;
	   size_t m,n;
	   n=A.coldim();
	   m=A.rowdim();
	   R.init(one,1);
	   R.init(zero,0);
	   R.init(H_col_sqr, 1);
	   typename std::vector<Integer>::const_iterator iter;
	   std::vector<Integer> e(n,zero),tmp(m);
	   for (size_t i=0;i<n;i++){
	   e[i]=one;
	   A.apply(tmp,e);
	   sqsum=zero;
	   for (iter=tmp.begin();iter!=tmp.end();++iter)
	   sqsum += (*iter)*(*iter);
	   R.mulin(H_col_sqr, sqsum);
	   if (i==0 || sqsum < short_col_sqr) 
	   short_col_sqr = sqsum;
	   e[i]=zero;
	   }
	   }
	*/

	/**
	 * ApplyBound: computes
	 *         bound_A <- max_i(max(sum_{j|a_ij > 0} a_ij, sum_{j|a_ij < 0} |a_ij|))
	 * this is useful because for all u, v >= 0:
	 *     [b has all entries in -u..v] 
         *       => [each entry of A.b is at most (u+v)*bound_A in absolute value]
	 */
	template <class Ring, class ItMatrix> //iterable matrix
	void ApplyBound(const Ring& R, typename Ring::Element& bound_A, const ItMatrix& A) {
		typedef typename Ring::Element Integer;
		Integer possum, negsum, zero;
		R.init(bound_A, 0);
		R.init(zero, 0);
		
		typename ItMatrix::ConstRowIterator rowIter;
		rowIter = A.rowBegin();

		for (; rowIter != A.rowEnd(); rowIter++) {
			typename ItMatrix::ConstRow::const_iterator elm;
			R.init(possum, 0);
			R.init(negsum, 0);
			for (elm = rowIter->begin(); elm != rowIter->end(); elm++) 
				if (*elm > zero) 
					R.addin(possum, *elm);
				else
					R.subin(negsum, *elm);
			
			if (possum > bound_A)
				R.assign(bound_A, possum);
			if (negsum > bound_A)
				R.assign(bound_A, negsum);
		}
	}

	template< class _Ring>
	class LiftingContainer {
		
	public:
		typedef _Ring Ring;
		typedef typename Ring::Element Integer;
				
		// return the length of container
		virtual const size_t length() const =0;

		// return the size of the solution
		virtual const size_t size() const = 0;

		// return the ring
		virtual const Ring& ring() const = 0;
				
		// return the prime
		virtual const Integer& prime () const = 0; 		

		virtual ~LiftingContainer () {}
		
	};

	template< class _Ring, class _IMatrix>
	class LiftingContainerBase : public LiftingContainer< _Ring> {

	public:		
		typedef _IMatrix                  IMatrix;
		typedef _Ring                        Ring;
		typedef typename _Ring::Element   Integer;		
		typedef std::vector<Integer>      IVector;
#ifdef RSTIMING
		mutable Timer ttSetup, tRingApply, tRingOther, ttRingOther, ttRingApply;
#endif

	protected:
						
		const IMatrix&            _A;
		Ring                      _R;
		Integer                   _p;
		IVector                   _b;
		VectorDomain<Ring>      _VDR;
		size_t               _length;
		Integer            _numbound;
		Integer            _denbound;
		MatrixApplyDomain<Ring,IMatrix>    _MAD;
		//BlasApply<Ring>          _BA;

	
	

		void convertPrime(Integer& e, const integer& p){
			_R.init(e,p);
		}
		

		void convertPrime(Integer& e, const std::vector<integer>& p){
			integer tmp=1;
			for (size_t i=0;i<p.size();++i)
				tmp*=integer(p[i]);			
			_R.init(e,tmp);
		}
	

	public:

		template <class Prime_Type, class Vector1>
		LiftingContainerBase (const Ring& R, const IMatrix& A, const Vector1& b, const Prime_Type& p):
			_A(A), _R(R), _VDR(R), _MAD(R,A) {
		
#ifdef RSTIMING
			ttSetup.start();
#endif
			linbox_check(A.rowdim() == b.size());
			int n,m;
			n=A.rowdim();
			m=A.coldim();

			//assert(m == n); //logic may not work otherwise
			linbox_check( m == n );
			// initialise the prime as an Integer
			//this->_R.init(_p,p);
			this->convertPrime(_p, p);
			//std::cout<<"padic base= "<<_p<<std::endl;
			
			
			// initialize res = b
			_b.resize(b.size());
			typename Vector1::const_iterator         b_iter    = b.begin();
			typename std::vector<Integer>::iterator  res_iter  = _b.begin() ;
			for (; b_iter != b.end(); ++res_iter, ++b_iter) 
				this->_R.init(*res_iter, *b_iter);
						
			Integer had_sq, short_sq;
			BoundBlackbox(this->_R, had_sq, short_sq, A);
			
			typename std::vector<Integer>::const_iterator iterb = _b.begin();
			Integer normb_sq;
			this->_R.init(normb_sq, 0);
			for (;iterb!=_b.end();++iterb)
				normb_sq += (*iterb)*(*iterb);

			LinBox::integer had_sqi, short_sqi, normb_sqi, N, D, L, prime;
			this->_R.convert(had_sqi, had_sq);
			this->_R.convert(short_sqi, short_sq);
			this->_R.convert(normb_sqi, normb_sq);
			this->_R.convert(prime,_p);
			D = sqrt(had_sqi) + 1;
			N = sqrt(had_sqi * normb_sqi / short_sqi) + 1;
			L = N * D * 2; 
			_length = logp(L,prime) + 1;   // round up instead of down
#ifdef DEBUG_LC                                   
			std::cout<<" norms computed, p = "<<_p<<"\n";
			std::cout<<" N = "<<N<<", D = "<<D<<", length = "<<_length<<"\n";
			std::cout<<"A:=\n";
			//_A.write(std::cout);
			std::cout<<"b:=\n";
			for (size_t i=0;i<_b.size();++i) std::cout<<_b[i]<<" , ";
			std::cout<<std::endl;
#endif
			this->_R.init(_numbound,N);
			this->_R.init(_denbound,D);
	
			_MAD.setup( prime );		
			
#ifdef DEBUG_LC		
			std::cout<<"lifting container initialized\n";			
#endif
#ifdef RSTIMING
			ttSetup.stop();
			ttRingOther.clear();
			ttRingApply.clear();
#endif
		}
		
		virtual IVector& nextdigit (IVector& , const IVector&) const = 0;

		class const_iterator {
		private:
			std::vector<Integer>          _res;
			const LiftingContainerBase    &_lc;
			size_t                   _position;
		public:
			const_iterator(const LiftingContainerBase& lc,size_t end=0)
				: _res(lc._b), _lc(lc), _position(end) {}					
			
			/**
			 * @returns False if the next digit cannot be computed
			 * (probably indicates modulus is bad)
			 */
			bool next (IVector& digit)  {
				
				linbox_check (digit.size() == _lc._A.rowdim());	  
				// compute next p-adic digit
				_lc.nextdigit(digit,_res);
#ifdef RSTIMING
				_lc.tRingApply.start();
#endif

#ifdef DEBUG_LC
				std::cout<<"\n residu "<<_position<<": ";
				for (size_t i=0;i<digit.size();++i)
					std::cout<<_res[i]<<",";
				std::cout<<"\n digit "<<_position<<": ";			
				for (size_t i=0;i<digit.size();++i)
					std::cout<<digit[i]<<",";
				std::cout<<"\n";
#endif 	
				/*  prepare for updating residu */

				// compute v2 = _A * digit				
				IVector v2 (_lc._A.coldim());
				_lc._MAD.applyV(v2,digit, _res);

#ifdef DEBUG_LC
				
				//_A.write(std::cout<<"\n _A :\n"); 
				std::cout<<"\n A * digit "<<_position<<": ";
				for (size_t i=0;i<v2.size();++i)
					std::cout<<v2[i]<<",";

#endif
#ifdef RSTIMING
				_lc.tRingApply.stop();
				_lc.ttRingApply += _lc.tRingApply;
				_lc.tRingOther.start();
#endif

				// update _res -= v2 
				_lc._VDR.subin (_res, v2);
				typename std::vector<Integer>::iterator p0;
				// update _res = _res / p
				int index=0;
				for ( p0 = _res.begin(); p0 != _res.end(); ++ p0, ++index){ 
#ifdef LC_CHECK_DIVISION
 					if (! _lc._R.isDivisor(*p0,_lc._p)) {
						std::cout<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<std::endl;
 						std::cout<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<std::endl;
 						return false;
 					}
#endif
					_lc._R.divin(*p0, _lc._p);
				}
				
				// increase position of the iterator
				_position++;
#ifdef RSTIMING				
				_lc.tRingOther.stop();
				_lc.ttRingOther += _lc.tRingOther;
#endif
				return true;
			}

			bool operator != (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					;//assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position != iterator._position;					
			}
	
			bool operator == (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					;//assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position == iterator._position;					
			}

		};

		/** Bit manipulation function for possible use in optimization
		 // efficiently pulls out continuous blocks of bits, from lsb to msb inclusive
		 // least significant bits start at index 0, so msb >= lsb
		 // if any bits with index >= 8*numBytes are asked for they will be zeroes
		 static long long bytesToBits(unsigned char * byteArray, size_t numBytes, size_t lsb, size_t msb) { 
		 linbox_check(msb >= lsb);
		 size_t lsbi = lsb >> 3;
		 size_t msbi = msb >> 3;
		 if (msbi == lsbi) 
		 if (msbi >= numBytes) 
		 return 0;
		 else
		 return (byteArray[lsbi] >> (lsb & 7)) & ((1 << (msb - lsb + 1)) - 1);

		 long long result = (msbi < numBytes) ? (byteArray[msbi] & ((1 << ((msb & 7)+1)) - 1)) : 0;
		 for (size_t i=msbi-1; i>lsbi; i--) {
		 result <<= 8;
		 result |= (i < numBytes) ? byteArray[i] : 0;
		 }
		 result <<= 8 - (lsb & 7);
		 result |= (lsbi < numBytes) ? (byteArray[lsbi] >> (lsb & 7)) : 0;

		 return result;
		 }
		*/

		const_iterator begin() const {
			return const_iterator(*this);
		}		

		const_iterator end() const {
			return const_iterator (*this,_length);
		}
	
		virtual const size_t length() const
		{ return _length;}
		
		// return the size of the solution
		virtual const size_t size() const 
		{ return _A.coldim();}

		// return the ring
		virtual const Ring& ring() const 
		{ return this->_R; }
		
		// return the prime
		virtual const Integer& prime () const 
		{ return _p; }
		
		// return the bound for the numerator
		const Integer numbound() const
		{return _numbound;}
		
		// return the bound for the denominator
		const Integer denbound() const
		{return _denbound;}
		
		// return the matrix
		const IMatrix& getMatrix() const
		{return _A; }

		// return the right hand side
		const IVector& getVector() const
		{return _b;}
		
	};

	template <class _Ring, class _Field, class _IMatrix, class _FMatrix>
	class DixonLiftingContainer : public LiftingContainerBase< _Ring, _IMatrix> {

	public:
		typedef _Field                               Field;
		typedef _Ring                                 Ring;
		typedef _IMatrix                           IMatrix;		
		typedef _FMatrix                           FMatrix;	
		typedef typename Field::Element            Element;
		typedef typename IMatrix::Element          Integer;
		typedef std::vector<Integer>               IVector;
		typedef std::vector<Element>               FVector;

	protected:
		
		const FMatrix&                  _Ap;	
		Field                            _F;
		const VectorDomain<Field>      _VDF;
		mutable FVector              _res_p;
		mutable FVector            _digit_p;		
		BlasApply<Field>                _BA;
		
	public:
#ifdef RSTIMING
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert;
#endif

		template <class Prime_Type, class VectorIn>
		DixonLiftingContainer (const Ring&       R, 
				       const Field&      F, 
				       const IMatrix&    A, 
				       const FMatrix&   Ap, 
				       const VectorIn&   b, 
				       const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _F(F), _VDF(F), 
			  _res_p(b.size()), _digit_p(A.coldim()), _BA(F) 
		{

			for (size_t i=0; i< _res_p.size(); ++i)
				_F.init(_res_p[i]);
			for (size_t i=0; i< _digit_p.size(); ++i)
				_F.init(_digit_p[i]);
			
			//
#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
#endif	
#ifdef DEBUG_LC
			std::cout<<"Primes: ";
			_F.write(std::cout);
			std::cout<<"\n Matrix: \n";
			A.write(std::cout);
			std::cout<<"\n Matrix mod p: \n";
			Ap.write(std::cout);
			std::cout<<"\n Matrix LCBASE: \n";
			LiftingContainerBase<Ring,IMatrix>::_A.write(std::cout);
#endif

		}
		
		
		virtual ~DixonLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit, const IVector& residu) const {		
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			LinBox::integer tmp;
			
			Hom<Ring, Field> hom(this->_R, _F);
			// res_p =  residu mod p
			//VectorHom::map (_res_p, residu, _F, this->_R);			
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					//_F. init (*iter_p, this->_R.convert(tmp,*iter));
					hom.image(*iter_p, *iter);//std::cout<<*iter_p<<"= "<< *iter<<" mod "<<this->_p<<"\n";}
			}			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
			tGetDigit.start();
#endif			
			
			// compute the solution by applying the inverse of A mod p
			//_BA.applyV(_digit_p,_Ap,_res_p);
			_Ap.apply(_digit_p, _res_p);
#ifdef RSTIMING
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			//VectorHom::map(digit, _digit_p, this->_R, _F);			
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					//this->_R.init(*iter, _F.convert(tmp,*iter_p));
					hom.preimage(*iter, *iter_p);
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif		
			return digit;
		}
		
	}; // end of class DixonLiftingContainerBase



	template <class _Ring, class _Field, class _IMatrix, class _FMatrix, class _FPolynomial>
	class WiedemannLiftingContainer : public LiftingContainerBase<_Ring, _IMatrix> {

	public:
		typedef _Field                                     Field;
		typedef _Ring                                       Ring;
		typedef _IMatrix                                 IMatrix;		
		typedef _FMatrix                                 FMatrix;
		typedef typename Field::Element                  Element;
		typedef typename Ring::Element                   Integer;
		typedef std::vector<Integer>                     IVector;
		typedef std::vector<Element>                     FVector;
		typedef _FPolynomial                         FPolynomial;
		typedef typename FPolynomial::iterator     FPolyIterator;

	protected:

		const FMatrix                  &_Ap;
		mutable FPolynomial        _MinPoly;
		Field                            _F;
		const VectorDomain<Field>      _VDF;
		mutable FVector              _res_p;
		mutable FVector            _digit_p;
		typename Field::RandIter      _rand;		
#ifdef RSTIMING
	public:
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert;
#endif			
	public:

		template <class Prime_Type, class VectorIn>
		WiedemannLiftingContainer (const Ring& R,
					   const Field& F, 
					   const IMatrix& A, 
					   const FMatrix& Ap, 
					   const FPolynomial& MinPoly,
					   const VectorIn& b, 
					   const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _MinPoly(MinPoly), _F(F), _VDF(F), _res_p(b.size()), _digit_p(A.coldim()), _rand(F) {

			// Normalize the minimal polynomial as f(x)=1- a1/a0 x - a2/a0 x^2 - ...
			FPolyIterator iter=_MinPoly.begin();
			while(++iter != _MinPoly.end ()){
				_F.divin (*iter, _MinPoly.front ());
				_F.negin (*iter);
			}
#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
#endif
		}

		virtual ~WiedemannLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit,const IVector& residu) const {
			
			LinBox::integer tmp;
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, this->_R.convert(tmp,*iter));
			}
#ifdef RSTIMING			
			tGetDigitConvert.stop();
			ttGetDigitConvert+=tGetDigitConvert;
			tGetDigit.start();
#endif		
			// compute the solution of system by Minimal polynomial application
			_VDF.mul (_digit_p, _res_p, _MinPoly.back ());      
			FVector z(_Ap.rowdim ());      
			
			
			for (size_t i = _MinPoly.size () - 1; --i > 0;) {												
				_Ap.apply (z, _digit_p);
				_VDF.axpy (_digit_p, _MinPoly[i], _res_p, z);
			}      

			      
			// check results
			FVector error(_Ap.coldim());
			_Ap.apply(error,_digit_p);

			bool nosolution = false;
			int nosolution_threshold=5;
			int nst=0;
			size_t minpoly_degree;
			// until the digit is incorrect update the minpoly and recompute the digit
			while (!_VDF.areEqual(error,_res_p) && !nosolution ){
				minpoly_degree=_MinPoly.size();
				FPolynomial Poly;
				unsigned long deg;
				unsigned long size= (_Ap.rowdim() - _MinPoly.size())<<1 ;
				BlackboxContainer<Field, FMatrix > Sequence(&_Ap,_F,error,size);
				MasseyDomain<Field,BlackboxContainer<Field, FMatrix > > MD(&Sequence);
				MD.minpoly(Poly,deg);
				if (_F.isZero(Poly.front())) {
					// here we should stop the execution but not yet implemented
					std::cout<<" the prime was not good \n, result will be wrong";
					break;
				}
				
				// denormalize the minimal polynomial
				FPolyIterator iter=_MinPoly.begin();
				while (++iter != _MinPoly.end()) {
					_F.mulin (*iter, _MinPoly.front());
					_F.negin (*iter);
				}	
		       	
				// update the minimal polynomial 
				Element zero;
				_F.init(zero,0);
				FPolynomial newMinPoly(_MinPoly.size()+Poly.size()-1,zero);
				for (size_t i=0; i < _MinPoly.size(); i++)
					for (size_t j=0 ; j < Poly.size(); j++)
						_F.axpyin(newMinPoly[i+j],_MinPoly[i],Poly[j]);	
				_MinPoly.clear();
				Poly.clear();
				_MinPoly=newMinPoly;

				// normalize the new minimal polynomial
				iter=_MinPoly.begin ();
				while (++iter != _MinPoly.end ()) {
					_F.divin (*iter, _MinPoly.front ());
					_F.negin (*iter);
				}
				
				_VDF.mul (_digit_p, _res_p, _MinPoly.back ());      
				FVector z(_Ap.rowdim ());      
				for (size_t i = _MinPoly.size () - 1; --i > 0;) {		
					_Ap.apply (z, _digit_p);
					_VDF.axpy (_digit_p, _MinPoly[i], _res_p, z);
				}    

				_Ap.apply(error,_digit_p);	 
				if (_MinPoly.size() > minpoly_degree){
					minpoly_degree = _MinPoly.size();
					nst=0;std::cout<<"updating minpoly\n";
				}
				else {
					if (nst < nosolution_threshold) nst++;
					else{
						nosolution=true;
						throw PreconditionFailed (__FUNCTION__, __LINE__, "system is inconsistent or the choosen prime leads to inconsistent resolution");
					}
				}
			}
#ifdef RSTIMING
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					this->_R.init(*iter, _F.convert(tmp,*iter_p));
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif
			return digit;
		}


	}; // end of class WiedemannLiftingContainerBase


	// Block Wiedemann LiftingContianer
	template <class _Ring, class _Field, class _IMatrix, class _FMatrix>
	class BlockWiedemannLiftingContainer : public LiftingContainerBase<_Ring, _IMatrix> {

	public:
		typedef _Field                                	            Field;
		typedef _Ring                                 	             Ring;
		typedef _IMatrix                              	          IMatrix;       
		typedef _FMatrix                              	          FMatrix;
		typedef typename Field::Element               	          Element;
		typedef typename Ring::Element                            Integer;
		typedef std::vector<Integer>                              IVector;
		typedef std::vector<Element>                              FVector;
		typedef BlasMatrix<Element>                           Coefficient;
		typedef BlasMatrix<Element>                                 Block;
		typedef std::vector<Coefficient>                 FBlockPolynomial;
		typedef BlackboxBlockContainerRecord<Field, FMatrix>     Sequence;


	protected:

		const FMatrix                       &_Ap;
		Field                                 _F;
		const VectorDomain<Field>           _VDF;
		mutable FVector                   _res_p;
		mutable FVector                 _digit_p;
		typename Field::RandIter           _rand;
		size_t                              _row;
		size_t                              _col;
		size_t                                _m;
		size_t                                _n;
		Block                                 _U;
		BlasMatrixDomain<Field>             _BMD;
		Sequence                           *_Seq;
		BlockMasseyDomain<Field,Sequence>  *_Dom;
#ifdef RSTIMING
	public:
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert, tMinPoly, ttMinPoly;
#endif			
	public:

		template <class Prime_Type, class VectorIn>
		BlockWiedemannLiftingContainer (const Ring                         &R,
						const Field                        &F, 
						const IMatrix                      &A, 
						const FMatrix                     &Ap, 
						const VectorIn                     &b, 
						const Prime_Type                   &p,
						const size_t                        m,
						const size_t                        n)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap),
			  _F(F), 
			  _VDF(F), 
			  _res_p(b.size()), 
			  _digit_p(A.coldim()), 
			  _rand(F), 
			  _row(Ap.rowdim()), 
			  _col(Ap.coldim()), 
			  _m(m), 
			  _n(n), 
			  _U(m-1,Ap.rowdim()), 
			  _BMD(F) {
			
	
		
		
			for (size_t i=0;i<_m-1;++i)
				for (size_t j=0;j< _row;++j)
					_rand.random(_U.refEntry(i,j));
			
			Coefficient V(_col,n);
			for (size_t i=0;i< _col;++i)
				for (size_t j=0; j< n;++j)
					_rand.random(V.refEntry(i,j));
			
		       
			std::cout<<"U:\n";
			_U.write(std::cout, _F);

			std::cout<<"V:\n";
			V.write(std::cout, _F);

			Block UAp(_m, _row);

			typename Block::ConstRowIterator    iter_U   = _U.rowBegin();
			typename Block::RowIterator         iter_UAp = UAp.rowBegin();
			++iter_UAp;
			for (; iter_U != _U.rowEnd(); ++iter_UAp, ++iter_U) 
				Ap.applyTranspose( *iter_UAp , *iter_U );
				
			for (size_t i=0;i<m;++i)
				_rand.random(UAp.refEntry(0,i));			
			
			
			_Seq = new Sequence (&Ap, _F, UAp,V);
			std::cout<<"Sequence:\n";
			for (size_t i=0;i<_Seq->getRep().size();++i)
				_Seq->getRep()[i].write(std::cout,_F)<<"\n";
			std::cout<<"\n";

			
			_Dom = new BlockMasseyDomain<Field,Sequence> (_Seq);					      					      
			

#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
			ttMinPoly.clear();
#endif
		}

		virtual ~BlockWiedemannLiftingContainer() {
#ifdef _BM_TIMING
			_Dom->printTimer();
#endif
#ifdef _BBC_TIMING
			_Seq->printTimer();
#endif
			delete _Seq;
			delete _Dom;
		}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit,const IVector& residu) const {
			
			LinBox::integer tmp;
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, this->_R.convert(tmp,*iter));
			}
#ifdef RSTIMING			
			tGetDigitConvert.stop();
			ttGetDigitConvert+=tGetDigitConvert;
			tGetDigit.start();
#endif		

			std::cout<<"residue:\n";
			for (size_t i=0;i<_res_p.size();++i)
				_F.write(std::cout,_res_p[i])<<",";
			std::cout<<"\n";

			

			// compute the Minimal polynomial of the modified Sequence
			_Seq->setU(_res_p,0);
			_Seq->recompute();
			std::cout<<"Modified Sequence:\n";
			for (size_t i=0;i<_Seq->getRep().size();++i)
				_Seq->getRep()[i].write(std::cout,_F)<<"\n";
			std::cout<<"\n";
			
			FBlockPolynomial minpoly;
			std::vector<size_t> degree(_m);

#ifdef RSTIMING
			tMinPoly.start();
#endif
			_Dom->left_minpoly_rec(minpoly,degree);
#ifdef RSTIMING
			tMinPoly.stop();
			ttMinPoly+=tMinPoly;
#endif		
			std::cout<<"Block Minpoly:\n";
			for (size_t i=0;i<minpoly.size();++i)
				minpoly[i].write(std::cout,_F)<<"\n";			
			std::cout<<"\n";

			size_t idx=0;
			if ( _F.isZero(minpoly[0].getEntry(0,0))) {
				size_t i=1;
				while ( _F.isZero(minpoly[0].getEntry(i,0)))
					++i;
				if (i == _m)
					throw LinboxError(" block minpoly: matrix seems to be singular - abort");
				else 
					idx=i	;			
			}

			size_t deg = degree[idx];		
			BlasMatrix<Element> idx_poly(deg+1,_m-1);
			for (size_t i=0;i<deg+1;++i) 
				for (size_t j=0;j<_m-1;++j)
					idx_poly.setEntry(i,j,minpoly[i].getEntry(idx,j+1));

			BlasMatrix<Element> Combi(deg+1,_row);
			_BMD.mul(Combi,idx_poly,_U);
					

			FVector lhs(_col),row(_row);	
			for (size_t i=0;i<_row;++i)
				row[i]= Combi.getEntry(deg,i);
					
			_Ap.applyTranspose(lhs,row);					
			FVector lhsbis(lhs);
			for (int i = deg-1 ; i >= 0;--i) {
				for (size_t j=0;j<_row;++j)
					row[j]= Combi.getEntry(i,j);
				_VDF.add (lhs,row,lhsbis);
				_Ap.applyTranspose (lhsbis, lhs);			
			}   
					
			FVector accu (lhs);
			_Ap.applyTranspose(lhs,_res_p);
			_VDF.mulin(lhs,minpoly[deg].getEntry(idx,0));
			lhsbis=lhs;
			for (size_t i = deg-1 ; i > 0;--i) {
				_VDF.axpy (lhs,minpoly[i].getEntry(idx,0) , _res_p, lhsbis);
				_Ap.applyTranspose (lhsbis, lhs);			
			}  
					
			_VDF.addin(accu,lhs);
			Element scaling;
			_F.init(scaling);
			_F.neg(scaling,minpoly[0].getEntry(idx,0));
			_F.invin(scaling);
			_VDF.mul(_digit_p,accu,scaling);
			
			
			// check results
			FVector error(_Ap.coldim());
			_Ap.apply(error,_digit_p);		
			if (!_VDF.areEqual(error,_res_p)){
				std::cout<<"BlockMinpoly error\n";
				throw LinboxError("BlockMinpoly error\n");
			}
			
			
#ifdef RSTIMING			
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					this->_R.init(*iter, _F.convert(tmp,*iter_p));
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif
			return digit;
		}


	}; // end of class WiedemannLiftingContainerBase


	template <class _Ring, class _Field, class _IMatrix, class _FMatrix, class _Block>
	class BlockHankelLiftingContainer : public LiftingContainerBase< _Ring, _IMatrix> {

	public:
		typedef _Field                               Field;
		typedef _Ring                                 Ring;
		typedef _IMatrix                           IMatrix;		
		typedef _FMatrix                           FMatrix;	
		typedef typename Field::Element            Element;
		typedef typename IMatrix::Element          Integer;
		typedef std::vector<Integer>               IVector;
		typedef std::vector<Element>               FVector;
		typedef _Block                               Block;

	protected:
		
		const FMatrix&                      _Ap;
		const Diagonal<Field>               &_D;
		const BlockHankelInverse<Field>  &_Hinv;
		Field                                _F;
		mutable FVector                  _res_p;
		mutable FVector                _digit_p;		
		std::vector<std::vector<Element> >   _u;
		std::vector<std::vector<Element> >   _v;
		size_t                           _block;
		size_t                        _numblock;
		VectorDomain<Field>                 _VD;
		BlasMatrixDomain<Field>            _BMD;
		Element                           _zero;
		
	public:
#ifdef RSTIMING
		mutable Timer tGetDigit, ttGetDigit, tGetDigitConvert, ttGetDigitConvert;
#endif
		mutable Timer tApplyU, tApplyV, tApplyH, tAcc;;


		template <class Prime_Type, class VectorIn>
		BlockHankelLiftingContainer (const Ring&        R, 
					     const Field&       F, 
					     const IMatrix&     A, 
					     const FMatrix&    Ap,
					     const Diagonal<Field> &D,
					     const BlockHankelInverse<Field> &Hinv,
					     const Block&      U,
					     const Block&      V,
					     const VectorIn&   b, 
					     const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _Hinv(Hinv), _F(F),
			  _res_p(b.size()), _digit_p(A.coldim()),  _block(U.rowdim()), _numblock(A.coldim()/_block) , _VD(F), _BMD(F), _D(D)
		{
			tApplyU.clear();
			tApplyH.clear();
			tApplyV.clear();
			for (size_t i=0; i< _res_p.size(); ++i)
				_F.init(_res_p[i]);
			for (size_t i=0; i< _digit_p.size(); ++i)
				_F.init(_digit_p[i]);
			
			size_t block= U.rowdim();

			_u.resize(_block, std::vector<Element>(_numblock));
			_v.resize(_block, std::vector<Element>(_numblock));

			for (size_t i=0;i<_block;++i)
				for (size_t j=0;j<_numblock;++j){
					_F.assign(_u[i][j], U.getEntry(0, i*_numblock+j));
					_F.assign(_v[i][j], V.getEntry(i*_numblock+j, i));
				}
			_F.init(_zero,0);

			//Ap.write(std::cout,F);
#ifdef RSTIMING
			ttGetDigit.clear();
			ttGetDigitConvert.clear();
#endif	
#ifdef DEBUG_LC
			std::cout<<"Primes: ";
			_F.write(std::cout);
#endif

		}
		
		
		virtual ~BlockHankelLiftingContainer() {
#ifdef RSTIMING
			std::cout<<"time apply U: "<<tApplyU<<"\n";
			std::cout<<"time apply H: "<<tApplyH<<"\n";
			std::cout<<"time apply V: "<<tApplyV<<"\n";
#endif
		}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit, const IVector& residu) const {		
#ifdef RSTIMING			
			tGetDigitConvert.start();
#endif
			//LinBox::integer tmp;
			
			Hom<Ring, Field> hom(this->_R, _F);
			// res_p =  residu mod p
			//VectorHom::map (_res_p, residu, _F, this->_R);			
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					//_F. init (*iter_p, this->_R.convert(tmp,*iter));
					hom.image(*iter_p, *iter);
			}			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
			tGetDigit.start();
#endif			
		
			// compute the solution of _Ap^(-1).residu mod p = [V^T AV^T ... A^k]^T . Hinv . [U^T U^TA ... U^TA^k]^T residue mod p
			// with k= numblock -1
#ifdef RSTIMING		
			tAcc.clear();
			tAcc.start();			
#endif
			/*
			std::cout<<"b:=<";
			for (size_t i=0;i<_res_p.size()-1;++i)
				_F.write(std::cout,_res_p[i])<<",";
			_F.write(std::cout,_res_p[_res_p.size()-1])<<">;\n";
			*/

		

			size_t n = _Ap.coldim();
			// compute z0 = [U^T U^T Ap^T ... U^T Ap^k]^T . residue mod p
			FVector z0(n), b0(n), b1(n);
			_D.apply(b0, _res_p);
			_res_p=b0;
			BlasMatrix<Element> Apib(n, _numblock);
			for (size_t i=0;i<n;++i){
				_F.assign(Apib.refEntry(i,0), _res_p[i]);
			}
		
			int swi=1;
			for (size_t j=1; j<_numblock; ++j)
				if (swi){
					_Ap.apply(b1, b0);
					for (size_t i=0;i<n;++i)
						_F.assign(Apib.refEntry(i,j), b1[i]);					
					swi=0;
				}
				else{
					_Ap.apply(b0, b1);
					for (size_t i=0;i<n;++i)
						_F.assign(Apib.refEntry(i,j), b0[i]);
					swi=1;
				}
		
			FVector tmp(_numblock);
			for (size_t i=0; i<_block; ++i){
				BlasMatrix<Element> T(Apib, i*_numblock, 0, _numblock, _numblock); 				
				_BMD.mul(tmp, _u[i], T);
				for (size_t j=0;j<_numblock;++j){
					this->_F.assign(z0[j*_block+i], tmp[j]);					
				}
			}
#ifdef RSTIMING
			tAcc.stop();
			tApplyU+=tAcc;
			tAcc.clear();
			tAcc.start();
#endif
			// compute z1 = Hinv.z0
			FVector z1(n);
			_Hinv.apply(z1, z0);
#ifdef RSTIMING
			tAcc.stop();
			tApplyH+=tAcc;
			tAcc.clear();
			tAcc.start();
#endif
			/*
			std::cout<<" Hinv U b mod p done\n";
			std::cout<<"\n y:=<";
			for (size_t i=0;i<_digit_p.size()-1;++i)
				_F.write(std::cout,z1[i])<<",";
			_F.write(std::cout,z1[_digit_p.size()-1])<<">;\n";
			*/
			// compute digit_p  = [V^T AV^T ... A^k]^T.z1			
			FVector b_bar(n), b_hat(_numblock);			
			for (size_t i=0;i<n;++i)
				_F.assign(_digit_p[i], _zero);			
						
			for (int i= _numblock-1;i>=0; --i){				
				_Ap.apply(b1, _digit_p);
				_digit_p=b1;
				for (size_t j=0;j<_block;++j){
					_VD.mul(b_hat, _v[j], z1[i*_block+j]);
					for (size_t k=0;k<_numblock;++k)
						_F.assign(b_bar[j*_numblock+k], b_hat[k]);
				}
				_VD.addin(_digit_p, b_bar);
			}
#ifdef RSTIMING
			tAcc.stop();
			tApplyV+=tAcc;
#endif
			/*
			std::cout<<" V Hinv U b mod p done\n";
			std::cout<<"\n x:=<";
			for (size_t i=0;i<_digit_p.size()-1;++i)
				_F.write(std::cout,_digit_p[i])<<",";
			_F.write(std::cout,_digit_p[_digit_p.size()-1])<<">;\n";
			*/

#ifdef RSTIMING
			tGetDigit.stop();
			ttGetDigit+=tGetDigit;
			tGetDigitConvert.start();
#endif
			// digit = digit_p
			//VectorHom::map(digit, _digit_p, this->_R, _F);			
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					//this->_R.init(*iter, _F.convert(tmp,*iter_p));
					hom.preimage(*iter, *iter_p);
			}
			
#ifdef RSTIMING
			tGetDigitConvert.stop();
			ttGetDigitConvert += tGetDigitConvert;
#endif		
			return digit;
		}
		
	}; // end of class BlockHankelLiftingContainer


	
	/*
	 * SparseLULiftingContainer
	 */

	template <class _Ring, class _Field, class _IMatrix, class _FMatrix>
	class SparseLULiftingContainer : public LiftingContainerBase< _Ring, _IMatrix> {

	public:
		typedef _Field                               Field;
		typedef _Ring                                 Ring;
		typedef _IMatrix                           IMatrix;		
		typedef _FMatrix                           FMatrix;	
		typedef typename Field::Element            Element;
		typedef typename IMatrix::Element          Integer;
		typedef std::vector<Integer>               IVector;
		typedef std::vector<Element>               FVector;

	protected:
		
		const FMatrix&                       _L;
		const FMatrix&                       _U;
		const Permutation<_Field>&           _Q;
		const Permutation<_Field>&           _P;
		unsigned long                     _rank;
		Field                                _F;
		mutable FVector                  _res_p;
		mutable FVector                _digit_p;
		GaussDomain<Field>                  _GD;

		
	public:


		template <class Prime_Type, class VectorIn>
		SparseLULiftingContainer (const Ring&        R, 
					 const Field&       F, 
					 const IMatrix&     A,
					 const FMatrix&     L,
					 const Permutation<_Field>& Q,
					 const FMatrix&     U,
					 const Permutation<_Field>& P,
					 unsigned long   rank,
					 const VectorIn&    b, 
					 const Prime_Type&  p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _L(L),_Q(Q),_U(U), _P(P), _rank(rank),
			  _F(F), _res_p(b.size()), _digit_p(A.coldim()), _GD(F)
		{
			for (size_t i=0; i< _res_p.size(); ++i)
				_F.init(_res_p[i]);
			for (size_t i=0; i< _digit_p.size(); ++i)
				_F.init(_digit_p[i]);						
		}
		
		
		virtual ~SparseLULiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit, const IVector& residu) const {		

			// compute residu mod p
			Hom<Ring, Field> hom(this->_R, _F);
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					hom.image(*iter_p, *iter);
			}			
		
			// solve the system mod p using L.Q.U.P Factorization
			_GD.solve(_digit_p, _rank, _Q,_L,_U,_P, _res_p);			
		
			
			// promote new solution mod p to integers
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					//this->_R.init(*iter, _F.convert(tmp,*iter_p));
					hom.preimage(*iter, *iter_p);
			}

			return digit;
		}
		
	}; // end of class SparseLULiftingContainer



} // end of namespace LinBox
#endif //__LINBOX_lifting_container_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
