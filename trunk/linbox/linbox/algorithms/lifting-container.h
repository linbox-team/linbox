/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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

#ifndef _LIFTING_CONTAINER_H
#define _LIFTING_CONTAINER_H
//#define DEBUG_LC

#include <vector>

#include <linbox-config.h>
#include <linbox/util/debug.h>
#include <linbox/blackbox/apply.h>
#include <linbox/algorithms/blackbox-container.h>
#include <linbox/algorithms/massey-domain.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/compose.h>

namespace LinBox {

	// returns the square of the 'length of the longest column' of the matrix
	// here 'longest' means 'having greatest euclidean length'
	// so it is L(A) = Max_{j=1..N} Sum_{i=1..M} (A_ij)^2

	// this length function L(A) is pseudo-multiplicative in that
	//                L(AB) <= L(A)*L(B)*k   
	// where k is the number of columns of A (= the number of rows of B)
     
	// note: it is not actually a norm! but it is exactly what we need for lifting bound computation

	template < class Ring, class IMatrix>
	typename Ring::Element& NormBlackbox (const Ring& R, typename Ring::Element& norm, const IMatrix& A) {
		typedef typename Ring::Element Integer;
		Integer max,one,zero,sqsum;
		size_t m,n;
		n=A.coldim();
		m=A.rowdim();
		R.init(one,1);
		R.init(zero,0);
		typename std::vector<Integer>::const_iterator iter;
		std::vector<Integer> e(n,zero),tmp(m);
		max=zero;
		for (size_t i=0;i<n;i++){
			e[i]=one;
			A.apply(tmp,e);
			sqsum=zero;
			for (iter=tmp.begin();iter!=tmp.end();++iter)
				sqsum += (*iter)*(*iter);
			if (max < sqsum ) max=sqsum;
			e[i]=zero;
		}
		
		return norm= max;
	}

	template <class Ring, class Blackbox1, class Blackbox2>
	typename Ring::Element& NormBlackbox (const Ring& R, typename Ring::Element& norm, const Compose<Blackbox1,Blackbox2>& A) {
		typename Ring::Element maxL,maxR,result;
		NormBlackbox(R,maxL,*(A.getLeftPtr()));
		NormBlackbox(R,maxR,*(A.getRightPtr()));
				
		R.init(result, A.getLeftPtr()->coldim());
		R.mulin(result, maxR);
		R.mulin(result, maxL);
		return R.assign(norm, result);
	}
	
	
	template< class _Ring>
	class LiftingContainer {
		
	public:
		typedef _Ring Ring;
		typedef typename Ring::Element Integer;
				
		// return the length of container
		virtual const size_t length() const =0;

		// return the size of the solution
		virtual const int size() const = 0;

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

	protected:
						
		const IMatrix&            _A;
		Ring                      _R;			
		Integer                   _p;		
		IVector                   _b;
		VectorDomain<Ring>      _VDR;
		size_t               _length;
		Integer            _numbound;
		Integer            _denbound;

		BlasApply<Ring>          _BA;
	public:

		template <class Prime_Type, class Vector1>
		LiftingContainerBase (const Ring& R, const IMatrix& A, const Vector1& b, const Prime_Type& p):
			_A(A), _R(R), _VDR(R), _BA(R) {
			linbox_check(A.rowdim() == b.size());
			int n,m;
			n=A.rowdim();
			m=A.coldim();
#ifdef DEBUG_LC2
			cout<<"row: "<<n<<" col: "<<m<<endl;
#endif
			assert(m == n); //logic may not work otherwise
			// initialise the prime as an Integer
			_R.init(_p,p);
			
			// initialize res = b
			_b.resize(b.size());
			typename Vector1::const_iterator         b_iter    = b.begin();
			typename std::vector<Integer>::iterator  res_iter  = _b.begin() ;
			for (; b_iter != b.end(); ++res_iter, ++b_iter) 
				_R.init(*res_iter, *b_iter);
						
			// compute the length of container according to Hadamard's bound and Cramer's rules.
			Integer Anormsq;
						
			NormBlackbox(_R,Anormsq,A);
			
			typename std::vector<Integer>::const_iterator iterb = _b.begin();
			Integer bnormsq;
			_R.init(bnormsq, 0);
			for (;iterb!=_b.end();++iterb)
				bnormsq += (*iterb)*(*iterb);

			LinBox::integer normA,normb,N,D,L,prime;
			_R.convert(normA,Anormsq);
			_R.convert(normb,bnormsq);
			_R.convert(prime,_p);
			D = pow(sqrt(normA)+1, n);
			N = D * (sqrt(normb) + 1);
			L = (N * D * 8) / 3;          // the factor we need is (sqrt(5)+3)/2 ~= 2.616 < 8/3
			_length= logp(L,prime) + 1;   // round up instead of down
#ifdef DEBUG_LC                                   
			cout<<" norms computed, ||A||^2 = "<<normA<<", ||b||^2 = "<<normb<<", p = "<<_p<<"\n";
			cout<<" N = "<<N<<", D = "<<D<<", length = "<<_length<<"\n";
#endif
			_R.init(_numbound,N);
			_R.init(_denbound,D);
	
#ifdef DEBUG_LC		
			cout<<"lifting container initialized\n";			
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

				// prepare for computing next digit				
				// update tmp_integerv = _A * digit				
				IVector v2 (_lc._A.coldim());

				// v2 = _A.digit								
				_lc._BA.applyV (v2, _lc._A, digit);

				// update _res -= v2 
				_lc._VDR.subin (_res, v2);
				typename std::vector<Integer>::iterator p0;
				// update _res = _res / p
				int index=0;
				for ( p0 = _res.begin(); p0 != _res.end(); ++ p0, ++index){ 
					if (! _lc._R.isDivisor(*p0,_lc._p)) {
						cerr<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<"i: "<<endl;
						cout<<"residue "<<*p0<<" not divisible by modulus "<<_lc._p<<"i: "<<endl;
						return false;
					}
					_lc._R.divin(*p0, _lc._p);
				}

				_position++;
				return true;
			}

			bool operator != (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position != iterator._position;					
			}
	
			bool operator == (const const_iterator& iterator) const {
				if ( &_lc != &iterator._lc)
					assert("try to compare differents LiftingContainer::const_iterator , abort\n");
				return _position == iterator._position;					
			}

		};

				
		const_iterator begin() const {
			return const_iterator(*this);
		}		

		const_iterator end() const {
			return const_iterator (*this,_length);
		}
	
		virtual const size_t length() const
		{ return _length;}
		
		// return the size of the solution
		virtual const int size() const 
		{ return _A.coldim();}

		// return the ring
		virtual const Ring& ring() const 
		{ return _R; }
		
		// return the prime
		virtual const Integer& prime () const 
		{ return _p; }
		
		// return the bound for the numerator
		const Integer numbound() const
		{return _numbound;}
		
		// return the bound for the denominator
		const Integer denbound() const
		{return _denbound;}
		
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
		template <class Prime_Type, class VectorIn>
		DixonLiftingContainer (const Ring&       R, 
				       const Field&      F, 
				       const IMatrix&    A, 
				       const FMatrix&   Ap, 
				       const VectorIn&   b, 
				       const Prime_Type& p)
			: LiftingContainerBase<Ring,IMatrix> (R,A,b,p), _Ap(Ap), _F(F), _VDF(F), _res_p(b.size()), _digit_p(A.coldim()), _BA(F) {}
		
		
		virtual ~DixonLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit, const IVector& residu) const {
			
			
			LinBox::integer tmp;

			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, _R.convert(tmp,*iter));
			}
			
			// compute the solution by applying the inverse of A mod p
			_BA.applyV(_digit_p,_Ap,_res_p);

			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					_R.init(*iter, _F.convert(tmp,*iter_p));
			}
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
		}

		virtual ~WiedemannLiftingContainer() {}

		// return the field
		const Field& field() const { return _F; }

	protected:
		
		virtual IVector& nextdigit(IVector& digit,const IVector& residu) const {
			
			LinBox::integer tmp;

			// res_p =  residu mod p
			{
				typename FVector::iterator iter_p = _res_p.begin();
				typename IVector::const_iterator iter = residu.begin();
				for ( ;iter != residu. end(); ++iter, ++iter_p)
					_F. init (*iter_p, _R.convert(tmp,*iter));
			}
			
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
					cerr<<" the prime was not good \n, result will be wrong";
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
					nst=0;cerr<<"updating minpoly\n";
				}
				else {
					if (nst < nosolution_threshold) nst++;
					else{
						nosolution=true;
						throw PreconditionFailed (__FUNCTION__, __LINE__, "system is inconsistent or the choosen prime leads to inconsistent resolution");
					}
				}
			}
		
			// digit = digit_p
			{
				typename FVector::const_iterator iter_p = _digit_p.begin(); 
				typename IVector::iterator iter = digit.begin();
				for ( ; iter_p!= _digit_p.end(); ++iter_p, ++iter)
					_R.init(*iter, _F.convert(tmp,*iter_p));
			}
			
			return digit;
		}


	}; // end of class WiedemannLiftingContainerBase




    } // end of namespace LinBox
#endif
