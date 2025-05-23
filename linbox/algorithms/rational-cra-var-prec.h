/* linbox/blackbox/rational-reconstruction-base.h
 * Copyright (C) 2009 Anna Marszalek
 *
 * Written by Anna Marszalek <aniau@astronet.pl>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_rational2_cra_H
#define __LINBOX_rational2_cra_H


#include "givaro/zring.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/algorithms/classic-rational-reconstruction.h"

namespace LinBox
{
	/**
     * \brief Chinese remainder of vector of rationals.
     *
     * VarPrec: Variable preconditioner = lift random combination of residues
	 *
	 * Compute the reconstruction of rational numbers
	 * Either by Early Termination see [Dumas, Saunder, Villard, JSC 32 (1/2), pp 71-99, 2001],
	 * Or via a bound on the size of the Integers.
	 */
	//typedef Givaro::ZRing<Integer> Integers;
	//typedef Integers::Element Integer;

	template<class RatCRABase, class RatRecon = RReconstruction<Givaro::ZRing<Integer>, ClassicMaxQRationalReconstruction<Givaro::ZRing<Integer> > > >
	struct RationalChineseRemainderVarPrec {


		typedef typename RatCRABase::Domain		Domain;
		typedef typename RatCRABase::DomainElement	DomainElement;
	protected:
		RatCRABase Builder_;
		RatRecon RR_;
	public:

		int IterCounter;

		template<class Param>
		RationalChineseRemainderVarPrec(const Param& b, const RatRecon& RR = RatRecon()) :
			Builder_(b), RR_(RR)
		{
			IterCounter = 0;
		}

		RationalChineseRemainderVarPrec(RatCRABase b, const RatRecon& RR = RatRecon()) :
			Builder_(b), RR_()
		{
			IterCounter = 0;
		}

		/** \brief The Rational CRA loop

		  Given a function to generate residues mod a single prime,
		  this loop produces the residues resulting from the Chinese
		  remainder process on sufficiently many primes to meet the
		  termination condition.

		  \param Iteration  Function object of two arguments, <code>Iteration(r,
		  p)</code>, given prime \p p it outputs residue(s) \p r.  This
		  loop may be parallelized.  \p Iteration must be reentrant, thread
		  safe.  For example, \p Iteration may be returning the coefficients of
		  the minimal polynomial of a matrix \c mod \p p.

		  @warning We won't detect bad primes.

		  \param genprime  RandIter object for generating primes.
		  \param[out] num  the rational numerator
		  \param[out] den  the rational denominator
		  */
		template<class Function, class RandPrimeIterator>
		Integer & operator() (Integer& num, Integer& den
				      , Function& Iteration, RandPrimeIterator& genprime)
		{
			{
				++genprime;
				Domain D(*genprime);
				DomainElement r; D.init(r);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.initialize( D, r );
				++IterCounter;
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			Integer f_in,m_in;
			Builder_.getPreconditioner(f_in,m_in);

			while( ! Builder_.terminated() ) {

				++genprime;
				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						Integer res; return Builder_.result(res);
					}
				}
				coprime = 0;
				Domain D(*genprime);
				DomainElement r; D.init(r);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );
				if (RR_.scheduled((size_t)IterCounter-1)) {
					Integer Mint ; Builder_.getModulus(Mint);
					Integer rint ; Builder_.getResidue(rint);
					if (RR_.RationalReconstruction(num,den,rint,Mint)) {
						Builder_.changePreconditioner(f_in*num,m_in*den);
						int k ; Builder_.getThreshold(k);
						if (this->operator()(k,num,den,Iteration,genprime)) break;
						else {
							Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
						}
					}
				}
				++IterCounter;
			}
			Integer g;
			Builder_.result(num,den);
			if (gcd(g,num,den) != 1) { num /=g; den/=g; }
			return num; //Builder_.result(num, den);
		}

		/*
		 * progress for k>=0 iterations
		 * run until terminated() if k<0
		 * no rational reconstruction!
		 */

		template<class Function, class RandPrimeIterator>
		bool operator() (const int k, Integer& num, Integer& den
				 , Function& Iteration, RandPrimeIterator& genprime)
		{

			if ((IterCounter==0) && (k != 0)) {
				++IterCounter;
				++genprime;
				Domain D(*genprime);
				DomainElement r; D.init(r);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.initialize( D, r );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			Integer f_in,m_in;
			Builder_.getPreconditioner(f_in,m_in);

			for ( int i=0; ((k< 0) && ! Builder_.terminated()) || (i<k) ; ++i ) {

				++genprime;

				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						//Integer res; return Builder_.result(res);
						return false;
					}
				}
				coprime = 0;
				Domain D(*genprime);
				DomainElement r; D.init(r);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );
				//if (RR_.scheduled(IterCounter-1)) {
				++IterCounter;
#if 0

				Integer M ; Builder_.getModulus(M);

				Integer r ; Builder_.getResidue(r);

				if (RR_.RationalReconstruction(num,den,r,M)) {

					if (Builder_.changePreconditioner(f_in*num,m_in*den)) break;

					else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results


				}
#endif

				//}
			}

			if (Builder_.terminated() ) {
				Integer g;
				//Builder_.getPreconditioner(p_div,p_mul);
				Builder_.result(num,den);
				//num = p_div*res;
				//den = p_mul;
				if (gcd(g,num,den) != 1) { num /=g; den/=g; }
				return true;
			}
			else return false;
		}

		template <template<class,class> class Vect, template<class> class Alloc,  class Function, class RandPrimeIterator>
		Vect<Integer,Alloc<Integer> > & operator() (Vect<Integer,Alloc<Integer> >& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime)
		{
            typedef Vect<DomainElement,Alloc<DomainElement> > DomVect;
			{
				++IterCounter;
				++genprime;
				Domain D(*genprime);
				DomVect r;
                Iteration(r, D);  // FIXME bad primes ignored
                Builder_.initialize( D, r );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			Vect<Integer, Alloc<Integer>> f_in,m_in;
			Builder_.getPreconditioner(f_in,m_in);

			//while( ! Builder_.terminated() )
			while (1) { // in case of terminated() - checks for RR of the whole vector
				//++IterCounter;
				++genprime;
				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return num;
					}
				}
				coprime = 0;


				Domain D(*genprime);
                DomVect r;
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );

				if (RR_.scheduled((size_t)IterCounter-1) || Builder_.terminated()) {
					Integer Mint ; Builder_.getModulus(Mint);
					if ( Builder_.terminated() ) {//early or full termination occurred, check reconstruction of the whole vector
						//early or full termination
						Vect<Integer,Alloc<Integer> > r_v ;
						Builder_.getResidue(r_v);
						if (RR_.RationalReconstruction(num,den,r_v,Mint) ) {
							Vect<Integer,Alloc<Integer> > vnum(num),vden(m_in.size(),den);
							for (int i=0; i < (int)vnum.size(); ++ i) {
								if (vnum[(size_t)i]==0) vnum[(size_t)i] = 1; // no prec
							}
							Builder_.productin(vnum, f_in);
							Builder_.productin(vden,m_in);
							Builder_.changePreconditioner(vnum,vden) ;
							int k ;
							Builder_.getThreshold(k);
							if (this->operator()(k,num,den,Iteration,genprime)) {
								break;
							}
							else {	// back to original preconditioners
								Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
								Builder_.changeVector();
							}
						}
						else {  //back to original preconditioners
							Builder_.changePreconditioner(f_in,m_in);
							Builder_.changeVector();
						}
					}
					else {
						//heuristics: reconstruction of vector
						Integer rint ;
						Builder_.getResidue(rint);
						Integer n,d;
						if (RR_.RationalReconstruction(n,d,rint,Mint)) {
							Vect<Integer,Alloc<Integer> > vden(m_in.size(),d);
							Builder_.productin(vden,m_in);
							Builder_.changePreconditioner(f_in,vden);
							int k; Builder_.getThreshold(k);
							if (this->operator()(k,num,den,Iteration,genprime)) { //prob. certify result of RR
								m_in = vden;
							}
							else {	//false result of RR
								Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
							}

						}
					}
				}
				++IterCounter;
			}
			Builder_.result(num,den);

			return num;
		}

		template<class Function, class RandPrimeIterator>
		BlasVector<Givaro::ZRing<Integer> > & operator() (BlasVector<Givaro::ZRing<Integer> >& num, Integer& den
						      , Function& Iteration, RandPrimeIterator& genprime)
		{
			{
				++IterCounter;
				++genprime;
				Domain D(*genprime);
				BlasVector<Domain> r(D);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.initialize( D, r );
			}

			int coprime =0;
			int maxnoncoprime = 1000;

			Givaro::ZRing<Integer> Z;
			BlasVector<Givaro::ZRing<Integer> > f_in(Z),m_in(Z);
			Builder_.getPreconditioner(f_in,m_in);

			//while( ! Builder_.terminated() )
			while (1) { // in case of terminated() - checks for RR of the whole vector
				//++IterCounter;
				++genprime;
				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return num;
					}
				}
				coprime = 0;


				Domain D(*genprime);
				BlasVector<Domain > r(D);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );

				if (RR_.scheduled((size_t)IterCounter-1) || Builder_.terminated()) {
					Integer Mint ; Builder_.getModulus(Mint);
					if ( Builder_.terminated() ) {//early or full termination occurred, check reconstruction of the whole vector
						//early or full termination
						BlasVector<Givaro::ZRing<Integer> > r_v(Z) ;
						Builder_.getResidue(r_v);
						if (RR_.RationalReconstruction(num,den,r_v,Mint) ) {
							BlasVector<Givaro::ZRing<Integer> > vnum(num),vden(Z,m_in.size(),den);
							for (int i=0; i < (int)vnum.size(); ++ i) {
								if (vnum[(size_t)i]==0) vnum[(size_t)i] = 1; // no prec
							}
							Builder_.productin(vnum, f_in);
							Builder_.productin(vden,m_in);
							Builder_.changePreconditioner(vnum,vden) ;
							int k ;
							Builder_.getThreshold(k);
							if (this->operator()(k,num,den,Iteration,genprime)) {
								break;
							}
							else {	// back to original preconditioners
								Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
								Builder_.changeVector();
							}
						}
						else {  //back to original preconditioners
							Builder_.changePreconditioner(f_in,m_in);
							Builder_.changeVector();
						}
					}
					else {
						//heuristics: reconstruction of vector
						Integer rint ;
						Builder_.getResidue(rint);
						Integer n,d;
						if (RR_.RationalReconstruction(n,d,rint,Mint)) {
							BlasVector<Givaro::ZRing<Integer> > vden(Z,m_in.size(),d);
							Builder_.productin(vden,m_in);
							Builder_.changePreconditioner(f_in,vden);
							int k; Builder_.getThreshold(k);
							if (this->operator()(k,num,den,Iteration,genprime)) { //prob. certify result of RR
								m_in = vden;
							}
							else {	//false result of RR
								Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
							}

						}
					}
				}
				++IterCounter;
			}
			Builder_.result(num,den);

			return num;
		}


		/*
		 * progress for k>=0 iterations
		 * run until terminated if k <0
		 */
		template<template <class, class> class Vect, template<class> class Alloc, class Function, class RandPrimeIterator>
		bool operator() (const int k, Vect<Integer, Alloc<Integer>  >& num, Integer& den
				 , Function& Iteration, RandPrimeIterator& genprime)
		{
			if ((IterCounter==0) && (k != 0)) {
				++IterCounter;
				++genprime;
				Domain D(*genprime);
				Vect<DomainElement, Alloc<DomainElement>  > r;
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.initialize( D, r );
			}
			int coprime =0;
			int maxnoncoprime = 1000;

			Vect<Integer, Alloc<Integer>  > f_in,m_in;
			Builder_.getPreconditioner(f_in,m_in);
			for (int i=0; ((k<0) && Builder_.terminated()) || (i <k); ++i ) {
				//++IterCounter;
				++genprime;
				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return false;
						//return Builder_.result(res);
					}
				}
				coprime = 0;
				Domain D(*genprime);
				Vect<DomainElement, Alloc<DomainElement>  > r;
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );
				//if (RR_.scheduled(IterCounter-1))
				++IterCounter;

#if 0
				Integer M ; Builder_.getModulus(M);
				if ( Builder_.terminated() ) {
					Vect<Integer> r ; Builder_.getResidue(r);
					if (RR_.RationalReconstruction(num,den,r,M) ) {
						Vect<Integer> vnum(num),vden(m_in.size(),den);
						Builder_.productin(vnum, f_in); Builder_.productin(vden,m_in);
						if (Builder_.changePreconditioner(vnum,vden)) break;
						else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
					}
				}
				else {
					Integer r ; Builder_.getResidue(r);
					Integer n,d;
					if (RR_.RationalReconstruction(n,d,r,M)) {
						Vect<Integer > vden(m_in.size(),d);
						Builder_.productin(vden,m_in);
						if (Builder_.changePreconditioner(f_in,vden)) break;
						else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
					}
				}
#endif
			}
			if (Builder_.terminated()) {
				//Vect<Integer> p_mul, p_div,res, g;
				//Builder_.getPreconditioner(p_div,p_mul);
				Builder_.result(num,den);
#if 0
				num = Builder_productin(res,p_div);
				typename Vect<Integer>::iterator itnum,itden,ittmp;
				den = 1; Integer denold = 1;
				for (itnum = num.begin(), itden=p_mul.begin(); itnum != num.end(); ++itnum,++itden) {
					lcm(den,den,*itden);
					if (denold != den) {
						Integer h = den/denold;
						ittmp = num.begin();
						for (;itnum != itnum; ++ittmp)  *ittmp *= h;
					}
					denold = den;
				}
				den = p_mul;
#endif

				return true;
			}
			else return false;
		}


		// temp. spec. for BlasVector
		template<class Function, class RandPrimeIterator>
		bool operator() (const int k, BlasVector<Givaro::ZRing<Integer>  >& num
				 , Integer& den, Function& Iteration, RandPrimeIterator& genprime)
		{
			if ((IterCounter==0) && (k != 0)) {
				++IterCounter;
				++genprime;
				Domain D(*genprime);
				BlasVector<Domain> r(D);
				Builder_.initialize( D, Iteration(r, D) );
			}
			int coprime =0;
			int maxnoncoprime = 1000;
			Givaro::ZRing<Integer> Z;

			BlasVector<Givaro::ZRing<Integer> > f_in(Z),m_in(Z);
			Builder_.getPreconditioner(f_in,m_in);
			for (int i=0; ((k<0) && Builder_.terminated()) || (i <k); ++i ) {
				//++IterCounter;
				++genprime;
				while(Builder_.noncoprime(*genprime) ) {
					++genprime;
					++coprime;
					if (coprime > maxnoncoprime) {
						std::cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
						return false;
						//return Builder_.result(res);
					}
				}
				coprime = 0;
				Domain D(*genprime);
				BlasVector<Domain > r(D);
				Iteration(r, D); // FIXME bad primes ignored
				Builder_.progress( D, r );
				//if (RR_.scheduled(IterCounter-1))
				++IterCounter;

#if 0
				Integer M ; Builder_.getModulus(M);
				if ( Builder_.terminated() ) {
					BlasVector<Givaro::ZRing<Integer> > r ; Builder_.getResidue(r);
					if (RR_.RationalReconstruction(num,den,r,M) ) {
						BlasVector<Givaro::ZRing<Integer> > vnum(num),vden(m_in.size(),den);
						Builder_.productin(vnum, f_in); Builder_.productin(vden,m_in);
						if (Builder_.changePreconditioner(vnum,vden)) break;
						else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
					}
				}
				else {
					Integer r ; Builder_.getResidue(r);
					Integer n,d;
					if (RR_.RationalReconstruction(n,d,r,M)) {
						BlasVector<Givaro::ZRing<Integer> > vden(m_in.size(),d);
						Builder_.productin(vden,m_in);
						if (Builder_.changePreconditioner(f_in,vden)) break;
						else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
					}
				}
#endif
			}
			if (Builder_.terminated()) {
				//BlasVector<Givaro::ZRing<Integer> > p_mul, p_div,res, g;
				//Builder_.getPreconditioner(p_div,p_mul);
				Builder_.result(num,den);
#if 0
				num = Builder_productin(res,p_div);
				typename BlasVector<Givaro::ZRing<Integer> >::iterator itnum,itden,ittmp;
				den = 1; Integer denold = 1;
				for (itnum = num.begin(), itden=p_mul.begin(); itnum != num.end(); ++itnum,++itden) {
					lcm(den,den,*itden);
					if (denold != den) {
						Integer h = den/denold;
						ittmp = num.begin();
						for (;itnum != itnum; ++ittmp)  *ittmp *= h;
					}
					denold = den;
				}
				den = p_mul;
#endif

				return true;
			}
			else return false;
		}


#ifdef __LB_CRA_TIMING__
		std::ostream& reportTimes(std::ostream& os)
		{
			//Builder_.reportTimes(os);
			return os <<  "Iterations:" << IterCounter << "\n" ;
		}
#endif
	};

#ifdef _LB_RCRATIMING

	class RCRATimer {
	public:
		mutable Timer ttInit, tt RRecon, ttIRecon, ttImaging, ttIteration, ttOther;
		void clear() const
		{
			ttInit.clear();
			ttRRecon.clear();
			ttIRecon.clear();
			ttImaging.clear();
			ttIteration.clear();
			ttother.clear();
		}
#if 0
		template<class RR, class LC>
		void update(RR& rr, LC& lc) const
		{
			ttSetup += lc.ttSetup;
			ttRecon += rr.ttRecon;
			ttGetDigit += lc.ttGetDigit;
			ttGetDigitConvert += lc.ttGetDigitConvert;
			ttRingOther += lc.ttRingOther;
			ttRingApply += lc.ttRingApply;
		}
#endif
	};
#endif

}

#endif // __LINBOX_rational2_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
