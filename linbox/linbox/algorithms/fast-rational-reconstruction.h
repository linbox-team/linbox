/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/classic-rational-reconstruction.h
 * Copyright (C) 2009 Anna Marszalek
 *
 * Written by Anna Marszalek <aniau@astronet.pl>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_fast_reconstruction_H
#define __LINBOX_fast_reconstruction_H

#define __FASTRR_DEFAULT_THRESHOLD 384

#include <iostream>
#include "linbox/algorithms/rational-reconstruction-base.h"
#include "linbox/util/timer.h"


namespace LinBox
{

	/*
	 * implements fast rational reconstruction by Wan & Pan algorithm [Wan & Pan 2002],
	 * Wang's bounds [Wang 1981] are used as default
	 */

	template <class Ring>
	class FastRationalReconstruction: public RReconstructionBase<Ring> {
	protected:
		size_t _threshold;
	public:
		const Ring _intRing;
		typedef typename Ring::Element Element;

		FastRationalReconstruction(const Ring& Z) :
			RReconstructionBase<Ring>(Z), _intRing(Z)
		{
			_threshold = __FASTRR_DEFAULT_THRESHOLD;
			if (_threshold < 2) _threshold = 2;
		}

		~FastRationalReconstruction() {}

		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) const
		{
			Element a_bound; _intRing.sqrt(a_bound, m/2);
			reconstructRational(a,b,x,m,a_bound);
			return (a < a_bound);
		}

		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound ) const
		{
			if (x < a_bound) {
				//case (i)
				a=x;
				b=1;
				return true;
			}
			else if (m-x < a_bound) {
				a = x-m;
				b = 1;
				return true;
			}

			Element bound = a_bound << 1;

			if (m/bound > 1) {
				fastReconstructRational(a,b,x,m,m/bound);

				if (_intRing.abs(a) < a_bound) {
					return true;
				}
				return false;
			}
			else {
				//either case (i) or false
				return false;
			}
		}

	protected:
		/*
		 * using mutable can be messiy, probabily need better solution
		 */
		mutable Element cur_ri;
		mutable Element cur_rinext;
		mutable Element cur_ainext;
		mutable Element cur_qinext;

		Element& powtwo(Element& h, const Element& log_h) const
		{
			h = 1;
			if (log_h <= 0) return h;
			if (log_h < ULONG_MAX) {
				h<<=log_h;
				return h;
			}
			else {
				Element n,m;
				quoRem(n,m,log_h,(Element)(ULONG_MAX-1));
				for (int i=0; i < n; ++i) {
					h <<=(long int)(ULONG_MAX-1);
				}
				h <= (long int)m;
				return h;
			}
		}

		Element& powtwo(Element& h, const size_t log_h) const
		{
			h = 1;
			if (log_h <= 0) return h;
			h<<=log_h;
			return h;
		}

		bool fastReconstructRational(Element& n, Element& d, const Element& x, const Element& m, const Element& d_bound ) const
		{

			size_t log_m = m.bitsize()-1;
			size_t log_bound = d_bound.bitsize()-1;
			Element ai, bi, ci, di;
			ai=1;bi=0;ci=0;di=1;
			Element bound;  _intRing.powtwo(bound, log_bound);

			cur_ri = m;
			cur_rinext = x;
			_intRing.quo(cur_ainext, cur_ri, cur_rinext);
			cur_qinext = cur_ainext;//  ri/ri_next

			if (!fastEEA (ai,bi,ci,di,m,log_m,x,bound, log_bound)) return false;

			int K=0;

			Element init_ainext = cur_ainext;
			Element init_qinext = cur_qinext;
			Element init_rinext = cur_rinext;
			Element init_ri = cur_ri;
			/* correction if $K <> bound = 2^log_bound */
			if (cur_rinext > 0) {
				while (cur_ainext <= d_bound) {
					++K;
					Element tmp;
					tmp = ai;
					ai = cur_ainext;
					bi = tmp;
					tmp = ci;
					_intRing.axpy(ci, cur_qinext,tmp, di);
					++this->C.mul_counter;
					di = tmp;

					tmp = cur_rinext;
					_intRing.axpy(cur_rinext, -cur_qinext, tmp, cur_ri);
					cur_ri = tmp;

					if (cur_rinext==0 ) {
						cur_qinext = m+1;
						cur_ainext = m+1;
						break;
					}

					++this->C.mul_counter;
					_intRing.quo(cur_qinext, cur_rinext, cur_ri);
					++this->C.div_counter;
					_intRing.axpy(cur_ainext, ai, cur_qinext, bi);
					++this->C.mul_counter;
				}
			}

			_intRing.mul(n,x, ai);
			_intRing.maxpyin(n,m,ci);
			//n = x_in*ai-m*ci;
			d = ai;
			return true;

		}

		void prevEEA(Element& aprev, Element& bprev, Element& cprev, Element& dprev,
			     const Element& ai, const Element& bi, const Element& ci, const Element& di) const
		{
			aprev = bi;
			cprev = di;
			if ((bi==1) && (di==1)) {
				bprev = 1; dprev = 0;
				Element tmp = cur_ri;
				cur_ri = cur_ri+cur_rinext;
				cur_rinext = tmp;
				cur_ainext = ai;
				cur_qinext = 1;
				return;
				//before last matrix we know next call to prevEEA is Id
			}
			Element qi;
			_intRing.quo(qi, ai, bi);
			++this->C.div_counter;

			_intRing.axpy(bprev, -qi, bi,ai);
			_intRing.axpy(dprev, -qi, di,ci);
			++this->C.mul_counter;++this->C.mul_counter;
			Element tmp;
			tmp = cur_ri;
			_intRing.axpy(cur_ri, tmp, qi, cur_rinext);
			//cur_ri = cur_ri *qi + cur_rinext;
			cur_rinext = tmp;
			cur_ainext = ai;
			cur_qinext = qi;
			++this->C.mul_counter;

		}

		/* extended Euclidean Algorithm */
		bool classicEEA(Element& ai, Element& bi, Element& ci, Element& di, const Element& r0, const Element& r1, const Element& bound, int K =0) const
		{

			Element ri, rinext;
			ri = r0; rinext = r1;
			ai =di = 1;
			bi =ci = 0;

			Element ainext,cinext;
			Element qinext;
			_intRing.quo(qinext,ri,rinext);
			++this->C.div_counter;
			while (1) {
				++K;
				//ainext = ai*qinext + bi;
				_intRing.axpy(ainext, ai,qinext,bi);
				++this->C.mul_counter;
				if (bound < ainext) {
					cur_ri = ri;
					cur_rinext = rinext;
					cur_ainext = ainext;
					cur_qinext = qinext;
					return true;
				}
				else {
					bi = ai;
					ai = ainext;
					Element temp = ci;
					_intRing.axpy(ci, temp,qinext,di);
					++this->C.mul_counter;
					di = temp;

					temp = ri;
					ri=rinext;
					_intRing.axpy(rinext, -qinext,ri,temp);//rinext = ri-qinext*rinext
					++this->C.mul_counter;
					if (rinext==0) {
						//ainext = infinity
						cur_ri = ri;
						cur_rinext = rinext;
						cur_ainext = r0+1;//infinity
						cur_qinext = cur_ainext;//infinity
						return true;
					}
					_intRing.quo(qinext,ri,rinext);
					++this->C.div_counter;
				}
			}
		}

		/* log(m)-1 <= d < log(m) ; d >=1
		 * 2^h=powh=bound>=2, h <= d
		 */
		bool fastEEA(Element& ai,Element& bi, Element& ci, Element& di, const Element& m, const size_t d, const Element& n, const Element& powh, const size_t& h) const
		{
			ai=Element(1);
			di=Element(1);
			bi=Element(0);
			ci=Element(0); //Q(0)=Id

			if (powh > m) std::cerr << "should not happen m < powh" << std::flush;
			if (m==n) {
				cur_ri = m;
				cur_rinext = 0;
				cur_ainext = m+1;//infinity
				cur_qinext = m+1;
				ai=1; bi=1; ci=1; di=0;
				return true;
			}

			if (powh < 1) {
				std::cerr << "wrong powh used"<< std::flush;
				cur_ri=m, cur_rinext=n;cur_ainext=m+1; cur_qinext = m+1;
				return false; //should not happen
			}
			if (n < 2) {
				//n==1 -> Q1=(m,1//1,0)
				//n==0 -> Q1=infinity
				if (n==0) std::cerr << "should not happen n=0" << std::flush;
				cur_ri = m;
				cur_rinext = 1;
				cur_ainext = m;
				cur_qinext = m;
				if (cur_ainext > powh) return true;
				else {
					cur_ri = 1;
					cur_rinext = 0;
					cur_ainext = m+1;
					cur_qinext = m+1;
				}
				return true;
			}

			if (h < 1) {
				//if (qinext==1) { return (1,1,1,0) or (1,0,0,1) }
				if ((n << 1) > m) {
					cur_ri = n;
					cur_rinext = m-n;
					_intRing.quo(cur_qinext, cur_ri, cur_rinext);
					++this->C.div_counter;
					cur_ainext = cur_qinext + 1;
					ai=1; bi=1; ci=1; di=0;
				}
				else {
					cur_ri = m;
					cur_rinext = n;
					_intRing.quo(cur_ainext,m,n);
					++this->C.div_counter;
					cur_qinext = cur_ainext;
				}
				return true;
			}

			if (n.bitsize() < _threshold) {       //what about m?
				return classicEEA(ai,bi,ci,di,m,n,powh);
			}

			// size_t log_n = n.bitsize()-1;

			if (2*h+1 < d) {
				size_t lambda = d-2*h-1;
				Element aistar, bistar, cistar, distar;
				aistar=distar=1;
				bistar=cistar=0;

				Element mstar = m >> (long unsigned int) lambda;
				Element nstar = n >> (long unsigned int) lambda;

				size_t log_mstar = 2*h+1;
				if (nstar > 0) if (!fastEEA(aistar, bistar, cistar, distar, mstar, log_mstar,nstar, powh, h)) return false;

				int K=2;//reatreat steps;

				if ((aistar > 1) && (distar > 0)) { //we have to go back 2 steps
					Element aprev,bprev,cprev,dprev;
					for (int i=0; i < K; ++i) {
						prevEEA(aprev,bprev,cprev,dprev,aistar,bistar,cistar,distar);
						aistar=aprev;bistar=bprev;cistar=cprev;distar=dprev;
					}
					ai = aistar; bi =bistar; ci = cistar; di= distar;
				}

				_intRing.mul(cur_ri,m,di);
				_intRing.maxpyin(cur_ri,n,bi);
				_intRing.mul(cur_rinext,n,ai);
				_intRing.maxpyin(cur_rinext,m,ci);
				this->C.mul_counter+=4;

				if (cur_ri < 0) {
					cur_ri = -cur_ri;
					cur_rinext = -cur_rinext;
				}
				if (cur_rinext>0) {
					_intRing.quo(cur_qinext,cur_ri,cur_rinext);
					++this->C.div_counter;
					_intRing.axpy(cur_ainext, ai,cur_qinext,bi);
					++this->C.mul_counter;
				}
				else {
					cur_ainext = m+1;//infinity
					cur_qinext = cur_ainext;
				}
			}
			else { //if (h <= d-1)  // modification of Wan&Pan
				Element a1,a2,b1,b2,c1,c2,d1,d2;
				a1=a2=d1=d2=1;
				b1=b2=c1=c2=0;

				Element sqrth;
				size_t logsqrth;
				logsqrth = h >> (int) 1;
				powtwo(sqrth, logsqrth);

				if (!fastEEA(a1,b1,c1,d1,m,d,n,sqrth, logsqrth)) return false;

				Element ri = cur_ri;
				Element rinext = cur_rinext;

				size_t log_m;

				if ((rinext > 0) && (cur_ainext <= powh)){
					log_m = rinext.bitsize()-1;
					Element m2, n2;
					m2 = rinext;
					_intRing.axpy(n2, -cur_qinext, rinext, ri);
					++this->C.mul_counter;

					/* compute Q(i+1) */
					Element tmp = a1;
					a1 = cur_ainext;
					b1 = tmp;
					tmp = c1;
					_intRing.axpy (tmp, cur_qinext, c1, d1);
					++this->C.mul_counter;
					d1 = c1;
					c1 = tmp;

					int k = (int)a1.bitsize()-1 ;
					int _k;
					if (h-k > 2)
						_k = (int)(h-k-2);
					else _k = 0;
					if (n2 >0) {
						if (a1 < powh) {
							if (!fastEEA(a2,b2,c2,d2,m2,log_m,n2, powtwo(sqrth,_k), _k)) return false;
						}
						else {
							ai = a1; bi = b1; ci=c1; di = d1;
							cur_ri = m2;
							cur_rinext = n2;
							_intRing.quo(cur_qinext,m2,n2);
							++this->C.div_counter;
							_intRing.axpy(cur_ainext,a1,cur_qinext,b1);
							++this->C.mul_counter;
							return true;
						}
					}
					else {
						ai = a1; bi = b1; ci=c1; di = d1;
						cur_ri = m2;
						cur_rinext = n2;
						cur_qinext = m+1;
						cur_ainext = m+1;
						return true;
					}
				}
				else {//ri_next == 0 || cur_ainext >powh
					ai = a1; bi = b1; ci=c1; di = d1;
					if (cur_rinext<=0) {
						cur_ainext = m +1;
						cur_qinext = cur_ainext;
					}
					return true;
				}

				_intRing.mul(ai,b1,c2);
				_intRing.axpyin(ai,a1,a2);
				//aistar = a1*a2 + b1*c2;
				_intRing.mul(bi,b1,d2);
				_intRing.axpyin(bi, a1,b2);
				//bistar = a1*b2 + b1*d2;
				_intRing.mul(ci,d1,c2);
				_intRing.axpyin(ci, c1,a2);
				//cistar = c1*a2 + d1*c2;
				_intRing.mul(di,d1,d2);
				_intRing.axpyin(di, c1,b2);
				//distar = c1*b2 + d1*d2;
				this->C.mul_counter+=8;

				_intRing.mul(cur_ri,m,di);
				_intRing.maxpyin(cur_ri,n,bi);
				_intRing.mul(cur_rinext,n,ai);
				_intRing.maxpyin(cur_rinext,m,ci);
				this->C.mul_counter+=4;
				//ri = m*di - n*bi;
				//rinext = -m*ci + n*ai;
				if (cur_ri < 0) {
					cur_ri = -cur_ri;
					cur_rinext = -cur_rinext;
				}
				if (cur_rinext>0) {
					_intRing.quo(cur_qinext,cur_ri,cur_rinext);
					++this->C.div_counter;
					_intRing.axpy(cur_ainext, ai,cur_qinext,bi);
					++this->C.mul_counter;
				}
				else {
					cur_ainext = m+1;//infinity
					cur_qinext = cur_ainext;
					return true;
				}


				Element aprev,bprev,cprev,dprev;
				aprev=dprev=1;
				bprev=cprev=0;

				if (ai > powh) {//at most 2 forward steps, at most 2 backward steps
					//backward loop (max 0 steps)
					int K=0;
					while (ai > powh) {//one step back
						++K;
						prevEEA(aprev,bprev,cprev,dprev,ai,bi,ci,di);
						ai=aprev;bi=bprev;ci=cprev;di=dprev;
					}
					std::cerr << "Error: " << K << " backward steps in step 2\n";
					std::cerr << "->End:" << cur_ri << " " << cur_rinext << " " <<ai<<"/"<< cur_ainext << "\n"<< std::flush;
					return false;
				}
				/* //modification of Wan&Pan
				   }
				   else {//h=d, h = d+1;
				   Element hh = powh;
				   size_t log_hh = h;
				   hh = powh >> 1;
				   while (log_hh > d-1) {
				   --log_hh;
				   hh >>= 1;
				   }
				   if (!fastEEA(ai,bi,ci,di,m,d,n,hh, d-1)) return false;

*/		}

				Element ri = cur_ri;
				Element rinext = cur_rinext;
				Element qinext = cur_qinext;
				if (rinext==0) {
					cur_ainext = m+1;//infinity
					cur_qinext = cur_ainext;
					return true;
				}

				if (qinext <=0) {
					std::cout << "ERROR sth went very very wrong:" ;
					std::cout << "m:" << m << " n:" << n << " h:" << powh << "\n";
					std::cout << ai << " " << bi << "\n" << ci << " " << di <<"\n"; //getchar();
					return false;
				}

				Element ainext, binext, cinext,dinext;
				ainext=dinext=1;
				binext=cinext=0;
				ainext = cur_ainext;
				int K=-1;
				while (1) {
					++K;
					if (powh < ainext) {
						cur_ri = ri;
						cur_rinext = rinext;
						cur_ainext = ainext;
						cur_qinext = qinext;
						return true;
					}
					else {
						bi = ai;
						ai = ainext;
						Element temp = ci;
						//ci = ci*qinext + di;
						_intRing.axpy(ci, ci,qinext,di);
						++this->C.mul_counter;
						di = temp;

						temp = ri;
						ri=rinext;
						//rinext=qinext*ri+temp;
						//rinext = temp - qinext*ri;
						_intRing.axpy(rinext, -qinext,ri,temp);
						++this->C.mul_counter;
						if (rinext==0) {
							cur_ri = ri;
							cur_rinext = rinext;
							cur_ainext = m+1;
							cur_qinext = m+1;
							return true;
						}
						_intRing.quo(qinext,ri,rinext);
						++this->C.div_counter;
						_intRing.axpy(ainext, ai,qinext,bi);
						++this->C.mul_counter;
					}
				}

				return false;
		}
	};

	/* structures to perform MaxQ Fast Rational Reconstruction,
	 * stores not confirmed quotients and corresponding matrices
	 */
	template <class Ring>
	class QMatrix {
	public:
		Ring _intRing;
		typedef typename Ring::Element Element;

		Element a,b,c,d;
		Element q;

		QMatrix(const Ring Z) :
			_intRing(Z)
		{
			a=1;b=0;c=0;d=1;q=0;
		}

		QMatrix(const QMatrix& Q) :
			_intRing(Q._intRing)
		{
			a = Q.a;
			b = Q.b;
			c = Q.c;
			d = Q.d;
			q = Q.q;
		}

		QMatrix(Ring Z,const Element& ai, const Element& bi, const Element& ci, const Element& di, const Element& qi=0) :
			_intRing(Z)
		{
			a = ai;
			b = bi;
			c = ci;
			d = di;
			q = qi;
		}

		QMatrix& operator=(const QMatrix& Q) {
			a = Q.a;
			b = Q.b;
			c = Q.c;
			d = Q.d;
			q = Q.q;
			return *this;
		}

		void leftmultiply(const QMatrix Q) {
			leftmultiply(Q.a,Q.b, Q.c,Q.d);
		}

		void leftmultiply(const Element& ai,const Element& bi,const Element& ci,const Element& di) {
			Element tmpa,tmpb,tmpc,tmpd;
			_intRing.mul(tmpa,ai,a);
			_intRing.axpyin(tmpa,bi,c);

			_intRing.mul(tmpb,ai,b);
			_intRing.axpyin(tmpb,bi,d);

			_intRing.mul(tmpc,ci,a);
			_intRing.axpyin(tmpc,di,c);

			_intRing.mul(tmpd,ci,b);
			_intRing.axpyin(tmpd,di,d);

			a = tmpa; b = tmpb; c = tmpc; d = tmpd;
		}

		QMatrix& max(const QMatrix& max1, const QMatrix max2) {
			if (max1.q >= max2.q) return QMatrix(max1);
			else return QMatrix(max2);
		}

		bool maxin(const QMatrix max2) {
			if (q >= max2.q) return true;
			else {
				a = max2.a;
				b = max2.b;
				c = max2.c;
				d = max2.d;
				q = max2.q;
				return false;
			}
		}
	};

	template <class Ring>
	class myQueue: public std::deque<QMatrix<Ring > > {
	public:
		typedef typename Ring::Element Element;
		typedef QMatrix<Ring> QMatrix_;

		size_t _maxSize;
		size_t _size;

		myQueue(const size_t K=0) {
			_maxSize = K;
			_size = 0;
		}

		bool pushpop(QMatrix_& top, const QMatrix_& bottom) {
			if (_size+1 < _maxSize) {
				push_back(bottom);
				return false;
				++_size;
			}
			else {
				if (!this->empty()) {
					top = this->front();
					this->pop_front();
					push_back(bottom);

				}
				else {
					top=bottom;
				}
				return true;
			}
		}

		QMatrix_& clearmax(QMatrix_& max1) {
			while (!this->empty()) {
				QMatrix_ max2(this->front());
				if (max2.q > max1.q) return max1=max2;
				this->pop_front();
			}
		}
	};

	/*
	 *  implements fast rational reconstruction by Wan & Pan algorithm [Wan & Pan 2002]
	 *  MQRR Alg. of Monagan [Monagan2004] is used - maximal quotient is found
	 *  can be changed to large quotient for better performance q > m.bitsize() +c is returned
	 */
	template <class Ring>
	class FastMaxQRationalReconstruction: public FastRationalReconstruction<Ring> {
	public:
		const Ring _intRing;
		typedef typename Ring::Element Element;

		FastMaxQRationalReconstruction(const Ring& Z) :
			FastRationalReconstruction<Ring>(Z), _intRing(Z)
		{}

		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) const
		{
			bool res = fastQMaxReconstructRational(a,b,x,m);
			return res;
		}

		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) const
		{
			bool res= false;
			return res = FastRationalReconstruction<Ring>::reconstructRational(a,b,x,m,a_bound);
		}

	protected:
		mutable Element cur_ri;
		mutable Element cur_rinext;
		mutable Element cur_ainext;
		mutable Element cur_qinext;
		mutable Element T;
		mutable int c;


		bool fastQMaxReconstructRational(Element& n, Element& d, const Element& x, const Element& m) const
		{

			T = m.bitsize();
			c = 5; //should be changed here to enhance probability of correctness

			size_t log_m = m.bitsize()-1; //true unless m = 2^k

			Element ai, bi, ci, di;
			ai=1;bi=0;ci=0;di=1;

			cur_ri = m;
			cur_rinext = x;
			_intRing.quo(cur_ainext, cur_ri, cur_rinext);
			cur_qinext = cur_ainext;//  ri/ri_next

			Element powh;
			myQueue<Ring > queueMax(0);
			QMatrix<Ring > maxQ(_intRing);

			if (!fastQMaxEEA (ai,bi,ci,di,m,log_m,x,powtwo(powh,log_m+1), log_m+1,queueMax,maxQ)) {
				return false;
			}

#if 0
			if (cur_rinext != 0) {
				std::cout << "bad bounds - should not happen\n" << std::flush;
			}
			if (!queueMax.empty()) {
				std::cout << "Queue is not empty\n - sth wrong" << std::flush;
			}
#endif
			_intRing.mul(n,x, maxQ.a);
			_intRing.maxpyin(n,m,maxQ.c);
			//n = x_in*ai-m*ci;
			d = maxQ.a;

			//Element T = m.bitsize();int c = 5;
			if (maxQ.q.bitsize() > T.bitsize() + c) return true;
			else return false;
		}

		bool classicQMaxEEA(Element& ai, Element& bi, Element& ci, Element& di, const Element& r0, const Element& r1,const Element& powh, myQueue<Ring >&  queueMax, QMatrix<Ring>& maxQ) const
		{
			if (maxQ.q.bitsize() > T.bitsize() + c) return true;
			Element ri, rinext;
			ri = r0; rinext = r1;
			ai =di = 1;
			bi =ci = 0;
			if (rinext==0) return true;
			Element ainext,cinext;
			Element qinext;
			_intRing.quo(qinext,ri,rinext);
			++this->C.div_counter;

			ainext =qinext;
			while (ainext <= powh) {

				QMatrix<Ring> newQ(_intRing,ai,bi,ci,di,qinext);
				QMatrix<Ring> top(_intRing);
				if (queueMax.pushpop(top, newQ)) {
					if (maxQ.q < top.q) {
						maxQ = top;
						if (maxQ.q.bitsize() > T.bitsize() + c) return true;
					}
				}

				Element tmpbi = bi;
				Element tmpdi = di;
				bi = ai;
				ai = ainext;
				Element temp = ci;
				_intRing.axpy(ci, temp,qinext,di);
				++this->C.mul_counter;
				di = temp;

				temp = ri;
				ri=rinext;
				_intRing.axpy(rinext, -qinext,ri,temp);//rinext = ri-qinext*rinext
				++this->C.mul_counter;

				if (rinext==0) {
					//ainext = infinity
					cur_ri = ri;
					cur_rinext = rinext;
					cur_ainext = r0+1;
					cur_qinext = cur_ainext;//infinity
					return true;
				}

				_intRing.quo(qinext,ri,rinext);
				++this->C.div_counter;
				_intRing.axpy(ainext, ai,qinext,bi);
				++this->C.mul_counter;

			}
			cur_ri = ri;
			cur_rinext = rinext;
			cur_ainext = ainext;
			cur_qinext = qinext;
			return true;
		}


		bool fastQMaxEEA(Element& ai,Element& bi, Element& ci, Element& di,
				 const Element& m, const size_t d, const Element& n,
				 const Element& powh, const size_t& h, myQueue<Ring >&  queueMax, QMatrix<Ring>& maxQ) const
		{
			if (maxQ.q.bitsize() > T.bitsize() + c) return true;
			ai=Element(1);
			di=Element(1);
			bi=Element(0);
			ci=Element(0); //Q(0)=Id

			if (m==n) {
				cur_ri = n;
				cur_rinext = 0;
				cur_ainext = m+1;//infinity
				cur_qinext = m+1;
				ai=1; bi=1; ci=1; di=0;
				QMatrix<Ring> newQ(_intRing,1,0,0,1,1);
				QMatrix<Ring> top(_intRing);
				if (queueMax.pushpop(top, newQ)) {
					if (maxQ.q < top.q) {
						maxQ = top;
						if (maxQ.q.bitsize() > T.bitsize() + c) return true;
					}
				}
				return true;
			}

			if (powh < 1) return false; //should not happen
			if (n < 2) {
				//n==1 -> Q1=(m,1//1,0)
				//n==0 -> Q1=infinity
				cur_ri = m;
				cur_rinext = 1;
				cur_ainext = m;//infinity
				cur_qinext = m;
				if (cur_ainext > powh) {
					//we do not have to treat identity;
					return true;
				}
				else {
					ai=m; bi=1; ci=1; di=0;
					cur_ri = 1;
					cur_rinext = 0;
					cur_ainext = m+1;
					cur_qinext = m+1;
					QMatrix<Ring> newQ(_intRing,1,0,0,1,m);
					QMatrix<Ring> top(_intRing);
					if (queueMax.pushpop(top, newQ)) {
						if (maxQ.q < top.q) {
							maxQ = top;
							if (maxQ.q.bitsize() > T.bitsize() + c) return true;
						}
					}

					return true;
				}
				return true;
			}

			if (h < 1) {
				//if (qinext==1) { return (1,1,1,0) or (1,0,0,1) }
				if ((n << 1) > m) {
					cur_ri = n;
					cur_rinext = m-n;
					_intRing.quo(cur_qinext, cur_ri, cur_rinext);
					++this->C.div_counter;
					cur_ainext = cur_qinext + 1;
					ai=1; bi=1; ci=1; di=0;
					QMatrix<Ring> newQ(_intRing,1,0,0,1,1);
					QMatrix<Ring> top(_intRing);
					if (queueMax.pushpop(top, newQ)) {
						if (maxQ.q < top.q) {
							maxQ = top;
							if (maxQ.q.bitsize() > T.bitsize() + c) return true;
						}
					}
				}
				else {
					//we do not have to treat identity;
					cur_ri = m;
					cur_rinext = n;
					_intRing.quo(cur_ainext,m,n);
					++this->C.div_counter;
					cur_qinext = cur_ainext;
				}
				return true;
			}

			if (n.bitsize() < FastRationalReconstruction<Ring>::_threshold) {       //what about m?
				return classicQMaxEEA(ai,bi,ci,di,m,n,powh,queueMax,maxQ);
			}

			if (2*h+1 < d) {
				size_t lambda = d-2*h-1;

				Element aistar, bistar, cistar, distar;
				aistar=distar=1;
				bistar=cistar=0;

				Element mstar = m >> (long unsigned int) lambda;
				Element nstar = n >> (long unsigned int) lambda;

				size_t log_mstar = 2*h+1;

				queueMax._maxSize +=2;
				if (nstar > 0) if (!fastQMaxEEA(aistar, bistar, cistar, distar, mstar, log_mstar,nstar, powh, h, queueMax, maxQ)) return false;
				if (maxQ.q.bitsize() > T.bitsize() + c) return true;
				if (queueMax._size > 1) {
					queueMax.pop_back();
					QMatrix<Ring> Q_i_2 (queueMax.back());
					queueMax.pop_back();//pop Q*(i-1) q*(i)
					--queueMax._size ;
					--queueMax._size ;
					ai = Q_i_2.a;
					bi = Q_i_2.b;
					ci = Q_i_2.c;
					di = Q_i_2.d;// Q(i-2)= Q*(i-2) q*(i-1)
				}
				else {
					queueMax.clear();
				}
				queueMax._maxSize -=2;

				_intRing.mul(cur_ri,m,di);
				_intRing.maxpyin(cur_ri,n,bi);
				_intRing.mul(cur_rinext,n,ai);
				_intRing.maxpyin(cur_rinext,m,ci);
				this->C.mul_counter+=4;

				if (cur_ri < 0) {
					cur_ri = -cur_ri;
					cur_rinext = -cur_rinext;
				}
				if (cur_rinext>0) {
					_intRing.quo(cur_qinext,cur_ri,cur_rinext);
					++this->C.div_counter;
					_intRing.axpy(cur_ainext, ai,cur_qinext,bi);
					++this->C.mul_counter;
				}
				else {//should never happen
					cur_ainext = m+1;//infinity
					cur_qinext = cur_ainext;
				}
			}
			else { //if (h <= d-1)  //modificition of Wan&Pan
				Element a1,a2,b1,b2,c1,c2,d1,d2;
				a1=a2=d1=d2=1;
				b1=b2=c1=c2=0;

				Element sqrth;
				size_t logsqrth;
				logsqrth = h >> (int) 1;
				powtwo(sqrth, logsqrth);

				if (!fastQMaxEEA(a1,b1,c1,d1,m,d,n,sqrth, logsqrth, queueMax, maxQ)) return false;
				if (maxQ.q.bitsize() > T.bitsize() + c) return true;

				ai = a1; bi = b1; ci=c1; di = d1;

				Element ri = cur_ri;
				Element rinext = cur_rinext;

				size_t log_m;

				myQueue<Ring> queueTmp (queueMax._maxSize);
				QMatrix<Ring> maxQTmp(_intRing);
				if ((rinext > 0) && (cur_ainext <= powh)){

					QMatrix<Ring> newQ(_intRing,ai,bi,ci,di,cur_qinext);
					QMatrix<Ring> top(_intRing);
					if (queueMax.pushpop(top, newQ)) {
						if (maxQ.q < top.q) {
							maxQ = top;
							if (maxQ.q.bitsize() > T.bitsize() + c) return true;
						}
					}

					log_m = rinext.bitsize()-1;
					Element m2, n2;
					m2 = rinext;
					_intRing.axpy(n2, -cur_qinext, rinext, ri);
					++this->C.mul_counter;
					/* compute Q(i+1) */
					Element tmp = a1;
					a1 = cur_ainext;
					b1 = tmp;
					tmp = c1;
					_intRing.axpy (tmp, cur_qinext, c1, d1);
					++this->C.mul_counter;
					d1 = c1;
					c1 = tmp;

					size_t k = a1.bitsize()-1 ;
					int _k;
					if (h-k>2)
						_k = (int)(h-k-2);
					else _k = 0;

					if (n2 >0) {
						if (a1 < powh) {
							if (!fastQMaxEEA(a2,b2,c2,d2,m2,log_m,n2, powtwo(sqrth,_k), _k, queueTmp,maxQTmp)) return false;
						}
						else {
							ai = a1; bi = b1; ci=c1; di = d1;
							cur_ri = m2;
							cur_rinext = n2;
							_intRing.quo(cur_qinext,m2,n2);
							++this->C.div_counter;
							_intRing.axpy(cur_ainext,a1,cur_qinext,b1);
							++this->C.mul_counter;
							return true;
						}
					}
					else {
						ai = a1; bi = b1; ci=c1; di = d1;
						cur_ri = m2;
						cur_rinext = n2;
						cur_qinext = m+1;
						cur_ainext = m+1;
						return true;
					}
				}
				else {//ri_next == 0 || cur_ainext >powh
					ai = a1; bi = b1; ci=c1; di = d1;
					if (cur_rinext<=0) {
						cur_ainext = m +1;
						cur_qinext = cur_ainext;
					}
					//do not add matrix
					return true;
				}

				_intRing.mul(ai,b1,c2);
				_intRing.axpyin(ai,a1,a2);
				//aistar = a1*a2 + b1*c2;
				_intRing.mul(bi,b1,d2);
				_intRing.axpyin(bi, a1,b2);
				//bistar = a1*b2 + b1*d2;
				_intRing.mul(ci,d1,c2);
				_intRing.axpyin(ci, c1,a2);
				//cistar = c1*a2 + d1*c2;
				_intRing.mul(di,d1,d2);
				_intRing.axpyin(di, c1,b2);
				//distar = c1*b2 + d1*d2;
				this->C.mul_counter+=8;

				_intRing.mul(cur_ri,m,di);
				_intRing.maxpyin(cur_ri,n,bi);
				_intRing.mul(cur_rinext,n,ai);
				_intRing.maxpyin(cur_rinext,m,ci);
				this->C.mul_counter+=4;

				if (cur_ri < 0) {
					cur_ri = -cur_ri;
					cur_rinext = -cur_rinext;
				}

				if (cur_rinext>0) {
					_intRing.quo(cur_qinext,cur_ri,cur_rinext);
					++this->C.div_counter;
					_intRing.axpy(cur_ainext, ai,cur_qinext,bi);
					++this->C.mul_counter;
				}
				else {
					cur_ainext = m+1;//infinity
					cur_qinext = cur_ainext;
				}

				//multiply queueTmp by a1b1c1d1
				QMatrix<Ring > Q_i(_intRing,a1,b1,c1,d1);
				//update maximum
				if (maxQ.q < maxQTmp.q) {
					maxQTmp.leftmultiply(Q_i);
					maxQ = maxQTmp;
					if (maxQ.q.bitsize() > T.bitsize() + c) return true;
				}
				int K=0;
				QMatrix<Ring > Q(_intRing);
				while (!queueTmp.empty()) {
					if (ai > powh) {
						++K;
						Q=queueTmp.back();
						queueTmp.pop_back();
						--queueTmp._size;
						Q.leftmultiply(Q_i);
						cur_ainext = ai;
						cur_qinext = Q.q;
						Element tmp = cur_ri;
						_intRing.axpy(cur_ri,tmp,Q.q,rinext);
						cur_rinext = cur_ri;
						ai = Q.a;
						bi = Q.b;
						ci = Q.c;
						di = Q.d;
						//qi = Q.q;
					}
					else {
						//update queue;
						Q=queueTmp.front();
						queueTmp.pop_front();
						--queueTmp._size;
						if (maxQ.q < Q.q) {
							Q.leftmultiply(Q_i);
						}
						QMatrix<Ring > top(_intRing);
						if (queueMax.pushpop(top, Q)) {
							if (maxQ.q < top.q) {
								maxQ = top;
								if (maxQ.q.bitsize() > T.bitsize() + c) return true;
							}
						}
					}
				}
				if (K >0) {
					std::cout << "Error:" << K << " backward steps - should not happen\n"<< std::flush;
					return false;
				}
				return true;
#if 0
				//modification of Wan &Pan
				else {//h=d, h = d+1;
					Element hh = powh;
					size_t log_hh = h;
					hh = powh >> 1;
					while (log_hh > d-1) {
						--log_hh;
						hh >>= 1;
					}
					if (!fastQMaxEEA(ai,bi,ci,di,m,d,n,hh, d-1,queueMax,maxQ)) return false;
				}

#endif
			}

			Element ri = cur_ri;
			Element rinext = cur_rinext;
			Element qinext = cur_qinext;
			Element qi;
			if (rinext==0) {
				return true;
			}

			if (qinext <=0) {
				std::cout << "ERROR sth went very very wrong:"<< std::flush ;
				std::cout << "m:" << m << " n:" << n << " h:" << powh << "\n"<< std::flush;
				std::cout << ai << " " << bi << "\n" << ci << " " << di <<"\n"<< std::flush; //getchar();
				return false;
			}

			Element ainext, binext, cinext,dinext;
			ainext=dinext=1;
			binext=cinext=0;
			ainext = cur_ainext;
			int K=-1;
			while (1) {
				++K;
				if (powh < ainext) {
					cur_ri = ri;
					cur_rinext = rinext;
					cur_ainext = ainext;
					cur_qinext = qinext;

					return true;
				}
				else {

					QMatrix<Ring > Q(_intRing,ai,bi,ci,di,qinext);
					QMatrix<Ring > top(_intRing);
					if (queueMax.pushpop(top, Q)) {
						if (maxQ.q < top.q) {
							maxQ = top;
							if (maxQ.q.bitsize() > T.bitsize() + c) return true;
						}
					}

					bi = ai;
					ai = ainext;
					Element temp = ci;
					//ci = ci*qinext + di;
					_intRing.axpy(ci, ci,qinext,di);
					++this->C.mul_counter;
					di = temp;
					//qi = qinext;

					temp = ri;
					ri=rinext;
					//rinext=qinext*ri+temp;
					//rinext = temp - qinext*ri;
					_intRing.axpy(rinext, -qinext,ri,temp);
					++this->C.mul_counter;

					if (rinext==0) {
						cur_ri = ri;
						cur_rinext = rinext;
						cur_ainext = m+1;
						cur_qinext = m+1;
						return true;
					}
					_intRing.quo(qinext,ri,rinext);
					++this->C.div_counter;

					_intRing.axpy(ainext, ai,qinext,bi);
					++this->C.mul_counter;
				}
			}

			return false;
		}
	};

}

#endif //__LINBOX_fast_reconstruction_H
