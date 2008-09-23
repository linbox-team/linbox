/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#ifndef __LINBOXX__CLASSIC_RECONSTRUCTION_H__
#define __LINBOXX__CLASSIC_RECONSTRUCTION_H__

#include <iostream>

#include "linbox/algorithms/rational-reconstruction-base.h"

using namespace std;

namespace LinBox {
	
	template <class Ring>
	class ClassicRationalReconstruction: public RationalReconstructionBase<Ring> 
	{
	protected:
		bool _reduce;
		bool _recursive;
	public:
		//Ring _Z;
		typedef typename Ring::Element Element;
		
		ClassicRationalReconstruction(const Ring& Z, const bool reduce = true, const bool recursive = false): RationalReconstructionBase<Ring>(Z) {
			_reduce = reduce; _recursive = recursive;
		}
		
		ClassicRationalReconstruction<Ring> (const ClassicRationalReconstruction<Ring>& RR): RationalReconstructionBase<Ring>(RR._Z) {
			_reduce = RR._reduce; _recursive = RR._recursive;
		}
		
		~ClassicRationalReconstruction() {}
		
		//Wang method 
		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
			Element a_bound; _Z.sqrt(a_bound, m/2);
			bool res = reconstructRational(a,b,x,m,a_bound);
			res = res && (b <= a_bound);
			return res;
		}
		
		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {
			bool res=false;
			
			if (x == 0) {
				a = 0;
				b = 1;
			} else {
				res = ratrecon(a,b,x,m,a_bound);
				if (_recursive) {
					for(Element newbound = a_bound + 1; (!res) && (newbound<x) ; ++newbound)
						res = ratrecon(a,b,x,m,newbound);
				}
			}
			if (!res) {
				a = x> m/2? x-m: x; 
				b = 1;
				if (a > 0) res = (a < a_bound);
				else res = (-a < a_bound);
			}
			
			return res;		
		}
		
		
	protected:
		
		bool ratrecon(Element& a,Element& b,const Element& x,const Element& m,const Element& a_bound) {
			
			Element  r0, t0, q, u;
			r0=m;
			t0=0;
			a=x;
			b=1;
			//Element s0,s; s0=1,s=0;//test time gcdex;
			while(a>=a_bound)
			{
				
				q = r0;
				_Z.divin(q,a);        // r0/num
				++this->C.div_counter;		    
				
				u = a;
				a = r0;  	
				r0 = u;	// r0 <-- num
				
				_Z.axmyin(a,u,q); // num <-- r0-q*num
				++this->C.mul_counter;
				//if (a == 0) return false;
				
				u = b;
				b = t0;  	
				t0 = u;	// t0 <-- den
				
				_Z.axmyin(b,u,q); // den <-- t0-q*den
				++this->C.mul_counter;
				
			} 
			
			//                        if (den < 0) {
			//                                _Z.negin(num);
			//                                _Z.negin(den);
			//                        }
			
			if ((a>0) && (_reduce)) {
	    
				// [GG, MCA, 1999] Theorem 5.26
				// (ii)
				Element gg;
				++this->C.gcd_counter;
				if (_Z.gcd(gg,a,b) != 1) {
					
					Element ganum, gar2;
					for( q = 1, ganum = r0-a, gar2 = r0 ; (ganum >= a_bound) || (gar2<a_bound); ++q ) {
						ganum -= a;
						gar2 -= a;
					}
					
					//_Z.axmyin(r0,q,a);
					r0 = ganum;
					_Z.axmyin(t0,q,b);
					++this->C.mul_counter;++this->C.mul_counter;			
					if (t0 < 0) {
						a = -r0;
						b = -t0;
					} else {
						a = r0;
						b = t0;
					}
					
					//                                if (t0 > m/k) {
					if (abs((double)b) > (double)m/(double)a_bound) {
						if (!_recursive) {
							std::cerr 
								<< "*** Error *** No rational reconstruction of " 
								<< x 
								<< " modulo " 
								<< m 
								<< " with denominator <= " 
								<< (m/a_bound)
								<< std::endl; 
						}
						return false;
					}
					if (_Z.gcd(gg,a,b) != 1) {
						if (!_recursive) 
							std::cerr 
								<< "*** Error *** There exists no rational reconstruction of " 
								<< x 
								<< " modulo " 
								<< m 
								<< " with |numerator| < " 
								<< a_bound
								<< std::endl
								<< "*** Error *** But " 
								<< a
								<< " = "
								<< b
								<< " * "
								<< x
								<< " modulo " 
								<< m 
								<< std::endl;
						return false;
					}
				}
			}
			// (i)
			if (b < 0) {
				_Z.negin(a);
				_Z.negin(b);
			}
			
			// std::cerr << "RatRecon End " << num << "/" << den << std::endl;
			return true;    
		}
	};
	
	template <class Ring>
	class ClassicMaxQRationalReconstruction:public ClassicRationalReconstruction<Ring> {
	public:
		typedef typename Ring::Element Element;
		
		ClassicMaxQRationalReconstruction(const Ring& Z, const bool reduce = true, const bool recursive = false):ClassicRationalReconstruction<Ring>(Z,reduce,recursive) {}
		
		ClassicMaxQRationalReconstruction(const ClassicMaxQRationalReconstruction<Ring>& RR): ClassicRationalReconstruction<Ring>(RR) {}
		
		~ClassicMaxQRationalReconstruction() {}
		
		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
			bool res = maxEEA(a,b,x,m);
			return res;    
		}
		
		bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {
			bool res= false;
			return res = ClassicRationalReconstruction<Ring>::reconstructRational(a,b,x,m,a_bound);    
		}  
		
	protected: 
		bool maxEEA(Element& a, Element& b, const Element& x, const Element& m) {
			
			Element qmax = 0, amax=x, bmax =1; 
			
			Element  r0, t0, q, u;
			r0=m;
			t0=0;
			a=x;
			b=1;
			//Element s0,s; s0=1,s=0;//test time gcdex;
			while(a>0)
			{
				q = r0;
				_Z.divin(q,a);        // r0/num
				++this->C.div_counter;
				//cout << r0<< " " << a << " " <<  q << ",";
				if (q > qmax) {
					//cout << r0<< " " << a << " " <<  q << ",";
					amax = a;
					bmax = b;
					qmax = q;
				}
				
				u = a;
				a = r0;  	
				r0 = u;	// r0 <-- num
				
				_Z.axmyin(a,u,q); // num <-- r0-q*num
				++this->C.mul_counter;
	//if (a == 0) return false;
				
				u = b;
				b = t0;  	
				t0 = u;	// t0 <-- den
				
				_Z.axmyin(b,u,q); // den <-- t0-q*den
				++this->C.mul_counter;
				
			} 
			
			a = amax;
			b = bmax;
			
			if (b < 0) {
				_Z.negin(a);
				_Z.negin(b);
			}
			
			Element gg;
			_Z.gcd(gg,a,b);
			++this->C.gcd_counter;
			if (gg > 1) return false; 
			else return true;
		}
	};
	
}
#endif

