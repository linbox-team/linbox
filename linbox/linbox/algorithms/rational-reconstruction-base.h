/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#include <iostream>
#include <queue>

#ifndef __LINBOXX__RECONSTRUCTION_BASE_H__
#define __LINBOXX__RECONSTRUCTION_BASE_H__

//#define __FASTRR_DEFAULT_THRESHOLD 9223372036854775807 //2^63-1
//#define __FASTRR_DEFAULT_THRESHOLD 2147483647 //2^31-1 
#define __FASTRR_DEFAULT_THRESHOLD 7 
//max int=2^32 4294967295
//max int=2^64 18446744073709551616
//__FASTRR_DEFAULT_THRESHOLD > 4

using namespace std;

namespace LinBox {

template <class Ring, class RRBase>
class RationalReconstruction
{
public:
	Ring _Z;
	RRBase _RR;

	typedef typename Ring::Element Element;
  
	RationalReconstruction(const Ring& Z):_Z(Z), _RR(Z) {}
	RationalReconstruction(const RRBase& RR): _Z(RR._Z), _RR(RR) {}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
		
		Element x_in(x);
		if (x<0) {
			if ((-x)>m)
				x_in %= m;
			if (x<0)
				x_in += m;
		} else {
			if (x>m)
				x_in %= m;
		}
		bool res= _RR.reconstructRational(a,b,x_in,m);
		_RR.write(cout);
		return res;
	}
  
	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {
		
		Element x_in(x);
		if (x<0) {
			if ((-x)>m)
				x_in %= m;
			if (x<0)
				x_in += m;
		} else {
			if (x>m)
				x_in %= m;
		}
		bool res= _RR.reconstructRational(a,b,x_in,m,a_bound);
		_RR.write(cout);
		return res;
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound, const Element& b_bound) {
		Element x_in(x);
		if (x<0) {
			if ((-x)>m)
				x_in %= m;
			if (x<0)
				x_in += m;
		} else {
			if (x>m)
				x_in %= m;
		}
		Element bound = x/b_bound;
		_RR.reconstructRational(a,b,x,m,(bound>a_bound?bound:a_bound));
		bool res=  (b > b_bound)? false: true;
		_RR.write(cout);
		return res;
	}

  //fastReconstruction();

  //classicReconstruction();

};

class OpCounter {
public:
	size_t div_counter;
	size_t mul_counter;
	size_t gcd_counter;

	OpCounter() {
		div_counter=0; mul_counter=0; gcd_counter=0;
	}

	void write(ostream& is) {
		is << div_counter << " divisions\n";
		is << mul_counter << " multiplications\n";
		is << gcd_counter << " gcds\n";
	}
};

template <class Ring>
class RationalReconstructionBase
{
public:
	Ring _Z;
	typedef typename Ring::Element Element;
	
	OpCounter C;

	RationalReconstructionBase(const Ring& Z): _Z(Z) {}
	RationalReconstructionBase(const RationalReconstructionBase<Ring>& RR): _Z(RR._Z) {}

	virtual bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m)=0;

	virtual bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound)=0;

	virtual ~RationalReconstructionBase() {}

	void write(ostream& is) {
		C.write(is);
	}
};
/*
template <class Ring>
class WangRationalReconstruction: public RationalReconstructionBase<Ring> 
{
protected:
	bool _reduce;
	bool _recursive;
public:
	Ring _Z;
	typedef typename Ring::Element Element;
	
	OpCounter C;

	WangRationalReconstruction(const Ring& Z, const bool reduce = true, const bool recursive = false): RationalReconstructionBase<Ring>(Z) {
		_reduce = reduce; _recursive = recursive;
	}

	WangRationalReconstruction<Ring> (const WangClassicRationalReconstruction<Ring>& RR): RationalReconstructionBase<Ring>(RR._Z) {
		_reduce = RR._reduce; _recursive = RR._recursive;
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m)=0 {
		Element a_bound; _Z.sqrt(a_bound, m/2);
		bool res = rationalReconstruction(a,b,x,m,a_bound);
		res = res && (b <= a_bound);
		return res;
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {
		bool res;

		if (x == 0) {
		    a = 0;
		    b = 1;
		} else {
			bool res = ratrecon(a,b,x,m,a_bound);
			if (_recursive)
			for(Element newbound = a_bound + 1; (!res) && (newbound<x) ; ++newbound)
			    res = ratrecon(a,b,x,m,newbound);
		}
		if (!res) {
			a = x> m/2? x-m: x; 
			b = 1;
			if (a > 0) res = (a < a_bound);
			else res = (-a < a_bound);
		}
		
		return res;		
	}

	virtual ~WangRationalReconstruction() {}
protected:

// m> x> 0 
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
++C.div_counter;		    

		    u = a;
		    a = r0;  	
		    r0 = u;	// r0 <-- num
                                                           
		    _Z.axmyin(a,u,q); // num <-- r0-q*num
++C.mul_counter;
		    //if (a == 0) return false;
		    
		    u = b;
		    b = t0;  	
		    t0 = u;	// t0 <-- den
		               
		    _Z.axmyin(b,u,q); // den <-- t0-q*den
++C.mul_counter;
		} 

//                        if (den < 0) {
//                                _Z.negin(num);
//                                _Z.negin(den);
//                        }

		if ((a>0) && (_reduce)) {

			// [GG, MCA, 1999] Theorem 5.26
			// (ii)
		    Element gg;
++C.gcd_counter;
		    if (_Z.gcd(gg,a,b) != 1) {
			
			Element ganum, gar2;
			for( q = 1, ganum = r0-a, gar2 = r0 ; (ganum >= a_bound) || (gar2<a_bound); ++q ) {
			    ganum -= a;
			    gar2 -= a;
			}
					
			//_Z.axmyin(r0,q,a);
			r0 = ganum;
			_Z.axmyin(t0,q,b);
++C.mul_counter;++C.mul_counter;			
			if (t0 < 0) {
			    a = -r0;
			    b = -t0;
			} else {
			    a = r0;
			    b = t0;
			}
			
//                                if (t0 > m/k) {
			if ((double)b > (double)m/(double)a_bound) {
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
*/
template <class Ring>
class QMatrix {
public:
	Ring _Z;
	typedef typename Ring::Element Element;

	Element a,b,c,d;
	Element q;

	QMatrix(const Ring Z): _Z(Z) {
		a=1;b=0;c=0;d=1;q=0;
	}

	QMatrix(const QMatrix& Q): _Z(Q._Z) {
		a = Q.a;
		b = Q.b;
		c = Q.c;
		d = Q.d;
		q = Q.q;
	}

	QMatrix(Ring Z,const Element& ai, const Element& bi, const Element& ci, const Element& di, const Element& qi=0): _Z(Z) {
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
		_Z.mul(tmpa,ai,a);
		_Z.axpyin(tmpa,bi,c);

		_Z.mul(tmpb,ai,b);
		_Z.axpyin(tmpb,bi,d);

		_Z.mul(tmpc,ci,a);
		_Z.axpyin(tmpc,di,c);
		
		_Z.mul(tmpd,ci,b);
		_Z.axpyin(tmpd,di,d);		

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
class myQueue: public deque<QMatrix<Ring > > {
public:
	typedef typename Ring::Element Element;
	typedef QMatrix<Ring> QMatrix;

	size_t _maxSize;
	size_t _size;

	myQueue(const size_t K=0) {
		_maxSize = K;
		_size = 0;
	}

	bool pushpop(QMatrix& top, const QMatrix& bottom) {
		if (_size+1 < _maxSize) {
			push_back(bottom);
			return false;
			++_size;
		} else {
			if (!empty()) {
				top = front();
				pop_front();
				push_back(bottom);
				
			} else {
				top=bottom; 
			}
			return true;					
		}
	}
	
	QMatrix& clearmax(QMatrix& max1) {
		while (!empty()) {
			QMatrix max2(front());
			if (max2.q > max1.q) return max1=max2;
			pop_front();
		}
	}
};

template <class Ring>
class MaxQFastRationalReconstruction: public RationalReconstructionBase<Ring>
{
protected:
	size_t _threshold;

typedef typename Ring::Element Element;
public:

	MaxQFastRationalReconstruction(const Ring& Z): RationalReconstructionBase<Ring>(Z) {		
		_threshold = __FASTRR_DEFAULT_THRESHOLD;
		if (_threshold <5) _threshold=5; 
	}

	~MaxQFastRationalReconstruction() {}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {  
		//reconstructRational(a,b,x,m,1);
		return fastQMaxReconstructRational(a,b,x,m);
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {
		bool res = reconstructRational(a,b,x,m);
		return res && (a < a_bound);
	}

protected:

	Element cur_ri;
	Element cur_rinext;
	Element cur_ainext;
	Element cur_qinext;

	//QMatrix max;

	Element& powtwo(Element& h, const Element& log_h) {
		h = 1;
		//cout << "max" << ULONG_MAX << "x" << x << "?" << (x < ULONG_MAX);
		if (log_h < ULONG_MAX) {
			h<<=log_h;
			//cout << "z"<< z;
			return h;
		} else {
			Element n,m;
			quoRem(n,m,log_h,(Element)(ULONG_MAX-1));
			for (int i=0; i < n; ++i) {
				h <<=(long int)(ULONG_MAX-1);
			}
			h <= (long int)m;
			return h;
		}
	}

	Element& powtwo(Element& h, const size_t log_h) {
		h = 1;
		h<<=log_h;
		return h;
	}

	bool fastQMaxReconstructRational(Element& n, Element& d, const Element& x, const Element& m) {
	       	
		size_t log_m = m.bitsize()-1; //true unless m = 2^k
		
		Element ai, bi, ci, di;
		ai=1;bi=0;ci=0;di=1;
		
		cur_ri = m;
		cur_rinext = x;
		_Z.quo(cur_ainext, cur_ri, cur_rinext);
		cur_qinext = cur_ainext;//  ri/ri_next

		Element powh;
		myQueue<Ring > queueMax(0);
		QMatrix<Ring > maxQ(_Z);
		
		if (!fastQMaxEEA (ai,bi,ci,di,m,log_m,x,powtwo(powh,log_m+1), log_m+1,queueMax,maxQ)) {
			cout << "false\n" << flush;
			return false;
		}
		
		if (cur_rinext != 0) {
			cout << "bad bounds - should not happen\n" << flush;
		}

		if (!queueMax.empty()) {
			cout << "Queue is empty\n - sth wrong" << flush;
		}
		
		_Z.mul(n,x, maxQ.a);
		_Z.axmyin(n,m,maxQ.c);
		//n = x_in*ai-m*ci;
		d = maxQ.a;

/*
		_Z.mul(n, maxQ.d, m);
		_Z.axmyin(n, maxQ.b, n);
		d = maxQ.b;	
*/
		Element _gcd;
		if (_Z.gcd(_gcd, n,d) > 1) {
		/* solution not implemented */
			cout << "Solution pair is not coprime\n";
			return false;
		}

		return true;
	}		

	bool classicQMaxEEA(Element& ai, Element& bi, Element& ci, Element& di, const Element& r0, const Element& r1,const Element& powh, myQueue<Ring >&  queueMax, QMatrix<Ring>& maxQ) {
		Element ri, rinext;
		ri = r0; rinext = r1;
		ai =di = 1;
		bi =ci = 0;
		if (rinext==0) return true;
		Element ainext,cinext;
		Element qinext;
		_Z.quo(qinext,ri,rinext);
++C.div_counter;

		ainext =qinext;	
		//_Z.axpy(ainext, ai,qinext,bi);
		while (ainext <= powh) {

			//_Z.quo(qinext,ri,rinext);
++C.div_counter;
                        QMatrix<Ring> newQ(_Z,ai,bi,ci,di,qinext);
                        QMatrix<Ring> top(_Z);
                        if (queueMax.pushpop(top, newQ)) {
				cout << "1new max " << top.q << "\n" << flush;
                                if (maxQ.q < top.q) maxQ = top;
                        }
			cout << "EEA" << qinext << ",";
			//++i;
			//ainext = ai*qinext + bi;
			//_Z.axpy(ainext, ai,qinext,bi);
++C.mul_counter;

                        Element tmpbi = bi;
			Element tmpdi = di; 
			bi = ai;
			ai = ainext;
			Element temp = ci;
			_Z.axpy(ci, temp,qinext,di);
++C.mul_counter;
			di = temp;

			temp = ri;
			ri=rinext;
			_Z.axpy(rinext, -qinext,ri,temp);//rinext = ri-qinext*rinext

++C.mul_counter;
			if (rinext==0) {
				//ainext = infinity
			        cur_ri = ri;
			        cur_rinext = rinext;
			        cur_ainext = r0+1;
				cur_qinext = cur_ainext;//infinity
				cout << "->1E:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "\n"<< flush;
				return true;
			}

			_Z.quo(qinext,ri,rinext);
++C.div_counter;			
			_Z.axpy(ainext, ai,qinext,bi);
			
		}
		cur_ri = ri;
		cur_rinext = rinext;
		cur_ainext = ainext;
		cur_qinext = qinext;	
		cout << "->2E:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "\n"<< flush;
		return true;
	}
	
	bool fastQMaxEEA(Element& ai,Element& bi, Element& ci, Element& di, const Element& m, const size_t d, const Element& n, const Element& powh, const size_t& h, myQueue<Ring >&  queueMax, QMatrix<Ring>& maxQ) {
		
		cout << m << " " << n << " " << powh << "\n" << flush;
		
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
			QMatrix<Ring> newQ(_Z,1,0,0,1,1);
			//QMatrix<Ring> newQ(_Z,1,1,1,0,1);
			QMatrix<Ring> top(_Z);
			if (queueMax.pushpop(top, newQ)) {
				cout << "2new max " << top.q << "\n"<< flush;
				if (maxQ.q < top.q) maxQ = top;
			}
			cout << "m=n1" << ",";
			return true;
		}

		if (powh < 1) return false; //should not happen
		if (n < 2) {
			cout << "n=1\n"<< flush;
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
				QMatrix<Ring> newQ(_Z,1,0,0,1,m);
				//QMatrix<Ring> newQ(_Z,m,1,1,0,m);
				QMatrix<Ring> top(_Z);
				if (queueMax.pushpop(top, newQ)) {
				//	cout << "3new max " << top.q << "\n"<< flush;
					if (maxQ.q < top.q) maxQ = top;
				}
				cout << "n=1" << m << ",";
				
				return true;
			}
		        return true;
		}

		if (h < 1) {
			cout << "h < 1" << flush;
			//quo(qinext,m,n); 
			//if (qinext==1) { return (1,1,1,0) or (1,0,0,1)
			if ((n << 1) > m) {
				cur_ri = n;
				cur_rinext = m-n;
				_Z.quo(cur_qinext, cur_ri, cur_rinext);
++C.div_counter;
				cur_ainext = cur_qinext + 1;
				ai=1; bi=1; ci=1; di=0;
				QMatrix<Ring> newQ(_Z,1,0,0,1,1);
				//QMatrix<Ring> newQ(_Z,1,1,1,0,1);
				QMatrix<Ring> top(_Z);
				if (queueMax.pushpop(top, newQ)) {
					cout << "4new max " << top.q << "\n"<< flush;
					if (maxQ.q < top.q) maxQ = top;
				}
				cout << "h=01" << ",";

			} else {
				//we do not have to treat identity;
				cur_ri = m;
				cur_rinext = n;
				_Z.quo(cur_ainext,m,n);
++C.div_counter;
				cur_qinext = cur_ainext;
			}
			return true;
		}

		if (n < _threshold) {       //what about m?
			return classicQMaxEEA(ai,bi,ci,di,m,n,powh,queueMax,maxQ);
		}

		size_t log_n = n.bitsize()-1;

		//		cout << d << " " << log_n << " "  << h << "\n"<< flush;

		if (2*h+1 < d) {
			cout << "choice1"<< flush;
			//if (2*h+1 <= log_n) {
		  if (true) {	
			cout << ".1\n"<< flush;
			size_t lambda = d-2*h-1;

			Element aistar, bistar, cistar, distar;
			aistar=distar=1;
			bistar=cistar=0;

			Element mstar = m >> (long unsigned int) lambda;
			Element nstar = n >> (long unsigned int) lambda;

			size_t log_mstar = 2*h+1;

			queueMax._maxSize +=2;
			if (nstar > 0) if (!fastQMaxEEA(aistar, bistar, cistar, distar, mstar, log_mstar,nstar, powh, h, queueMax, maxQ)) return false;

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
			}  else {
				queueMax.clear();
			}
		       	queueMax._maxSize -=2;
			
			_Z.mul(cur_ri,m,di);
			_Z.axmyin(cur_ri,n,bi);
			_Z.mul(cur_rinext,n,ai);
			_Z.axmyin(cur_rinext,m,ci);
C.mul_counter+=4;
			//ri = m*di - n*bi;
			//rinext = -m*ci + n*ai; 
 
			if (cur_ri < 0) {
				cur_ri = -cur_ri;
				cur_rinext = -cur_rinext;
			}
			if (cur_rinext>0) {
				_Z.quo(cur_qinext,cur_ri,cur_rinext);
++C.div_counter;
				_Z.axpy(cur_ainext, ai,cur_qinext,bi); 
++C.mul_counter;
			} else {//should never happen
				cur_ainext = m+1;//infinity
				cur_qinext = cur_ainext;
			}

			cout << "\n2 backward steps\n"<< flush;
			cout << "->1:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "->\n"<< flush;

		   }
		} else if (h <= d-1) {
			cout << "choice2\n"<< flush;
			Element a1,a2,b1,b2,c1,c2,d1,d2;
                        a1=a2=d1=d2=1;
                        b1=b2=c1=c2=0;

                        Element sqrth;
			size_t logsqrth;
                                //sqrt(sqrth,powh);
                        logsqrth = h >> (int) 1;
			powtwo(sqrth, logsqrth);
                                //logtwo(logsqrth, sqrth);
                                //_sqrth = h/2;

                                //if (!EEA(sgn, ri, rinext, qinext, a1,b1,c1,d1,m,d,n,sqrth, logsqrth)) return 0;
			if (!fastQMaxEEA(a1,b1,c1,d1,m,d,n,sqrth, logsqrth, queueMax, maxQ)) return 0;

			ai = a1; bi = b1; ci=c1; di = d1;

			Element ri = cur_ri;
			Element rinext = cur_rinext;

			size_t log_m;
			//if (ri !=m) //logtwo(log_m,ri);

			myQueue<Ring> queueTmp (queueMax._maxSize+2);
			QMatrix<Ring> maxQTmp(_Z);
			if ((rinext > 0) && (cur_ainext <= powh)){

                                QMatrix<Ring> newQ(_Z,ai,bi,ci,di,cur_qinext);
                                QMatrix<Ring> top(_Z);
                                if (queueMax.pushpop(top, newQ)) {
                                        cout << "6new max " << top.q << "\n";
                                        if (maxQ.q < top.q) maxQ = top;
                                }
				cout << "2step1:" << cur_qinext << ",";


				log_m = rinext.bitsize()-1;
				Element m2, n2;
				m2 = rinext;
				_Z.axpy(n2, -cur_qinext, rinext, ri);
++C.mul_counter;
				/* compute Q(i+1) */
				Element tmp = a1;
				a1 = cur_ainext;
				b1 = tmp;
				tmp = c1;
				_Z.axpy (tmp, cur_qinext, c1, d1);
++C.mul_counter;
				d1 = c1;
				c1 = tmp;

				//add matrix
				//_Z.quo(cur_qinext, m2,n2);
				//_Z.axpy(cur_ainext, cur_qinext, ai, bi); not needed
				
				ai = a1; bi = b1; ci=c1; di = d1;

				size_t k = a1.bitsize()-1 ;
				if (n2 >0) {
					if (a1 < powh) {
						if (!fastQMaxEEA(a2,b2,c2,d2,m2,log_m,n2, powtwo(sqrth,h-k-1), h-k-1, queueTmp,maxQTmp)) return 0;
					} else {
						cur_ri = m2;
						cur_rinext = n2;
						_Z.quo(cur_qinext,m2,n2);
++C.div_counter;
						_Z.axpy(cur_ainext,a1,cur_qinext,b1);
++C.mul_counter;
						return 1;
					}
				} else {
					cur_ri = m2;
					cur_rinext = n2;
					cur_qinext = m+1;
					cur_ainext = m+1;
					return 1;
				}
			} else {
				//do not add matrix
				ai = a1; bi = b1; ci=c1; di = d1;
				cout << "End2\n"<< flush;
				return 1;
			}

			_Z.mul(ai,b1,c2);
			_Z.axpyin(ai,a1,a2);
                               //aistar = a1*a2 + b1*c2;
			_Z.mul(bi,b1,d2);
			_Z.axpyin(bi, a1,b2);
                               //bistar = a1*b2 + b1*d2;
			_Z.mul(ci,d1,c2);
			_Z.axpyin(ci, c1,a2);
                               //cistar = c1*a2 + d1*c2;
			_Z.mul(di,d1,d2);
			_Z.axpyin(di, c1,b2);
                               //distar = c1*b2 + d1*d2;
C.mul_counter+=8;

			_Z.mul(cur_ri,m,di);
			_Z.axmyin(cur_ri,n,bi);
			_Z.mul(cur_rinext,n,ai);
			_Z.axmyin(cur_rinext,m,ci);
C.mul_counter+=4;
			//ri = m*di - n*bi;
			//rinext = -m*ci + n*ai;
			if (cur_ri < 0) {
				cur_ri = -cur_ri;
				cur_rinext = -cur_rinext;
			}
			cout << "->2:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "->\n"<< flush;
			if (cur_rinext>0) {
				_Z.quo(cur_qinext,cur_ri,cur_rinext);
++C.div_counter;
				_Z.axpy(cur_ainext, ai,cur_qinext,bi); 
++C.mul_counter;
			} else {
				cur_ainext = m+1;//infinity
				cur_qinext = cur_ainext;
			}

//multiply queueTmp by a1b1c1d1
                        QMatrix<Ring > Q_i(_Z,a1,b1,c1,d1);
			//update maximum
			if (maxQ.q < maxQTmp.q) {
				maxQTmp.leftmultiply(Q_i);
				maxQ = maxQTmp;
			}
			int K=0;
			//if (!queueTmp.empty()) { 
			//Element qi;
			QMatrix<Ring > Q(_Z);
				//Q = queueTmp.back();
				//queueTmp.pop_back();
				//Q.leftmultiply(Q_i);
				//qi = Q.q;
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
					_Z.axpy(cur_ri,tmp,Q.q,rinext);
					cur_rinext = cur_ri; 
					ai = Q.a;
					bi = Q.b;
					ci = Q.c;
					di = Q.d;
					//qi = Q.q;
				} else {
					//update queue;
					Q=queueTmp.front();
					queueTmp.pop_front();
					--queueTmp._size;
					if (maxQ.q < Q.q) {
						Q.leftmultiply(Q_i);
					}
					QMatrix<Ring > top(_Z);
					if (queueMax.pushpop(top, Q)) {
						cout << "8new max " << top.q << "\n"<< flush;
						if (maxQ.q < top.q) maxQ = top;
					}
					cout << "queue" << Q.q << ",";

					//else: queueMax may contain invalid matrices, but if Q.q is max then Q is correct
				}
			}
			cout << "\n" << K << " backward steps\n"<< flush;
			cout << "->End:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n"<< flush;
                        return 1;
			//}                               

		} else {//h=d, h = d+1;
			cout << "choice3\n"<< flush;
			Element hh = powh;
			size_t log_hh = h;
			hh = powh >> 1;
			while (log_hh > d-1) {
				--log_hh;
				hh >>= 1;
			}
			if (!fastQMaxEEA(ai,bi,ci,di,m,d,n,hh, d-1,queueMax,maxQ)) return 0;

				cout << "->3:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "->\n"<< flush;
		}
		
		Element ri = cur_ri;
		Element rinext = cur_rinext;
		Element qinext = cur_qinext;
		Element qi;
		//_Z.mul(ri,m,di);
		//_Z.axmyin(ri,n,bi);
		//_Z.mul(rinext,n,ai);
		//_Z.axmyin(rinext,m,ci);
                ////ri = m*di - n*bi;
                ////rinext = -m*ci + n*ai;
                if (rinext==0) {
			cout << "\n0 forward steps\n";
			cout << "->End1:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
			return 1;
		}
                //if (ri < 0) {
		//	ri = -ri;
                //        rinext = -rinext;
                //}
                //quo(qinext,ri,rinext);
				
		if (qinext <=0) {
			cout << "ERROR sth went very very wrong:"<< flush ;
                        cout << "m:" << m << " n:" << n << " h:" << powh << "\n"<< flush;
                        cout << ai << " " << bi << "\n" << ci << " " << di <<"\n"<< flush; //getchar();
                        return 0;
                }

                Element ainext, binext, cinext,dinext;
                ainext=dinext=1;
                binext=cinext=0;
		ainext = cur_ainext;
		int K=-1;
                while (1) {
			++K;
			//_Z.axpy(ainext, ai,qinext,bi);
			//ainext = cur_ainext;
			if (powh < ainext) {
				cur_ri = ri;
				cur_rinext = rinext;
				cur_ainext = ainext;
				cur_qinext = qinext;
				cout << "\n" << K << " forward steps\n"<< flush;
				cout << "->End2:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n"<< flush;

				return 1;
			}
			//else if ((ah==1) && ((int)ah << (int)h) != ainext) return 1;
			else {

                                QMatrix<Ring > Q(_Z,ai,bi,ci,di,qinext);
                                QMatrix<Ring > top(_Z);
                                if (queueMax.pushpop(top, Q)) {
                                        cout << "9new max " << top.q << "\n"<< flush;
                                        if (maxQ.q < top.q) maxQ = top;
                                }
				cout << "End" << qinext << ",";


				bi = ai;
				ai = ainext;
				Element temp = ci;
				//ci = ci*qinext + di;
				_Z.axpy(ci, ci,qinext,di);
++C.mul_counter;
				di = temp;
				//qi = qinext;

				temp = ri;
				ri=rinext;
				//rinext=qinext*ri+temp;
				//rinext = temp - qinext*ri;
				_Z.axpy(rinext, -qinext,ri,temp);
++C.mul_counter;

				if (rinext==0) {
					cur_ri = ri;
					cur_rinext = rinext;
					cur_ainext = m+1;
					cur_qinext = m+1;
					cout << "\n" << K << " forward steps\n"<< flush;
					cout << "->End3:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n"<< flush;
					return 1;
				}
                                        //++counter;
				_Z.quo(qinext,ri,rinext);
++C.div_counter;


                                _Z.axpy(ainext, ai,qinext,bi);
++C.mul_counter; 
			}
		}

                return 0;
	}

};

template <class Ring>
class MaxQClassicRationalReconstruction: public RationalReconstructionBase<Ring>
{
typedef typename Ring::Element Element;
public:

	MaxQClassicRationalReconstruction(const Ring& Z):RationalReconstructionBase<Ring>(Z) {}

	~MaxQClassicRationalReconstruction() {}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
		bool res = maxExGcd(a,b,x,m);
		return res;
			
	}	

	/* for sake of completeness, should not be used */
	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {

		bool res =maxExGcd(a,b,x,m);
		return res && (a < a_bound);				
	}

protected:
	/* returns a,b corresponding to the maximal quotient of Euclide GCD algorithm 
	 * returns true if a.b are coprime false otherwise
	 */
	bool maxExGcd(Element& a, Element& b, const Element& x, const Element& m) {
		
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
++C.div_counter;
			cout << r0<< " " << a << " " <<  q << ",";
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
++C.mul_counter;
		    //if (a == 0) return false;
		    
		    u = b;
		    b = t0;  	
		    t0 = u;	// t0 <-- den
		               
		    _Z.axmyin(b,u,q); // den <-- t0-q*den
++C.mul_counter;

		} 

		a = amax;
		b = bmax;
		
                if (b < 0) {
			_Z.negin(a);
                        _Z.negin(b);
                }

		Element gg;
		_Z.gcd(gg,a,b);
++C.gcd_counter;
		if (gg > 1) return false; 
		else return true;
	}
};


/* specialization for GMP integers - bitsize(), << , >>*/  
template <class Ring>
class WangFastRationalReconstruction: public RationalReconstructionBase<Ring>
{

typedef typename Ring::Element Element;
protected:
	size_t _threshold;

public:
	WangFastRationalReconstruction(const Ring& Z): RationalReconstructionBase<Ring>(Z) {
		_threshold = __FASTRR_DEFAULT_THRESHOLD;
		if (_threshold < 5) _threshold = 5;
	}

	~WangFastRationalReconstruction() {}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
		Element a_bound; _Z.sqrt(a_bound, m/2);
		//if (b_bound >0)
		reconstructRational(a,b,x,m,a_bound);
		//else return false;
		return (a < a_bound);
		//b_bound == a_bound
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound ) {
		if (x < a_bound) {
			//case(i)
			a=x;
			b=1;
			return true;
		} else if (m-x < a_bound) {
			a = x-m;
			b = 1;
			return true;
		}
		
		Element bound = a_bound << 1;

		if (m/bound > 1) {
			fastReconstructRational(a,b,x,m,m/bound);
		
			if (_Z.abs(a) < a_bound) {
				return true;
			}
			return false;
		} else {
			//either case (i) or false
			return false;
		}
	}

	/* m/2 >= d_bound > 1 
	 * m/2 >= bound >=2
	 */ 
protected:

	Element cur_ri;
	Element cur_rinext;
	Element cur_ainext;
	Element cur_qinext;

	Element& powtwo(Element& h, const Element& log_h) {
		h = 1;
		//cout << "max" << ULONG_MAX << "x" << x << "?" << (x < ULONG_MAX);
		if (log_h < ULONG_MAX) {
			h<<=log_h;
			//cout << "z"<< z;
			return h;
		} else {
			Element n,m;
			quoRem(n,m,log_h,(Element)(ULONG_MAX-1));
			for (int i=0; i < n; ++i) {
				h <<=(long int)(ULONG_MAX-1);
			}
			h <= (long int)m;
			return h;
		}
	}

	Element& powtwo(Element& h, const size_t log_h) {
		h = 1;
		h<<=log_h;
		return h;
	}

	bool fastReconstructRational(Element& n, Element& d, const Element& x, const Element& m, const Element& d_bound ) {
	       	
		size_t log_m = m.bitsize()-1; //true unless m = 2^k
		size_t log_bound = d_bound.bitsize()-1; // true unless d_bound = 2^l, repairable anyway
		Element ai, bi, ci, di;
		ai=1;bi=0;ci=0;di=1;
		Element bound;  _Z.powtwo(bound, log_bound);

		cur_ri = m;
		cur_rinext = x;
		_Z.quo(cur_ainext, cur_ri, cur_rinext);
		cur_qinext = cur_ainext;//  ri/ri_next

		if (!fastEEA (ai,bi,ci,di,m,log_m,x,bound, log_bound)) return false;

		int K=0;

		if (cur_rinext > 0) {
			while (cur_ainext <= d_bound) {
				++K;
				Element tmp;
				tmp = ai;
				ai = cur_ainext;
				bi = tmp;
				tmp = ci;
				_Z.axpy(ci, cur_qinext,tmp, di);
++C.mul_counter;
				di = tmp;

				tmp = cur_rinext;
				_Z.axpy(cur_rinext, -cur_qinext, tmp, cur_ri);
				if (cur_rinext==0 ) {
					break;
				}
				cur_ri = tmp;
++C.mul_counter;
				_Z.quo(cur_qinext, cur_rinext, cur_ri);
++C.div_counter;
				_Z.axpy(cur_ainext, ai, cur_qinext, bi);
++C.mul_counter;			
			}
		}
		
		_Z.mul(n,x, ai);
		_Z.axmyin(n,m,ci);
		//n = x_in*ai-m*ci;
		d = ai;

		Element _gcd;
		if (_Z.gcd(_gcd, n,d) > 1) {
		/* solution not implemented */
			cout << "Solution pair is not coprime\n";
			return false;
		}

		return true;
	}

	//void prevEEA(Element& aprev, Element& bprev,const Element& ai, const Element& bi) {	
	//}

	void prevEEA(Element& aprev, Element& bprev, Element& cprev, Element& dprev, 
		     const Element& ai, const Element& bi, const Element& ci, const Element& di) {
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
		_Z.quo(qi, ai, bi);
++C.div_counter;
		
		_Z.axpy(bprev, -qi, bi,ai);
		_Z.axpy(dprev, -qi, di,ci);
++C.mul_counter;++C.mul_counter;
		Element tmp;
		tmp = cur_ri;
		_Z.axpy(cur_ri, tmp, qi, cur_rinext);
		//cur_ri = cur_ri *qi + cur_rinext;
		cur_rinext = tmp;
		cur_ainext = ai;
		cur_qinext = qi;
++C.mul_counter;

	}

	/* extended Euclidean Algorithm */
	bool classicEEA(Element& ai, Element& bi, Element& ci, Element& di, const Element& r0, const Element& r1, const Element& bound) {

		Element ri, rinext;
		ri = r0; rinext = r1;
		ai =di = 1;
		bi =ci = 0;

		Element ainext,cinext;
		Element qinext;
		_Z.quo(qinext,ri,rinext);
++C.div_counter;
		while (1) {
			//++i;
			//ainext = ai*qinext + bi;
			_Z.axpy(ainext, ai,qinext,bi);
++C.mul_counter;
			if (bound < ainext) {
				cur_ri = ri;
				cur_rinext = rinext;
				cur_ainext = ainext;
				cur_qinext = qinext;
				//cout << "->E:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "\n";
				return true;
			}
			else {
			       	bi = ai;
				ai = ainext;
				Element temp = ci;
				_Z.axpy(ci, temp,qinext,di);
++C.mul_counter;
				di = temp;

				temp = ri;
				ri=rinext;
				_Z.axpy(rinext, -qinext,ri,temp);//rinext = ri-qinext*rinext
++C.mul_counter;
				if (rinext==0) {
					//ainext = infinity
				        cur_ri = ri;
				        cur_rinext = rinext;
				        cur_ainext = r0+1;
					cur_qinext = cur_ainext;//infinity
					//cout << "->E:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "\n";
					return true;
				}
				_Z.quo(qinext,ri,rinext);
++C.div_counter;				
			}
		}		
	}
	
	/* log(m)-1 <= d < log(m) ; d >=1
	 * 2^h=powh=bound>=2, h <= d
	 */
	bool fastEEA(Element& ai,Element& bi, Element& ci, Element& di, const Element& m, const size_t d, const Element& n, const Element& powh, const size_t& h) {
		
		//cout << m << " " << n << " " << powh << "\n" << flush;
		
		ai=Element(1);
		di=Element(1);
		bi=Element(0);
		ci=Element(0); //Q(0)=Id

		if (m==n) {
			cur_ri = m;
			cur_rinext = 0;
			cur_ainext = m+1;//infinity
			cur_qinext = m+1;
			ai=1; bi=1; ci=1; di=0;
			return true;
		}

		if (powh < 1) return false; //should not happen
		if (n < 2) {
			//cout << "n=1\n";
			                 //n==1 -> Q1=(m,1//1,0)
		                         //n==0 -> Q1=infinity
			cur_ri = m;
			cur_rinext = n;
			cur_ainext = m;//infinity
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
			//quo(qinext,m,n); 
			//if (qinext==1) { return (1,1,1,0) or (1,0,0,1)
			if ((n << 1) > m) {
				cur_ri = n;
				cur_rinext = m-n;
				_Z.quo(cur_qinext, cur_ri, cur_rinext);
++C.div_counter;
				cur_ainext = cur_qinext + 1;
				ai=1; bi=1; ci=1; di=0;
			} else {
				cur_ri = m;
				cur_rinext = n;
				_Z.quo(cur_ainext,m,n);
++C.div_counter;
				cur_qinext = cur_ainext;
			}
			return true;
		}

		if (n < _threshold) {       //what about m?
			return classicEEA(ai,bi,ci,di,m,n,powh);
		}

		size_t log_n = n.bitsize()-1;

		//cout << d << " " << log_n << " "  << h << "\n";

		if (2*h+1 < d) {
			//cout << "choice1";
		  //if (2*h+1 <= log_n) {
		  if (true) {	
			//cout << ".1\n";
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
			
			_Z.mul(cur_ri,m,di);
			_Z.axmyin(cur_ri,n,bi);
			_Z.mul(cur_rinext,n,ai);
			_Z.axmyin(cur_rinext,m,ci);
C.mul_counter+=4;
			//ri = m*di - n*bi;
			//rinext = -m*ci + n*ai; 
 
			if (cur_ri < 0) {
				cur_ri = -cur_ri;
				cur_rinext = -cur_rinext;
			}
			if (cur_rinext>0) {
				_Z.quo(cur_qinext,cur_ri,cur_rinext);
++C.div_counter;
				_Z.axpy(cur_ainext, ai,cur_qinext,bi); 
++C.mul_counter;
			} else {
				cur_ainext = m+1;//infinity
				cur_qinext = cur_ainext;
			}

			//cout << "\n2 backward steps\n";
			//cout << "->1:" << cur_ri << " " << cur_rinext <<" " <<  cur_ainext << "->\n";

		   } else {
			//cout << ".2\n";
			Element a1,b1,c1,d1;
			Element qi; 
			_Z.quo(qi,m,n);
++C.div_counter;
			Element ri = n;
			Element rinext;
			_Z.axpy(rinext,-qi,n,m);
			size_t k = qi.bitsize()-1;
			Element hh;
			if (qi<powh) {
				if (!fastEEA(a1,b1,c1,d1,ri,log_n,rinext,powtwo(hh,h-k-1),h-k-1)) return 0;
			} else {
				if (qi == powh) {
					ai=qi; bi= 1; ci=1;di=0;
					cur_ri = ri;
					cur_rinext = rinext;
					_Z.quo(cur_qinext,ri,rinext); 
++C.div_counter;
					_Z.axpy(cur_ainext,qi,cur_qinext,1);
++C.mul_counter;
					return 1;
				} else {
					ai=1; bi=0; ci=0; di=1;
					cur_ri = m;
					cur_rinext = n;
					cur_qinext = qi;
					cur_ainext = qi;
					return 1;
				}
			}
			_Z.axpy(ai,qi,a1,c1);
			_Z.axpy(bi,qi,b1,d1);
			ci = a1;
			di = b1;
C.mul_counter +=2;
			_Z.axpy(cur_ainext, ai,cur_qinext,bi);
++C.mul_counter; 
                        Element aprev,bprev,cprev,dprev;  

			//cout << "->2.5:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "->\n";

                        if (ai > powh) {//at most 2 forward steps, at most 2 backward steps 
					 //backward loop (max 2 steps)
				int K=0;
				while (ai > powh) {//one step back
					++K;
					prevEEA(aprev,bprev,cprev,dprev,ai,bi,ci,di);
                                        ai=aprev;bi=bprev;ci=cprev;di=dprev;
				}
				//cout << "\n" << K << " backward steps\n";
				//cout << "->End:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
                                return 1;
                        } 
		   }
		} else if (h <= d-1) {
			//cout << "choice2\n";
			Element a1,a2,b1,b2,c1,c2,d1,d2;
                        a1=a2=d1=d2=1;
                        b1=b2=c1=c2=0;

                        Element sqrth;
			size_t logsqrth;
                                //sqrt(sqrth,powh);
                        logsqrth = h >> (int) 1;
			powtwo(sqrth, logsqrth);
                                //logtwo(logsqrth, sqrth);
                                //_sqrth = h/2;

                                //if (!EEA(sgn, ri, rinext, qinext, a1,b1,c1,d1,m,d,n,sqrth, logsqrth)) return 0;
			if (!fastEEA(a1,b1,c1,d1,m,d,n,sqrth, logsqrth)) return 0;

			Element ri = cur_ri;
			Element rinext = cur_rinext;

			//ri = m*d1;
			//Integer::axmyin(ri,n,b1);
			////ri = m*d1 - n*b1;

			//rinext = n*a1;
			//Integer::axmyin(rinext,m,c1);
			////rinext = -m*c1 + n*a1;
			if (ri < 0) {
				ri = -ri;
				rinext = -rinext;
			} 				

			size_t log_m;
			//if (ri !=m) //logtwo(log_m,ri);
			if ((rinext > 0) && (cur_ainext <= powh)){
				log_m = rinext.bitsize()-1;
				Element m2, n2;
				m2 = rinext;
				_Z.axpy(n2, -cur_qinext, rinext, ri);
++C.mul_counter;
				/* compute Q(i+1) */
				Element tmp = a1;
				a1 = cur_ainext;
				b1 = tmp;
				tmp = c1;
				_Z.axpy (tmp, cur_qinext, c1, d1);
++C.mul_counter;
				d1 = c1;
				c1 = tmp;
					
				size_t k = a1.bitsize()-1 ;
				if (rinext >0) {
					if (a1 < powh) {
						if (!fastEEA(a2,b2,c2,d2,m2,log_m,n2, powtwo(sqrth,h-k-1), h-k-1)) return 0;
					} else {
						ai = a1; bi = b1; ci=c1; di = d1;
						cur_ri = m2;
						cur_rinext = n2;
						_Z.quo(cur_qinext,m2,n2);
++C.div_counter;
						_Z.axpy(cur_ainext,a1,cur_qinext,b1);
++C.mul_counter;
						return 1;
					}
				} else {
					ai = a1; bi = b1; ci=c1; di = d1; 
					cur_ri = m2;
					cur_rinext = n2;
					cur_qinext = m+1;
					cur_ainext = m+1;
					return 1;
				}
			} else {
				ai = a1; bi = b1; ci=c1; di = d1;
				//cout << "End2\n";
				return 1;
			}

			_Z.mul(ai,b1,c2);
			_Z.axpyin(ai,a1,a2);
                               //aistar = a1*a2 + b1*c2;
			_Z.mul(bi,b1,d2);
			_Z.axpyin(bi, a1,b2);
                               //bistar = a1*b2 + b1*d2;
			_Z.mul(ci,d1,c2);
			_Z.axpyin(ci, c1,a2);
                               //cistar = c1*a2 + d1*c2;
			_Z.mul(di,d1,d2);
			_Z.axpyin(di, c1,b2);
                               //distar = c1*b2 + d1*d2;
C.mul_counter+=8;

			_Z.mul(cur_ri,m,di);
			_Z.axmyin(cur_ri,n,bi);
			_Z.mul(cur_rinext,n,ai);
			_Z.axmyin(cur_rinext,m,ci);
C.mul_counter+=4;
			//ri = m*di - n*bi;
			//rinext = -m*ci + n*ai;
			if (cur_ri < 0) {
				cur_ri = -cur_ri;
				cur_rinext = -cur_rinext;
			}
			if (cur_rinext>0) {
				_Z.quo(cur_qinext,cur_ri,cur_rinext);
++C.div_counter;
				_Z.axpy(cur_ainext, ai,cur_qinext,bi); 
++C.mul_counter;
			} else {
				cur_ainext = m+1;//infinity
				cur_qinext = cur_ainext;
			}

			//cout << "->2:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "->\n";

                        Element aprev,bprev,cprev,dprev;
                        aprev=dprev=1;
                        bprev=cprev=0;

			//ai=aistar; bi = bistar; ci=cistar; di=distar;

                        if (ai > powh) {//at most 2 forward steps, at most 2 backward steps 
					//backward loop (max 2 steps)
				int K=0;
				while (ai > powh) {//one step back
					++K;
					prevEEA(aprev,bprev,cprev,dprev,ai,bi,ci,di);
                                        ai=aprev;bi=bprev;ci=cprev;di=dprev;
				}
				//cout << "\n" << K << " backward steps\n";
				//cout << "->End:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
                                return 1;
                        }                               

		} else {//h=d, h = d+1;
			//cout << "choice3\n";
			Element hh = powh;
			size_t log_hh = h;
			hh = powh >> 1;
			while (log_hh > d-1) {
				--log_hh;
				hh >>= 1;
			}
			if (!fastEEA(ai,bi,ci,di,m,d,n,hh, d-1)) return 0;
			
			//cout << "->3:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "->\n";
		}
		
		Element ri = cur_ri;
		Element rinext = cur_rinext;
		Element qinext = cur_qinext;
		//_Z.mul(ri,m,di);
		//_Z.axmyin(ri,n,bi);
		//_Z.mul(rinext,n,ai);
		//_Z.axmyin(rinext,m,ci);
                ////ri = m*di - n*bi;
                ////rinext = -m*ci + n*ai;
                if (rinext==0) {
			//cout << "\n0 forward steps\n";
			//cout << "->End1:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
			return 1;
		}
                //if (ri < 0) {
		//	ri = -ri;
                //        rinext = -rinext;
                //}
                //quo(qinext,ri,rinext);
				
		if (qinext <=0) {
			//cout << "ERROR sth went very very wrong:" ;
                        //cout << "m:" << m << " n:" << n << " h:" << powh << "\n";
                        //cout << ai << " " << bi << "\n" << ci << " " << di <<"\n"; //getchar();
                        return 0;
                }

                Element ainext, binext, cinext,dinext;
                ainext=dinext=1;
                binext=cinext=0;
		ainext = cur_ainext;
		int K=-1;
                while (1) {
			++K;
			//_Z.axpy(ainext, ai,qinext,bi);
			//ainext = cur_ainext;
			if (powh < ainext) {
				cur_ri = ri;
				cur_rinext = rinext;
				cur_ainext = ainext;
				cur_qinext = qinext;
				//cout << "\n" << K << " forward steps\n";
				//cout << "->End2:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
				return 1;
			}
			//else if ((ah==1) && ((int)ah << (int)h) != ainext) return 1;
			else {

				bi = ai;
				ai = ainext;
				Element temp = ci;
				//ci = ci*qinext + di;
				_Z.axpy(ci, ci,qinext,di);
++C.mul_counter;
				di = temp;

				temp = ri;
				ri=rinext;
				//rinext=qinext*ri+temp;
				//rinext = temp - qinext*ri;
				_Z.axpy(rinext, -qinext,ri,temp);
++C.mul_counter;
				if (rinext==0) {
					cur_ri = ri;
					cur_rinext = rinext;
					cur_ainext = m+1;
					cur_qinext = m+1;
					//cout << "\n" << K << " forward steps\n";
					//cout << "->End3:" << cur_ri << " " << cur_rinext << " " << cur_ainext << "\n";
					return 1;
				}
                                        //++counter;
				_Z.quo(qinext,ri,rinext);
++C.div_counter;
                                _Z.axpy(ainext, ai,qinext,bi);
++C.mul_counter; 
			}
		}

                return 0;
	}
};

template <class Ring>
class WangClassicRationalReconstruction: public RationalReconstructionBase<Ring>
{
	bool _reduce;
	bool _recursive;
	typedef typename Ring::Element Element;
public:

	WangClassicRationalReconstruction(const Ring& Z, const bool reduce = true, const bool recursive = false): RationalReconstructionBase<Ring>(Z) {
		_reduce = reduce; _recursive = recursive;
	}

	WangClassicRationalReconstruction<Ring> (const WangClassicRationalReconstruction<Ring>& RR): RationalReconstructionBase<Ring>(RR._Z) {
		_reduce = RR._reduce; _recursive = RR._recursive;
	}	

	~WangClassicRationalReconstruction() {}
	       

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m) {
		Element a_bound; _Z.sqrt(a_bound, m/2);
		bool res = reconstructRational(a,b,x,m,a_bound);
		res = res && (b <= a_bound);
		return res;
	}

	bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound) {	
		bool res;

		if (x == 0) {
		    a = 0;
		    b = 1;
		} else {
			bool res = ratrecon(a,b,x,m,a_bound);
			if (_recursive)
			for(Element newbound = a_bound + 1; (!res) && (newbound<x) ; ++newbound)
			    res = ratrecon(a,b,x,m,newbound);
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
	/* m> x> 0 */
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
++C.div_counter;		    

		    u = a;
		    a = r0;  	
		    r0 = u;	// r0 <-- num
                                                           
		    _Z.axmyin(a,u,q); // num <-- r0-q*num
++C.mul_counter;
		    //if (a == 0) return false;
		    
		    u = b;
		    b = t0;  	
		    t0 = u;	// t0 <-- den
		               
		    _Z.axmyin(b,u,q); // den <-- t0-q*den
++C.mul_counter;
		} 

//                        if (den < 0) {
//                                _Z.negin(num);
//                                _Z.negin(den);
//                        }

		if ((a>0) && (_reduce)) {

			// [GG, MCA, 1999] Theorem 5.26
			// (ii)
		    Element gg;
++C.gcd_counter;
		    if (_Z.gcd(gg,a,b) != 1) {
			
			Element ganum, gar2;
			for( q = 1, ganum = r0-a, gar2 = r0 ; (ganum >= a_bound) || (gar2<a_bound); ++q ) {
			    ganum -= a;
			    gar2 -= a;
			}
					
			//_Z.axmyin(r0,q,a);
			r0 = ganum;
			_Z.axmyin(t0,q,b);
++C.mul_counter;++C.mul_counter;			
			if (t0 < 0) {
			    a = -r0;
			    b = -t0;
			} else {
			    a = r0;
			    b = t0;
			}
			
//                                if (t0 > m/k) {
			if ((double)b > (double)m/(double)a_bound) {
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


} //namespace LinBox
#endif

