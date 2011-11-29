/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#ifndef __LINBOXX__RECONSTRUCTION_BASE_H__
#define __LINBOXX__RECONSTRUCTION_BASE_H__

#include <iostream>
#include <deque>

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
	OpCounter C;
	typedef typename Ring::Element Element;
	
	//OpCounter C;

	RationalReconstructionBase(const Ring& Z): _Z(Z) {}
	RationalReconstructionBase(const RationalReconstructionBase<Ring>& RR): _Z(RR._Z) {}

	virtual bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m)=0;

	virtual bool reconstructRational(Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound)=0;

	virtual ~RationalReconstructionBase() {}

	void write(ostream& is) {
		C.write(is);
	}
};

} //namespace LinBox
#endif

