#ifndef _ELEMENT_BP_PADIC_ABSTRACT
#define _ELEMENT_BP_PADIC_ABSTRACT

#include "LinBox/element_abstract.h"
#include <vector>
#include <iostream>

namespace LinBox{
  class Element_bp_Padic: public Element_abstract
    {
      friend std::ostream& operator<<(std::ostream& out, const Element_bp_Padic& x)
	{
	  for(std::vector<integer>::const_iterator p=x.rep.begin();p!=x.rep.end();++p)
	    out<<*p<<" ";
	  out<<x.exp;
	  return out;
	}

    public:

      Element_bp_Padic(){}

      Element_bp_Padic(unsigned int k){rep.resize(k);}

      Element_bp_Padic(const Element_bp_Padic& x):rep(x.rep),exp(x.exp){}

      virtual Element_abstract* clone() const
	{
	  Element_bp_Padic* p=new Element_bp_Padic;
	  p->rep=rep;
	  p->exp=exp;
	  return p;
	}

      Element_bp_Padic& operator= (const Element_bp_Padic& a)
      {
	this->exp=a.exp;
	this->rep=a.rep;
	return *this;
      }

      virtual Element_abstract& operator=(const Element_abstract& a)
      {
	this->rep=static_cast<const Element_bp_Padic&>(a).rep;
	this->exp=static_cast<const Element_bp_Padic&>(a).exp;
	return *this;
      }
  
      virtual ~Element_bp_Padic()
	{
	  rep.clear();
	}
 
    public:
      std::vector<integer> rep;
      int exp;
    };
}
#endif
