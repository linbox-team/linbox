#ifndef _ELEMENT_SPHA_PADIC_ABSTRACT
#define _ELEMENT_SPHA_PADIC_ABSTRACT

#include "LinBox/element_abstract.h"
#include <vector>
#include <iostream>
#include "LinBox/integer.h"

//this is for small prime and high accuracy.

namespace LinBox{
  class Element_spha_Padic: public Element_abstract
    {
      friend std::ostream& operator<<(std::ostream& out, const Element_spha_Padic& x)
	{
	out<<x.rep<<"\t"<<x.exp;
	  return out;
	}

    public:

      Element_spha_Padic(){}

      Element_spha_Padic(const Element_spha_Padic& x):rep(x.rep),exp(x.exp){}

      virtual Element_abstract* clone() const
	{
	  Element_spha_Padic* p=new Element_spha_Padic;
	  p->rep=rep;
	  p->exp=exp;
	  return p;
	}

      Element_spha_Padic& operator= (const Element_spha_Padic& a)
      {
	this->exp=a.exp;
	this->rep=a.rep;
	return *this;
      }

      virtual Element_abstract& operator=(const Element_abstract& a)
      {
	this->rep=static_cast<const Element_spha_Padic&>(a).rep;
	this->exp=static_cast<const Element_spha_Padic&>(a).exp;
	return *this;
      }
  
      virtual ~Element_spha_Padic(){}
 
    public:
      integer rep;
      int exp;
    };
}
#endif
