#ifndef _ELEMENT_SPLA_PADIC_ABSTRACT
#define _ELEMENT_SPLA__PADIC_ABSTRACT

#include "LinBox/element_abstract.h"
#include <iostream>

namespace LinBox{
  class Element_spla_Padic: public Element_abstract
    {
      friend std::ostream& operator<<(std::ostream& out, const Element_spla_Padic& x)
	{
	  out<<x.rep<<" ";
	  out<<x.exp;
	  return out;
	}

    public:

      Element_spla_Padic(){}
 
      Element_spla_Padic(const Element_spla_Padic& x):rep(x.rep),exp(x.exp){}

      virtual Element_abstract* clone() const
	{
	  Element_spla_Padic* p=new Element_spla_Padic;
	  p->rep=rep;
	  p->exp=exp;
	  return p;
	}

      Element_spla_Padic& operator= (const Element_spla_Padic& a)
      {
	this->exp=a.exp;
	this->rep=a.rep;
	return *this;
      }

      virtual Element_abstract& operator=(const Element_abstract& a)
      {
	this->rep=static_cast<const Element_spla_Padic&>(a).rep;
	this->exp=static_cast<const Element_spla_Padic&>(a).exp;
	return *this;
      }
  
      virtual ~Element_spla_Padic()
	{}
	  
     public:
      long int rep;
      int exp;
    };
}
#endif
