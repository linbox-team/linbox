#ifndef _ELEMENT_PADIC_ABSTRACT
#define _ELEMENT_PADIC_ABSTRACT

#include "LinBox/element_abstract.h"
#include <vector>
#include <iostream>

namespace LinBox{
  class Element_Padic: public Element_abstract
    {
      friend std::ostream& operator<<(std::ostream& out, const Element_Padic& x)
	{
	  for(std::vector<long int>::const_iterator p=x.rep.begin();p!=x.rep.end();++p)
	    out<<*p<<" ";
	  out<<x.exp;
	  return out;
	}

    public:

      Element_Padic(){}
 
      Element_Padic(unsigned int k){rep.resize(k);}

      Element_Padic(const Element_Padic& x):rep(x.rep),exp(x.exp){}

      virtual Element_abstract* clone() const
	{
	  Element_Padic* p=new Element_Padic;
	  p->rep=rep;
	  p->exp=exp;
	  return p;
	}

      Element_Padic& operator= (const Element_Padic& a)
      {
	this->exp=a.exp;
	this->rep=a.rep;
	return *this;
      }

      virtual Element_abstract& operator=(const Element_abstract& a)
      {
	this->rep=static_cast<const Element_Padic&>(a).rep;
	this->exp=static_cast<const Element_Padic&>(a).exp;
	return *this;
      }
  
      virtual ~Element_Padic()
	{
	  rep.clear();
	}
 
    public:
      std::vector<long int> rep;
      long int exp;
    };
}
#endif
