#ifndef _FIELD_PADIC_ABSTRACT
#define _FIELD_PADIC_ABSTRACT

#include <iostream>
#include "LinBox/field_abstract.h"
#include "LinBox/element_padic.h"
#include <cassert>

namespace LinBox{

  class RandIter_Padic;

  class Field_Padic:public Field_abstract{
  public:
    // Element type.
    typedef Element_Padic element;

    typedef RandIter_Padic randIter;
 
    
    // virtual Destructor
    virtual ~Field_Padic(){}

   
    Field_Padic(unsigned long int prime,unsigned int accurate):base(prime), n_digit(accurate){}

    Field_Padic(const Field_Padic& f):base(f.base),n_digit(f.n_digit){}

    virtual Field_abstract* clone() const{ return new Field_Padic(*this);}
    
    virtual Field_abstract& operator= (const Field_abstract& F){return *this;}
   
    virtual Element_abstract& init(Element_abstract& x, const integer& value) const
      {
	element& ref=static_cast<element&>(x);
	ref.exp=0;
	ref.rep.resize(n_digit);
	integer y=value;
	if(y==0) 
	  {
	    ref.rep.clear();
	    ref.rep.insert(ref.rep.begin(),n_digit,0);
	    ref.exp=n_digit;
	  }
	else 
	  {
	    bool signal=true;
	    if(y<0) {signal=false;y=-y;}
	    while(y%base==0) {y=y/base;++ref.exp;}
	    for(int i=0;i<n_digit;++i) {ref.rep[i]=y%base;y=y/base;}
	    if(!signal)
	      {
		ref.rep[0]=base-ref.rep[0];
		for(int i=1;i<n_digit;++i)
		  {ref.rep[i]=base-1-ref.rep[i];}
	      }
	  }
	
	return x;
      }


    /** Conversion of field element to an integer.
     *If y's exponent is non-negative, return its representation, 
     *else just return x unchanged.
     */
    virtual integer& convert(integer& x, const Element_abstract& y) const
      {
	if(static_cast<const element&>(y).exp>=0) 
	    {
	      x=0;
	      int i;
	      for(i=n_digit-1;i>=0;--i)
		x=x*base+static_cast<const element&>(y).rep[i];
	      for(i=0;i<static_cast<const element&>(y).exp;++i)
		x*=base;
	    }
	
	return x;
      }

    /** Assignment of one field element to another.
     * @return reference to x
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    virtual Element_abstract& assign(Element_abstract& x, const Element_abstract& y) const
      {
	static_cast<element&>(x).rep=static_cast<const element&>(y).rep;
	static_cast<element&>(x).exp=static_cast<const element&>(y).exp;
	return x;
      }
  
    virtual integer& cardinality(integer& c) const {c=-1; return c;}
   
    virtual integer& characteristic(integer& c) const {c=0; return c;}

    virtual bool areEqual(const Element_abstract& x, const Element_abstract& y) const
      {
	return ((static_cast<const element&>(x).rep==static_cast<const element&>(y).rep)&&(static_cast<const element&>(x).exp==static_cast<const element&>(y).exp));
      }

    virtual Element_abstract& add(Element_abstract& x,const Element_abstract& y, const Element_abstract& z) const
      {	
	if(isZero(y))
	  return assign(x,z);

	if(isZero(z))
	  return assign(x,y);
	
	element& refx=static_cast<element&>(x);
	element  refy=static_cast<const element&>(y);
	element  refz=static_cast<const element&>(z);

	refx.rep.resize(n_digit);

	if(refy.exp>refz.exp)
	  {
	    int i;
	    for(i=n_digit-1;i>=refy.exp-refz.exp;--i)
	      refy.rep[i]=refy.rep[i-refy.exp+refz.exp];
	    for(i=0;(i<refy.exp-refz.exp)&&(i<n_digit);++i)
	      refy.rep[i]=0;
	    refy.exp=refz.exp;
	  }
      
	else 
	  if(refz.exp>refy.exp) 
	    {
	      int i;
	      for(i=n_digit-1;i>=refz.exp-refy.exp;--i)
		refz.rep[i]=refz.rep[i-refz.exp+refy.exp];
	      for(i=0;(i<refz.exp-refy.exp)&&(i<n_digit);++i)
		refz.rep[i]=0;
	      refz.exp=refy.exp;
	    }

        int r=0;
	for(int i=0;i<n_digit;++i)
	  {
	    refx.rep[i]=refy.rep[i]+refz.rep[i]+r;
	    r=refx.rep[i]/base;
	    refx.rep[i]%=base;
	  }
      
	refx.exp=refy.exp;
	int k=0;
	while((refx.rep[k]==0)&&(k<n_digit)) ++k;
      
	if(k>0)
	  {
	    int i;
	    for(i=0;i<n_digit-k;++i)
	      refx.rep[i]=refx.rep[i+k];
	    for(i=n_digit-k;i<n_digit;++i)
	      refx.rep[i]=0;
	    refx.exp+=k;
	  }

	return x;
      }

  
    virtual Element_abstract& sub(Element_abstract& x, const Element_abstract& y, const Element_abstract& z) const
      {
	if(areEqual(y, z))
	  return init(x, 0);

	element& refx=static_cast<element&>(x);
	element refy=static_cast<const element&>(y);
	element temp(n_digit);
	neg(temp,z);
      
	add(refx,refy,temp);
	return x;
      }

  
    virtual Element_abstract& neg(Element_abstract& x, const Element_abstract& y) const
      {
	if(!isZero(y))
	  { 
	    element& refx=static_cast<element&>(x);
	    element refy=static_cast<const element&>(y);
	    refx.rep.resize(n_digit);
	    refx.exp=refy.exp;
	    refx.rep[0]=base-refy.rep[0];
	    for(int i=1;i<n_digit;++i)
	      refx.rep[i]=base-1-refy.rep[i];
	  }

	else
	  init(x, 0);
      
	return x;
      }

 
    virtual Element_abstract& mul(Element_abstract& x,const Element_abstract& y, const Element_abstract& z) const
      {
	if(isZero(y) || isZero(z))
	  return init(x,0);	 
	 
	if(isOne(y))
	  return assign(x,z);

	if(isOne(z))
	  return assign(x,y);
	  
	element& refx=static_cast<element&>(x);
        element refy=static_cast<const element&>(y);
	element refz=static_cast<const element&>(z);

	refx.rep.resize(n_digit);

	for(int k=0;k<n_digit;++k)
	  {
	    refx.rep[k]=0;
	    for(int i=0;i<=k;++i)
	      refx.rep[k]+=refy.rep[i]*refz.rep[k-i];
	  }
	
	refx.exp=refy.exp+refz.exp;

	unsigned long int r=0;
	for(int i=0;i<n_digit;++i)
	  {
	    refx.rep[i]=refx.rep[i]+r;
	    r=refx.rep[i]/base;
	    refx.rep[i]%=base;
	  }

      return x;
      }

   
    virtual Element_abstract& div(Element_abstract& x,const Element_abstract& y, const Element_abstract& z) const
      {
	if(isZero(y)&&(!isZero(z)))
	  return init(x,0);

	if(isOne(z))
	  return assign(x,y);

	if(areEqual(y,z))
	  return init(x,1);

	element& refx=static_cast<element&>(x);
	element refy=static_cast<const element&>(y);
	element refz=static_cast<const element&>(z);
	refx.exp=refy.exp-refz.exp;
	
	refx.rep.resize(n_digit);
	
	assert(!isZero(z));
	
	long int inverse=inv_help(refz.rep[0]);

	for(int k=0;k<n_digit;++k)
	  {
	    refx.rep[k]=(refy.rep[k]*inverse)%base;
	    
	    int lend=0;
	    for(int i=k;i<n_digit;++i)
	      {
		refy.rep[i]+=-refz.rep[i]*refx.rep[k]-lend;

		if(refy.rep[i]>=0) lend=0;
		else 
		  {
		    refy.rep[i]=-refy.rep[i];
		    lend=refy.rep[i]/base;
		    refy.rep[i]=refy.rep[i]%base;
		    if(refy.rep[i]!=0) {refy.rep[i]=base-refy.rep[i];++lend;}
		  }
	      }
	  }
	return x;
      }

    
    virtual Element_abstract& inv(Element_abstract& x, const Element_abstract& y) const
      {
	if(isOne(y))
	  return assign(x,y);

	element temp;
	init(temp, 1);
	return div(x,temp, y);
      }

    virtual Element_abstract& axpy(Element_abstract& r, const Element_abstract& a, const Element_abstract& x, const Element_abstract& y) const
      {
	element tmp(n_digit);
	mul(tmp,a,x);
	return add(r,tmp,y);
      }

  
    virtual bool isZero(const Element_abstract& x) const
      {
	for(int i=0; i<n_digit;++i)
	  if(static_cast<const element&>(x).rep[i]!=0) return false;
	return true;
      }

    virtual bool isOne(const Element_abstract& x) const
      {
	if((static_cast<const element&>(x).rep[0]!=1)||(static_cast<const element&>(x).exp!=0)) return false;
	else
	  for(int i=1;i<n_digit;++i)
	    if(static_cast<const element&>(x).rep[i]!=0) return false;
	return  true;
      }

    virtual Element_abstract& addin(Element_abstract& x, const Element_abstract& y) const
      {
	return add(x,x,y);
      }

  
    virtual Element_abstract& subin(Element_abstract& x, const Element_abstract& y) const
      { 
	return sub(x,x,y);
      }

 
    virtual Element_abstract& mulin(Element_abstract& x, const Element_abstract& y) const
      {
	return mul(x,x,y);
      }

    virtual Element_abstract& divin(Element_abstract& x, const Element_abstract& y) const
      {
	return div(x,x,y);
      }

    virtual Element_abstract& negin(Element_abstract& x) const
      {
	return neg(x,x);
      }

  
    virtual Element_abstract& invin(Element_abstract& x) const
      {
	return inv(x,x);
      }

  
    virtual Element_abstract& axpyin(Element_abstract& r,const Element_abstract& a,const Element_abstract& x) const
      {
	return axpy(r,a,x,r);
      }

   
    virtual std::ostream& write(std::ostream& os) const{return os;}

  
    virtual std::istream& read(std::istream& is) {return is;}

  
    virtual std::ostream& write(std::ostream& os, const Element_abstract& x) const
      {integer a; 
      os<<convert(a,x); 
      return os;}

  
    virtual std::istream& read(std::istream& is, Element_abstract& x) const
      {int data;
      is>>data;
      init(x,data);
      return is;}

    //two more help function to get information about the field.
    unsigned long int get_base()
      {
	return base;
      }

    unsigned int get_digit()
      {
	return n_digit;
      }

    void dec_exp(element& x, int i) {x.exp-=i;}
  private:
    
    const unsigned long int base;
    const unsigned long int n_digit;
    long int inv_help(const long int& y) const
      {
	assert(y!=0);
	for(long int int_try=1;int_try<base;++int_try)
	  if((int_try*y)%base==1) return int_try;
      }
	
  };
}
#endif
