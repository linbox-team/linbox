#ifndef _FIELD_SPLA_PADIC_ABSTRACT
#define _FIELD_SPLA_PADIC_ABSTRACT

#include <iostream>
#include "LinBox/field_abstract.h"
#include "LinBox/element_spla_padic.h"
#include <cassert>

namespace LinBox{

  class RandIter_Padic;

  class Field_spla_Padic:public Field_abstract{
  public:
    
    typedef Element_spla_Padic element;

    typedef RandIter_Padic randIter;
 
    virtual ~Field_spla_Padic(){}

    Field_spla_Padic(long int p, int accurate):prime(p), n_digit(accurate),base((long int)pow(p, accurate)){}

    Field_spla_Padic(const Field_spla_Padic& f):prime(f.prime),base(f.base),n_digit(f.n_digit){}

    virtual Field_abstract* clone() const{ return new Field_spla_Padic(*this);}
    
    virtual Field_abstract& operator= (const Field_abstract& F){return *this;}
   
    virtual Element_abstract& init(Element_abstract& x, const integer& value) const
      {
	element& ref=static_cast<element&>(x);
	ref.exp=0;
	integer y=value;
	if(y==0) 
	  {
	    ref.rep=0;
	  }
	else 
	  {
	    bool signal=true;
	    if(y<0) {signal=false;y=-y;}
	    while(y%prime==0) {y=y/prime;++ref.exp;}
	    ref.rep=y%base;
	    if(!signal)
		ref.rep=base-ref.rep;
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
	      x=static_cast<const element&>(y).rep;
	      for(int i=0;i<static_cast<const element&>(y).exp;++i)
		x*=prime;
	    }
	
	return x;
      }

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

	if(refy.exp>refz.exp)
	  {
	    int diff=refy.exp-refz.exp;
	    if(diff>=n_digit)
	      refy.rep=0;
	    else
	      {
		for(int i=0;i<diff;++i)
		  refy.rep=refy.rep*prime;
		refy.rep=refy.rep%base;
	      }
	    refy.exp=refz.exp;
	  }
      
	else 
	  if(refz.exp>refy.exp) 
	    {
	      int diff=refz.exp-refy.exp;
	      if(diff>=n_digit)
		refz.rep=0;
	      else
		{
		  for(int i=0;i<diff;++i)
		    refz.rep=refz.rep*prime;
		  refz.rep=refz.rep%base;
		}
	      refz.exp=refy.exp;
	    }

       
	refx.rep=refy.rep+refz.rep;
	if(refx.rep>=base)
	  refx.rep=refx.rep-base;

	refx.exp=refy.exp;
	
	if(refx.rep!=0)
	  while((refx.rep%prime)==0)
	    {
	      refx.exp++;
	      refx.rep/=prime;
	    }

	return x;
      }

  
    virtual Element_abstract& sub(Element_abstract& x, const Element_abstract& y, const Element_abstract& z) const
      {
	if(areEqual(y, z))
	  return init(x, 0);

	element& refx=static_cast<element&>(x);
	element refy=static_cast<const element&>(y);
	element temp;
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
	    refx.exp=refy.exp;
	    refx.rep=base-refy.rep;
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
	
	refx.rep=(refy.rep*refz.rep)%base;
	refx.exp=refy.exp+refz.exp;

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
		
	long int inverse;
	inv_help(inverse, refz.rep);
	refx.rep=(refy.rep*inverse)%base;
	refx.exp=refy.exp-refz.exp;
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
	element tmp;
	mul(tmp,a,x);
	return add(r,tmp,y);
      }

  
    virtual bool isZero(const Element_abstract& x) const
      {
	return (static_cast<const element&>(x).rep==0);
      }
    
    virtual bool isOne(const Element_abstract& x) const
      {
	return ((static_cast<const element&>(x).rep==1)&&(static_cast<const element&>(x).exp==0));
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
    long int get_prime() const
      {
	return prime;
      }

    int get_digit()
      {
	return n_digit;
      }

    void dec_exp(element& x, int i) const {x.exp-=i;}
  
  private:
    
    const long int prime;
    const long int base;
    const long int n_digit;
    
    long int& inv_help(long int& x, const long int& y) const
      {
	assert(y!=0);
	long int gcd,t;
	EEA(gcd,x,t,y,base);
	while(x<0)
	  x=base+x;
      
	return x;

	    
      }
	
    void EEA(long int& gcd,long int& s, long int& t, const long int& f, const long int& g) const
      {
	long int s0=1;
	long int s1=0;
	long int t0=0;
	long int t1=1;
	long int r0=f;
	long int r1=g;
	while(r1!=0)
	  {
	    int q=r0/r1;
	    int tmp=r0;
	    r0=r1;
	    r1=tmp-q*r0;
	    tmp=s0;
	    s0=s1;
	    s1=tmp-q*s0;
	    tmp=t0;
	    t0=t1;
	    t1=tmp-q*t0;
	  }
	gcd=r0;
	s=s0;
	t=t0;
      }
  };
}
#endif
