/* File: src/wrapper/by_library/Lidia/lidia_gfq.h
 * Author: Pascal Giorgi for the linbox group
 */


#ifndef _LIDIA_GFQ_
#define _LIDIA_GFQ_

//-----------------------------------
// Files of C/C++ library
#include <iostream>

//-----------------------------------
// Files of LiDIA library
#include "LiDIA/gf_element.h"

//-----------------------------------
// Files of LinBox library
#include "LinBox/integer.h"
#include "LinBox/lidia_randiter_gfq.h"

//------------------------------------


// Namespace in which all LinBox library code resides
namespace LinBox
{
using namespace LiDIA;


 /** This class define Galois Field GF(p^k) with p prime and inherits from 
  *  galois_field of LiDIA.
  */
     

  class lidia_gfq  : public galois_field 
  {
  public:

    /** Element type.
     *  This type is inherited from the LiDIA class gf_element
     */
    typedef gf_element  element;
    
  
    /** Random element generator which is define in the wrapper LIDIA_randiter
     */
    typedef lidia_randIter_gfq<lidia_gfq>  randIter;



    /** Default constructor of the field
     */
    lidia_gfq() {}



    /** Constructor from two integer p, k.
     *  A GF(p^k) field is construct throught 
     *  the constructor of LiDIA galois_field
     *  We need a double cast to pass integer arguments to the LiDIA constructor
     */
    lidia_gfq(const integer& p , const integer& k) :
      galois_field(static_cast<bigint>(static_cast<int>(p)), 
		   static_cast<lidia_size_t>(static_cast<int>(k))) {}
 


    /** Destructor
     */
    ~lidia_gfq() {}




    /** Assignment operator.
     * Assigns unparam_field object F to field.
     * @param  F unparam_field object.
     */
    lidia_gfq& operator=(const lidia_gfq& F)
      {return *this;}
    



    /** Initialization of field element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output field element x has already been 
     * constructed, but that it is not already initialized.
     * We also need to define the element over the field.
     * So what we always initialize the element with the zero field value.
     * If an integer different from zero is passed to the function the element
     * is initialized to a constant polynom of Z/pZ
     * @return reference to field element.
     * @param x field element to contain output (reference returned).
     * @param y integer.
     */
    element& init(element& x , const integer& y=0) const
      {
	 x.assign_zero(*this);
	 if (y!=0) 	
	  {
	    Fp_polynomial Pol;
	    integer p;
	    characteristic(p);
	    Pol.set_modulus(static_cast<bigint>(static_cast<int>(p)));
	    Pol.set_max_degree((x.get_field()).degree());
	    
	    
	    integer rem, quo,tmp=y;
	    for(lidia_size_t i=0;i<(x.get_field()).degree();i++)
	      {		
		quo=tmp/p;
		rem=tmp%p;
		tmp=quo;
		Pol.set_coefficient(static_cast<bigint>(static_cast<int>(rem)),i);
	      }
	    element * e=new element(x.get_field());	   
	    e->set_polynomial_rep(Pol);
	    x.assign(*e);
	  }

	return x;
      }
    
    

    /** Conversion of field base element to an integer.
     * This function assumes the output field base element x has already been
     * constructed, but that it is not already initialized.
     * As elements are represented by polynom the convert function return 
     * the valuation of polynom in characteristic by the Horner Method.
     * That keeps unicity of each element.
     * @return reference to an integer.
     * @param x integer to contain output (reference returned).
     * @param y constant field base element.
     */
    integer& convert(integer& x , const element& y ) const
      {
	bigint fx(0) , X((y.get_field()).characteristic());
	bigint tmp;
	
	 
	for(int i=static_cast<int>((y.get_field()).degree());i>0;i--)
	  {
	    (y.polynomial_rep()).get_coefficient(tmp,i);
	    fx=fx+tmp;
	    fx=fx*X;
	  }
	(y.polynomial_rep()).get_coefficient(tmp,0);
	fx= fx + tmp;

	long i;
	fx.longify(i);
	
	return x=*(new integer(i));
      }




    /** Assignment of one field element to another.
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x
     * @param  x field element (reference returned).
     * @param  y field element.
     */
    element& assign(element& x, const element& y) const
      {
	return x=y;
      }




    /** Cardinality.
     * Return integer representing cardinality of the field.
     * Returns  p^k.
     * @return constant reference to integer representing cardinality 
     *	       of the field.
     */
    integer& cardinality(integer& c) const 
      {
	long tmp;
	(number_of_elements()).longify(tmp);
	return c=*(new integer(tmp));
 
      }



    /** Characteristic.
     * Return integer representing characteristic of the field.
     * Returns p.
     * @return constant reference to integer representing characteristic 
     * 	       of the field.
     */
    integer& characteristic(integer& c) const
      {
	galois_field F(*this);
	long tmp;
	(F.characteristic()).longify(tmp);
	return c=*(new integer(tmp));
      }


    //@} Object Management
    
    /** @name Arithmetic Operations 
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized field elements will
     * give undefined results.
     */
    //@{
     
    /** Equality of two elements.
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return boolean true if equal, false if not.
     * @param  x field element
     * @param  y field element
     */
     bool areEqual(const element& x, const element& y) const
       {
	 return x==y;
       }


 /** Addition.
     * x = y + z
     * This function assumes all the field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     * @param  z field element.
     */
     element& add(element& x, const element& y, const element& z) const
       { return x = y + z; }
    
     /** Subtraction.
      * x = y - z
      * This function assumes all the field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      * @param  z field element.
      */
     element& sub(element& x, const element& y, const element& z) const
       { return x = y - z; }
     
     /** Multiplication.
      * x = y * z
      * This function assumes all the field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      * @param  z field element.
      */
     element& mul(element& x, const element& y, const element& z) const
       { return x = y * z; }
     
     /** Division.
      * x = y / z
      * This function assumes all the field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      * @param  z field element.
      */
     element& div(element& x, const element& y, const element& z) const
       { return x = y / z; }
     
     
     /** Additive Inverse (Negation).
      * x = - y
      * This function assumes both field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      */
     element& neg(element& x, const element& y) const
       {
	 negate(x,y);
	 return x;
       }



     /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both field elements have already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     * @param  y field element.
     */
     element& inv(element& x, const element& y) const
       {
	 invert(x,y);
	 return x;
       }
     
     

     /** Natural AXPY.
     * r  = a * x + y
     * This function assumes all field elements have already been 
     * constructed and initialized.
     * @return reference to r.
     * @param  r field element (reference returned).
     * @param  a field element.
     * @param  x field element.
     * @param  y field element.
     */
     element& axpy(element& r, 
		   const element& a,
		   const element& x, 
		   const element& y) const
       {
	 return r=a*x+y;
       }

      

     /** Zero equality.
     * Test if field element is equal to zero of field.
     * This function assumes the field element has already been 
     * constructed and initialized.
     * @return boolean true if equals zero of field, false if not.
     * @param  x field element.
     */
     bool isZero(const element& x) const 
       {
	 x.is_zero();
       }



     /** One equality.
     * Test if field element is equal to one of field.
     * This function assumes the field element has already been 
     * constructed and initialized.
     * @return boolean true if equals one of field, false if not.
     * @param  x field element.
     */
     bool isOne(const element& x) const 
       {
	 x.is_one();
       }


     
     /** Inplace Addition.
      * x += y
      * This function assumes both field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      */
     element& addin(element& x, const element& y) const { return x += y; }
     
     /** Inplace Subtraction.
      * x -= y
      * This function assumes both field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      */
     element& subin(element& x, const element& y) const { return x -= y; }
     
     /** Inplace Multiplication.
      * x *= y
      * This function assumes both field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      */
     element& mulin(element& x, const element& y) const { return x *= y; }
    
     /** Inplace Division.
      * x /= y
      * This function assumes both field elements have already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      * @param  y field element.
      */
     element& divin(element& x, const element& y) const { return x /= y; }
     
     
     /** Inplace Additive Inverse (Inplace Negation).
      * x = - x
      * This function assumes the field element has already been 
      * constructed and initialized.
      * @return reference to x.
      * @param  x field element (reference returned).
      */
     element& negin(element& x) const
       {
	 x.negate();
	 return x;
       }
     
     

     /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the field elementhas already been 
     * constructed and initialized.
     * @return reference to x.
     * @param  x field element (reference returned).
     */
     element& invin(element& x) const
       {
	 x.invert();
	 return x;
       }



     /** Inplace AXPY.
      * r  += a * x
      * This function assumes all field elements have already been 
      * constructed and initialized.
      * @return reference to r.
      * @param  r field element (reference returned).
      * @param  a field element.
      * @param  x field element.
      */
     element& axpyin(element& r, const element& a, const element& x) const
       {
	 r+=a*x;
       }

     //@} Inplace Arithmetic Operations


 /** @name Input/Output Operations */
    //@{
    
    /** Print field.
     * @return output stream to which field is written.
     * @param  os  output stream to which field is written.
     */
     ostream& write(ostream& os) const
       {
	 integer c;
	 return os<<"corps de Galois GF("<<
	    characteristic(c)<<"^"<<degree()<<")";
       }
     
     
     /** Read field.
      * @return input stream from which field is read.
      * @param  is  input stream from which field is read.
      */
     istream& read(istream& is) const
       {
	 return is ;
       }
     

     /** Print field element like a polynom.
     * @return output stream to which field element is written.
     * @param  os  output stream to which field element is written.
     * @param  x   field element.
     */
     ostream& write(ostream& os,const element& e) const
	{
	  return os<<e;
	}
     
     
     
     /** Read field element.
      * @return input stream from which field element is read.
      * @param  is  input stream from which field element is read.
      * @param  x   field element.
      */
     istream& read(istream& is, element& e) const
       {
	 return is>>e ;
       }

     //@}
      
 }; // class lidia_gfq

} // namespace LinBox


#endif // LIDIA_FIELD
