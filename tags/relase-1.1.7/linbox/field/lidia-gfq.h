/* linbox/field/lidia-gfq.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */


#ifndef __LINBOX_lidia_gfq_H
#define __LINBOX_lidia_gfq_H

//-----------------------------------
// Files of C/C++ library
#include <iostream>

//-----------------------------------
// Files of LiDIA library
#include "LiDIA/gf_element.h"
#include "LiDIA/bigint.h"

//-----------------------------------
// Files of LinBox library
#include "linbox/integer.h"
#include <linbox/field/field-interface.h>
#include "linbox/randiter/lidia-gfq.h"
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <string>
#include <vector>

#endif

//------------------------------------


// Namespace in which all LinBox library code resides
namespace LinBox
{

	template <class Ring>
	struct ClassifyRing;

	class LidiaGfq;

	template<>
	struct ClassifyRing<LidiaGfq> {
		typedef RingCategories::ModularTag categoryTag;
	};

	/** @brief defines the  Galois Field GF(p<sup>k</sup>).

	* Inherits from galois\_field of LiDIA.
	 *  p must be prime.
	\ingroup field
	 */
	class LidiaGfq  : public LiDIA::galois_field, public FieldInterface 
	{
	public:

		/** Element type.
		 *  This type is inherited from the LiDIA class LiDIA::gf_element
		 */
		typedef LiDIA::gf_element  Element;
    
    
		/** Random element generator which is define in the wrapper LIDIA_randiter
		 */
		typedef LidiaGfqRandIter<LidiaGfq>  RandIter;
    


		/** Default constructor of the field
		 */
		LidiaGfq() {}


		/** @brief Constructor of GF(p^k).

		 *  A GF(p^k) field is constructed through 
		 *  the constructor of LiDIA LiDIA::galois_field
		 *  We need a double cast to pass integer arguments to the LiDIA constructor
		 */
		LidiaGfq(const integer& p , const integer& k=1) :
			LiDIA::galois_field(static_cast<LiDIA::bigint>(double(p)), 
				     static_cast<LiDIA::lidia_size_t>(int(k))) {}
     

		/** Copy constructor 
		 */
		LidiaGfq(const LidiaGfq& F) : LiDIA::galois_field(F) {}

#ifdef __LINBOX_XMLENABLED
		// XML LinBox::Reader constructor
		LidiaGfq(LinBox::Reader &R) : LiDIA::galois_field()
		{
			integer p, k;
			if(!R.expectTagName("field") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagName("finite") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(p)) return;
			R.upToParent();
			R.upToParent();

			if(R.getNextChild()) {
				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("extension") || !R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagNum(k)) return;

				R.upToParent();
				R.upToParent();
				R.getPrevChild();
			}
			else {
				k = Integer(1);
			}
			R.upToParent();

			LidiaGfq oth(p, k);
			*this = oth;

			return;

		}
#endif

    
		/** Destructor
		 */
		~LidiaGfq() {}




		/** Assignment operator.
		 * Assigns unparam_field object F to field.
		 * @param  F unparam_field object.
		 */
		LidiaGfq& operator=(const LidiaGfq& F)
		{
			return *this = F;
		}

		/** @name Object management
		 */

		//@{


		/** Initialization of field Element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field Element x has already been 
		 * constructed, but that it is not already initialized.
		 * We also need to define the Element over the field.
		 * So what we always initialize the Element with the zero field value.
		 * If an integer different from zero is passed to the function the Element
		 * is initialized to a constant polynom of Z/pZ
		 * @return reference to field Element.
		 * @param x field Element to contain output (reference returned).
		 * @param y integer.
		 */
		Element& init(Element& x , const integer& y=0) const
		{
			x.assign_zero(*this);
			if (y==1)
				x.assign_one(*this);
			else
				{
					if (y!=0)
						if (this->degree() > 1)						
							{
								LiDIA::Fp_polynomial Pol;
								LiDIA::bigint p=static_cast<const LiDIA::galois_field&>(*this).characteristic();
								Pol.set_modulus(p);
								//Pol.set_max_degree((x.get_field()).degree());
								
								//integer rem, quo,tmp=y;
								LiDIA::bigint rem, quo ,tmp;
								LiDIA::string_to_bigint(std::string(y).data(),tmp);
								LiDIA::lidia_size_t deg = this->degree();
								for(LiDIA::lidia_size_t i=0;i<deg;i++)
								 	{		
										quo=tmp/p;
										rem=tmp%p;
										tmp=quo;
										if (rem != LiDIA::bigint(0))
											Pol.set_coefficient(rem,i);
									}
								//Element e(x.get_field());	   
								x.set_polynomial_rep(Pol);
								//x.assign(e);
							}
						else {
							LiDIA::bigint tmp;
							LiDIA::string_to_bigint(std::string(y).data(),tmp);
							x.assign(tmp);
						}
				}
			return x;				
		}
    
    

		/** Conversion of field base Element to an integer.
		 * This function assumes the output field base Element x has already been
		 * constructed, but that it is not already initialized.
		 * As Elements are represented by polynom the convert function return 
		 * the valuation of polynom in characteristic by the Horner Method.
		 * That keeps unicity of each Element.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base Element.
		 */
		integer& convert(integer& x , const Element& y ) const
		{
			LiDIA::bigint fx(0) , X((y.get_field()).characteristic());
			LiDIA::bigint tmp;
	
	 
			for(int i=(int)(y.get_field()).degree();i>0;i--)
				{
					(y.polynomial_rep()).get_coefficient(tmp,i);
					fx=fx+tmp;
					fx=fx*X;
				}
			(y.polynomial_rep()).get_coefficient(tmp,0);
			fx= fx + tmp;

			//long i;
			//fx.longify(i);
			unsigned int size= fx.bit_length();

			char s[(size>>2)+4];
			LiDIA::bigint_to_string(fx,s);

			return x=integer(s);
		}




		/** Assignment of one field Element to another.
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& assign(Element& x, const Element& y) const
		{
			x.assign(y);
			return x;
		}




		/** Cardinality.
		 * Return integer representing cardinality of the field.
		 * Returns  p^k.
		 * @return constant reference to integer representing cardinality 
		 *	       of the field.
		 */
		integer& cardinality(integer& c) const 
		{
			LiDIA::bigint tmp=number_of_elements();
			unsigned int size= tmp.bit_length();
			char s[(size>>2)+4];
			LiDIA::bigint_to_string(tmp,s);

			return c=integer(s);
 
		}



		/** Characteristic.
		 * Return integer representing characteristic of the field.
		 * Returns p.
		 * @return constant reference to integer representing characteristic 
		 * 	       of the field.
		 */
		integer& characteristic(integer& c) const
		{
			
			LiDIA::bigint tmp = static_cast<const LiDIA::galois_field&>(*this).characteristic();
			unsigned int size= tmp.bit_length();
			char s[(size>>2)+4];
			LiDIA::bigint_to_string(tmp,s);

			return c= integer(s);
		}


		//@} Object Management
    
		/** @name Arithmetic Operations 
		 * x <- y op z; x <- op y
		 * These operations require all Elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field Elements will
		 * give undefined results.
		 */
		//@{
     
		/** Equality of two Elements.
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field Element
		 * @param  y field Element
		 */
		bool areEqual(const Element& x, const Element& y) const
		{
			return x==y;
		}


		/** Addition.
		 * x = y + z
		 * This function assumes all the field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 * @param  z field Element.
		 */
		Element& add(Element& x, const Element& y, const Element& z) const
		{ 
			LiDIA::add(x,y,z);
			return x;
		}
    
		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 * @param  z field Element.
		 */
		Element& sub(Element& x, const Element& y, const Element& z) const
		{
			LiDIA::subtract(x,y,z);
			return x;
		}
     
		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 * @param  z field Element.
		 */
		Element& mul(Element& x, const Element& y, const Element& z) const
		{
			LiDIA::multiply(x,y,z);
			return x; 
		}
     
		/** Division.
		 * x = y / z
		 * This function assumes all the field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 * @param  z field Element.
		 */
		Element& div(Element& x, const Element& y, const Element& z) const
		{
			LiDIA::divide(x,y,z);
			return x;
		}
     
     
		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& neg(Element& x, const Element& y) const
		{
			LiDIA::negate(x,y);
			return x;
		}



		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& inv(Element& x, const Element& y) const
		{
			LiDIA::invert(x,y);
			return x;
		}
     
     

		/** Natural AXPY.
		 * r  = a * x + y
		 * This function assumes all field Elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field Element (reference returned).
		 * @param  a field Element.
		 * @param  x field Element.
		 * @param  y field Element.
		 */
		Element& axpy(Element& r, 
			      const Element& a,
			      const Element& x, 
			      const Element& y) const
		{
			return r=a*x+y;
		}

      

		/** Zero equality.
		 * Test if field Element is equal to zero of field.
		 * This function assumes the field Element has already been 
		 * constructed and initialized.
		 * @return boolean true if equals zero of field, false if not.
		 * @param  x field Element.
		 */
		bool isZero(const Element& x) const 
		{
			return x.is_zero();
		}



		/** One equality.
		 * Test if field Element is equal to one of field.
		 * This function assumes the field Element has already been 
		 * constructed and initialized.
		 * @return boolean true if equals one of field, false if not.
		 * @param  x field Element.
		 */
		bool isOne(const Element& x) const 
		{
			return x.is_one();
		}


     
		/** Inplace Addition.
		 * x += y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& addin(Element& x, const Element& y) const {
			return x+=y; 
		}
     
		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& subin(Element& x, const Element& y) const {
			return x-=y; 
		}
     
		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& mulin(Element& x, const Element& y) const {
			return x*=y; 
		}
    
		/** Inplace Division.
		 * x /= y
		 * This function assumes both field Elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 * @param  y field Element.
		 */
		Element& divin(Element& x, const Element& y) const { 
			return x/=y;
		}
     
     
		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field Element has already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 */
		Element& negin(Element& x) const
		{
			x.negate();
			return x;
		}
     
     

		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field Elementhas already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field Element (reference returned).
		 */
		Element& invin(Element& x) const
		{
			x.invert();
			return x;
		}



		/** Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field Elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field Element (reference returned).
		 * @param  a field Element.
		 * @param  x field Element.
		 */
		Element& axpyin(Element& r, const Element& a, const Element& x) const
		{			 
			return  r+=a*x;
		}

		//@} Inplace Arithmetic Operations

#ifndef __LINBOX_XMLENABLED
		/** @name Input/Output Operations */
		//@{
    
		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream& write(std::ostream& os) const
		{
			integer c;
			return os<<"corps de Galois GF("<<
				characteristic(c)<<"^"<<degree()<<")";
		}
     
     
		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream& read(std::istream& is) const
		{
			return is ;
		}
     

		/** Print field Element like a polynom.
		 * @return output stream to which field Element is written.
		 * @param  os  output stream to which field Element is written.
		 * @param  x   field Element.
		 */
		std::ostream& write(std::ostream& os,const Element& e) const
		{
			/*integer tmp; 
			  (*this).convert(tmp,e);
			  return os<<tmp;
			*/
			return os<<e;
		}
     
     
     
		/** Read field Element.
		 * @return input stream from which field Element is read.
		 * @param  is  input stream from which field Element is read.
		 * @param  x   field Element.
		 */
		std::istream& read(std::istream& is, Element& e) const
		{ 
			// return is>>e;
			integer tmp;
			is>> tmp;
			(*this).init(e,tmp);
			return is;
	 
		}

		//@}
      
#else
		std::ostream &write(std::ostream &os) const
		{
			LinBox::Writer W;
			if( toTag(W)) 
				W.write(os);
			
			return os;
		}

		bool toTag(LinBox::Writer &W) const
		{
			LiDIA::bigint card, charac;
			LiDIA::lidia_size_t deg;
			std::string s;

			W.setTagName("field");
			W.setAttribute("implDetail", "lidia-gfq");

			card = LiDIA::galois_field::number_of_elements();
			W.setAttribute("cardinality", LinBox::Writer::numToString(s, card));
			W.addTagChild();
			W.setTagName("finite");

			W.addTagChild();
			W.setTagName("characteristic");
			charac = LiDIA::galois_field::characteristic();
			W.addNum(charac);
			W.upToParent();

			W.addTagChild();
			W.setTagName("extension");
			deg = LiDIA::galois_field::degree();
			W.addNum(deg);
			W.upToParent();

			W.upToParent();

			return true;
		}

		std::ostream &write(std::ostream &os, const Element &e) const 
		{
			LinBox::Writer W;
			if (toTag(W, e))
				W.write(os);

			return os;
		}

		bool toTag(LinBox::Writer &W, const Element &e) const
		{
			std::string s;
			LiDIA::Fp_polynomial poly = e.polynomial_rep();
			LiDIA::bigint accum = 0, base = LiDIA::galois_field::characteristic();
			LiDIA::lidia_size_t i;

			for(i = poly.degree(); i >= 0; i--) {
				accum *= base;
				accum += poly[i];
			}

			W.setTagName("cn");
			W.addDataChild(LinBox::Writer::numToString(s, accum));

			return true;
		}

		std::istream &read(std::istream &is, Element &e) const
		{
			LinBox::Reader R(is);
			if( !fromTag(R, e)) {
				is.setstate(std::istream::failbit);
				if(!R.initalized())
					is.setstate(std::istream::badbit);
			}

			return is;
		}

		bool fromTag(LinBox::Reader &R, Element &e) const
		{
			LiDIA::bigint total, base = LiDIA::galois_field::characteristic();
			LiDIA::lidia_size_t deg = LiDIA::galois_field::degree(), i = 0;

			if(!R.expectTagName("cn") || !R.expectTextNum(total)) return false;
			LiDIA::Fp_polynomial f;
			f.set_modulus(base);
			f.set_max_degree(deg); // ensure the poly is correct size

			// initalize the polynomial
			while(total > 0) {
				f[i] = total % base;
				total /= base;
				++i;
			}


			e.set_polynomial_rep(f);

			return true;
		}
#endif			
			
			
		static inline integer getMaxModulus()
		{ return integer( "9007199254740881" ); } // prevprime(2^53)

	}; // class lidia-gfq

	integer& FieldTraits< LidiaGfq >::maxExponent( integer& i )
	{ return i = integer(2147483647); } // 2^31-1

} // namespace LinBox


#endif // __LINBOX_lidia_gfq_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
