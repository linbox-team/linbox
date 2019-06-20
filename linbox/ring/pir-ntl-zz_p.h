/* Copyright (C) 2010 LinBox
 * written by Zhendong Wan
 *
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_pir_ntl_zz_p_H
#define __LINBOX_pir_ntl_zz_p_H

#include <givaro/zring.h>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include "linbox/ring/ntl/ntl-zz_p.h"
#include "linbox/vector/vector-domain.h"
#include <sstream>
#include "linbox/integer.h"
#include "linbox/field/field-traits.h"

namespace LinBox
{

	template<class Field>
	class FieldAXPY;

	template <class Ring>
	struct ClassifyRing;

	class PIR_ntl_ZZ_p;

	template <>
	struct ClassifyRing<PIR_ntl_ZZ_p> {
		typedef RingCategories::ModularTag categoryTag;
	};

	/** \brief extend Wrapper of ZZ_p from NTL.  Add PIR functions
	  \ingroup field
	  */

	class PIR_ntl_ZZ_p :  public NTL_ZZ_p {
		//public Givaro::ZRing<NTL::ZZ_p> {

		// typedef Givaro::ZRing<NTL::ZZ_p> Father_t ;
		typedef NTL_ZZ_p Father_t ;
	public:
		typedef Father_t::Element Element;

		template <class Int_t>
		PIR_ntl_ZZ_p(const Int_t& d) :
			Father_t(d)
		{
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		}

		PIR_ntl_ZZ_p (const NTL::ZZ& d) :
		Father_t(d)
		{
			// NTL::ZZ_p::init(d);
		}

		PIR_ntl_ZZ_p (const integer& d, int e = 1 )  :
			Father_t(d,(size_t)e)
		{

			// if(e != 1) throw PreconditionFailed(LB_FILE_LOC,"exponent must be 1");
			linbox_check(e == 1);

			// NTL::ZZ_p::init (NTL::to_ZZ(((std::string)d). c_str()));

		}

		inline integer& cardinality (integer& c) const {

			std::stringstream io;

			io << NTL::ZZ_p::modulus();

			io >> c;

			return c;
		}

		inline NTL::ZZ& cardinality (NTL::ZZ& c) const {

			return c = NTL::ZZ_p::modulus();
		}

		inline integer& characteristic (integer& c) const {

			std::stringstream  io;

			io << NTL::ZZ_p::modulus();

			io >> c;

			return c;
		}

		std::ostream& write (std::ostream& out) const {
			return out << "PIR_NTL_ZZ_p Ring";
		}

		std::istream& read (std::istream& in)  {
			return in;
		}


        inline Element& init (Element& x) const {return x;}
            
		/** @brief
		 *  Init x from y.
		 */
        
		template<class Element2>
		inline Element& init (Element& x,  const Element2& y) const {

			NTL::conv (x, y);

			return x;
		}

		/** @brief
		 *   Init from a NTL::ZZ_p
		 */
		inline Element& init (Element& x, const Element& y) const {

			x = y;

			return x;
		}

		/** @brief
		 *  I don't  know how to init from integer.
		 */
		inline Element& init (Element& x, const integer& y) const {


			NTL::conv(x, NTL::to_ZZ( (static_cast<const std::string>(y)).c_str() ) );
			return x;
		}

		/** @brief
		 *  Convert y to an Element.
		 */
		integer& convert (integer& x, const Element& y) const {
			bool neg=false;
			if (NTL::sign(NTL::rep(y)) <0)
				neg=true;
			long b = NTL::NumBytes(NTL::rep(y));
			unsigned char* byteArray; byteArray = new unsigned char[(size_t)b ];
			NTL::BytesFromZZ(byteArray, NTL::rep(y), b);
			integer base(256);
			x= integer(0);
			for(long i = b - 1; i >= 0; --i) {
				x *= base;
				x += integer(byteArray[i]);
			}
			delete [] byteArray;
			if (neg)
				x=-x;
			return x;
		}


		/** @brief
		 *  x = y.
		 */
		inline Element&  assign (Element& x, const Element& y) const {
			return x = y;
		}

		/** @brief
		 *  Test if x == y
		 */
		inline bool areEqual (const Element& x ,const Element& y) const {
			return x == y;
		}

		/** @brief
		 *  Test if x == 0
		 */
		inline bool isZero (const Element& x) const {

			return NTL::IsZero (x);
		}

		/** @brief
		 *  Test if x == 1
		 */
		inline bool isOne (const Element& x) const {
			return NTL::IsOne (x);
		}

        inline bool isMOne (const Element& x) const
		{
			Element y ; neg(y,x);
			return isOne(y);
		}

		// arithmetic

		/** @brief
		 *  return x = y + z
		 */
		inline Element& add (Element& x, const Element& y, const Element& z) const {

			NTL::add (x, y, z);

			return x;
		}

		/** @brief
		 *  return x = y - z
		 */
		inline Element& sub (Element& x, const Element& y, const Element& z) const  {

			NTL::sub (x, y, z);

			return x;
		}

		/** @brief
		 *  return x = y * z
		 */
		template <class Int>
		inline Element& mul (Element& x, const Element& y, const Int& z) const {

			NTL::mul (x, y, z);

			return x;
		}

		/** @brief
		 *  If exists a, such that a * z =y,
		 *  return x = one of them.
		 *  Otherwise, throw an exception
		 */
		inline Element& div (Element& x, const Element& y, const Element& z) const {

			NTL::ZZ g, s, t;

			NTL::XGCD (g, s, t, NTL::rep(z), NTL::ZZ_p::modulus());

			NTL::ZZ q, r;

			NTL::DivRem (q, r, NTL::rep(y), g);

			if (NTL::IsZero (r)) {

				Element tmp1, tmp2;

				NTL::conv (tmp1, s);

				NTL::conv (tmp2, q);

				NTL::mul (x, tmp1, tmp2);

			}

			else
				throw PreconditionFailed(LB_FILE_LOC,"Div: not dividable");


			return x;

		}

		/** @brief
		 *  If y is a unit, return x = 1 / y,
		 *  otherwsie, throw an exception
		 */
		inline Element& inv (Element& x, const Element& y) const {

			NTL::inv (x, y);

			return x;
		}

		/** @brief
		 *  return x = -y;
		 */
		inline Element& neg (Element& x, const Element& y) const {

			NTL::negate (x, y);

			return x;
		}


		/** @brief
		 *  return r = a x + y
		 */

		template <class Int>
		inline const Element& axpy (Element& r, const Element& a, const Int& x, const Element& y) const  {

			NTL::mul (r, a, x);

			return r += y;
		}


		// inplace operator

		/** @brief
		 *  return x += y;
		 */
		inline Element& addin (Element& x, const Element& y) const {

			return x += y;
		}

		/** @brief
		 *  return x -= y;
		 */
		inline Element& subin (Element& x, const Element& y) const {

			return x -= y;
		}

		/** @brief
		 *  return x *= y;
		 */
		template<class Int>
		inline Element& mulin (Element& x, const Int& y) const {

			return x *= y;
		}

		/** @brief
		 *  If y divides x, return x /= y,
		 *  otherwise throw an exception
		 */
		inline const Element& divin (Element& x, const Element& y) const {

			div (x, x, y);

			return x;
		}

		/** @brief
		 *  If x is a unit, x = 1 / x,
		 *  otherwise, throw an exception.
		 */
		inline Element& invin (Element& x) const {

			return x = NTL::inv(x);
		}

		/** @brief
		 *  return x = -x;
		 */
		inline Element& negin (Element& x) const {

			NTL::negate (x, x);

			return x;
		}

		/** @brief
		 *  return r += a x
		 */
		template <class Int>
		inline Element& axpyin (Element& r, const Element& a, const Int& x) const {

			return r += a * x;
		}


		// IO

		/** @brief
		 *  out << y;
		 */
		std::ostream& write(std::ostream& out,const Element& y) const {

			out << y;

			return out;
		}


		/** @brief
		 *  read x from istream in
		static std::istream& read(std::istream& in, Element& x) {

			return in >> x;
		}
		 */
		using Father_t::read;


		/** some PIR function
		*/

		/** @brief
		 *  Test if x is a unit.
		 */
		inline bool isUnit(const Element& x) const {
            NTL::ZZ d;
            return !NTL::InvModStatus(d,rep(x),NTL::ZZ_p::modulus());
        }


		/** @brief
		 *  return g = gcd (a, b)
		 */
		inline static Element& gcd (Element& g, const Element& a, const Element& b) {

			NTL::ZZ d;

			NTL::GCD (d, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, d);

			return g;
		}

		/** @brief
		 *  return g = gcd (g, b)
		 */
		inline static Element& gcdin (Element& g, const Element& b) {

			gcd (g, g, b);

			return g;
		}

		/** @brief
		 *  g = gcd(a, b) = a*s + b*t.
		 *  and gcd (s, t) is a unit.
		 */
		inline static Element& xgcd (Element& g, Element& s, Element& t, const Element& a, const Element& b){

			NTL::ZZ g1, s1, t1;

			NTL::XGCD (g1, s1, t1, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, g1);

			NTL::conv (s, s1);

			NTL::conv (t, t1);

			return g;
		}

		/** @brief
		 *  g = gcd(a, b) = a*s + b*t.
		 *  and gcd (s, t) is a unit.
		 *  s * a1 + t * b1 = a unit.
		 */
		inline static Element& dxgcd (Element& g, Element& s, Element& t, Element& a1, Element& b1,
					      const Element& a, const Element& b){

			NTL::ZZ g1, s1, t1, a2, b2;

			NTL::XGCD (g1, s1, t1, NTL::rep(a), NTL::rep(b));

			NTL::conv (g, g1);

			NTL::conv (s, s1);

			NTL::conv (t, t1);

			if (NTL::IsZero (g1)) {

				a1 = s;

				b1 = t;
			}

			else {

				NTL::div (a2, NTL::rep(a), g1);

				NTL::div (b2, NTL::rep(b), g1);

				NTL::conv (a1, a2);

				NTL::conv (b1, b2);
			}

			return g;
		}

		/** @brief
		 *  Test if a | b.
		 */
		inline static bool isDivisor (const Element& a, const Element& b) {

			if ( NTL::IsZero (a) ) return false;

			else if (NTL::IsZero (b)) return true;

			else {
				NTL::ZZ g, r;

				NTL::GCD (g, NTL::rep(a), NTL::ZZ_p::modulus());

				NTL::rem (r, NTL::rep(b), g);

				return NTL::IsZero (r);
			}
		}

		/** @brief
		 *  a = normalization of b.
		 */

		inline static Element& normal (Element& a,  const Element& b) {

			NTL::ZZ a1;

			NTL::GCD (a1, NTL::rep(b), NTL::ZZ_p::modulus());

			NTL::conv (a, a1);

			return a;
		}

		/** @brief
		*/

		inline static Element& normalIn (Element& a) {


			NTL::ZZ a1;

			NTL::GCD (a1, NTL::rep(a), NTL::ZZ_p::modulus());

			NTL::conv (a, a1);

			return a;
		}

		inline static integer maxCardinality()
		{ return integer( "4294967295" ); } // 2^32 - 1

	};



	template <>
	class FieldAXPY<PIR_ntl_ZZ_p>  {
	public:
		typedef PIR_ntl_ZZ_p Field;
		typedef Field::Element Element;
		typedef NTL::ZZ Abnormal;

		/** Constructor.
		 * A faxpy object if constructed from a Field.
		 * Copies of this objects are stored in the faxpy object.
		 * @param F field F in which arithmetic is done
		 */
		FieldAXPY (const Field &F) :
		       	_field (F)
		{ _y = NTL::ZZ::zero(); }

		/** Copy constructor.
		 * @param faxpy
		 */
		FieldAXPY (const FieldAXPY<Field> &faxpy) :
		       	_field (faxpy._field), _y (faxpy._y)
	       	{}

		/** Assignment operator
		 * @param faxpy
		 */
		FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
		{ _y = faxpy._y; return *this; }

		/** Add a*x to y
		 * y += a*x.
		 * @param a constant reference to element a
		 * @param x constant reference to element x
		 */
		inline NTL::ZZ& mulacc (const Element &a, const Element &x)
		{
			return _y +=  NTL::rep(a) * NTL::rep(x) ;
		}

		inline NTL::ZZ& accumulate (const Element &t)
		{
			return _y += NTL::rep(t);
		}

		/** Retrieve y
		 *
		 * Performs the delayed modding out if necessary
		 */
		inline Element &get (Element &y) { NTL::conv (y,  _y); return y; }

		/** Assign method.
		 * Stores new field element for arithmetic.
		 * @return reference to self
		 * @param y constant reference to element a
		 */
		inline FieldAXPY &assign (const Element& y)
		{
			_y = NTL::rep(y);

			return *this;
		}

		inline void reset() {
			_y = NTL::ZZ::zero();
		}

		inline const Field& field() const { return _field; }
	private:

		/// Field in which arithmetic is done
		/// Not sure why it must be mutable, but the compiler complains otherwise
		const Field &_field;

		/// Field element for arithmetic
		NTL::ZZ _y;

	};

	template<class Field>
	class DotProductDomain;

	template <>
	class DotProductDomain<PIR_ntl_ZZ_p> : public  VectorDomainBase<PIR_ntl_ZZ_p> {

	public:
		typedef PIR_ntl_ZZ_p::Element Element;
		DotProductDomain (const PIR_ntl_ZZ_p& F) :
			VectorDomainBase<PIR_ntl_ZZ_p> (F) {}


	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;

			NTL::ZZ y;
			//NTL::ZZ t;

			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {

				y += NTL::rep (*i) * NTL::rep(*j);

			}

			//NTL::rem (t, y, NTL::ZZ_p::modulus());

			NTL::conv (res, y);

			return res;

		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;

			NTL::ZZ y;
			y = (long)0;
			//NTL::ZZ t;

			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {

				y += NTL::rep(*i_elt) * NTL::rep(v2[*i_idx]);

			}


			//NTL::rem (t, y, NTL::ZZ_p::modulus());

			NTL::conv (res, y);

			return res;
		}
	};

	// Specialization of MVProductDomain for  PIR_ntl_ZZ_p field

	template <class Field>
	class MVProductDomain;

	template <>
	class MVProductDomain<PIR_ntl_ZZ_p>
	{
	public:

		typedef PIR_ntl_ZZ_p::Element Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;

		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

	};

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIR_ntl_ZZ_p>::mulColDenseSpecialized
	(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<NTL::ZZ>::iterator l;
		std::vector<NTL::ZZ> _tmp(w.size());

		//NTL::ZZ t;

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)

				*l += NTL::rep (*k) + NTL::rep (*j);

		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l) {

			//NTL::rem (t, *l, NTL::ZZ_p::modulus());

			NTL::conv (*w_j, *l);
		}


		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIR_ntl_ZZ_p>::mulColDenseSpecialized
	(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
	{
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<NTL::ZZ>::iterator l;
		std::vector<NTL::ZZ> _tmp(w.size());

		//NTL::ZZ t;

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {

				_tmp[k->first] += NTL::rep (k->second) * NTL::rep (*j);

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l) {

			//NTL::rem (t, *l, NTL::ZZ_p::modulus());

			NTL::conv (*w_j, *l);
		}

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIR_ntl_ZZ_p >::mulColDenseSpecialized
	(const VectorDomain<PIR_ntl_ZZ_p > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<NTL::ZZ>::iterator l;

		std::vector<NTL::ZZ> _tmp(w.size());

		//NTL::ZZ t;

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {

				_tmp[k->first] += NTL::rep(k -> second) * NTL::rep(*j);

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l) {

			//NTL::rem (t, *l, NTL::ZZ_p::modulus());

			NTL::conv (*w_j, *l);
		}


		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIR_ntl_ZZ_p>::mulColDenseSpecialized
	(const VectorDomain<PIR_ntl_ZZ_p> &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::first_type::const_iterator k_idx;
		typename Matrix::Column::second_type::const_iterator k_elt;
		std::vector<NTL::ZZ>::iterator l;

		std::vector<NTL::ZZ> _tmp(w.size());
		//NTL::ZZ t;

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)

				_tmp[*k_idx] += NTL::rep(*k_elt) * NTL::rep(*j);

		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l) {

			//NTL::rem (t, *l, NTL::ZZ_p::modulus());

			NTL::conv (*w_j, *l);
		}

		return w;
	}



}

#endif //__LINBOX_pir_ntl_zz_p_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
