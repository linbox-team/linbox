/* Copyright (C) 2010 LinBox
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


#ifndef __LINBOX_pir_modular_H
#define __LINBOX_pir_modular_H

#include <givaro/modular.h>

#include "linbox/field/field-traits.h"

// Namespace in which all LinBox code resides
namespace LinBox
{


	template <intType>
	struct ClassifyRing<PIRModular<intType> >  {
		typedef RingCategories::ModularTag categoryTag;
	};

	/// \ingroup ring
	template <intType>
	class PIRModular<intType> : public Givaro::Modular<intType> {

	public:

		using Parent_t = Givaro::Modular<intType>;
	
		using Givaro::Modular<intType>::Element;
		//typedef typename Givaro::Modular<intType>::Element Element;

		using Givaro::Modular<intType>::RandIter;

		//No default cstor

		PIRModular (intType value, uint32_t exp = 1) :
			Givaro::Modular<intType>(Givaro::power(value, exp)), 
			_irred(value), 
			_exponent(exp)
		{}

        using Parent_t:: zero;
        using Parent_t:: one;
        using Parent_t:: mOne;

        using Parent_t:: minElement;
        using Parent_t:: maxElement;

        using Parent_t:: characteristic;
        using Parent_t:: cardinality;
        
        using Parent_t:: maxCardinality;
        using Parent_t:: minCardinality;

        // ----- Checkers
        using Parent_t:: isZero;
        using Parent_t:: isOne ;
        using Parent_t:: isMOne;
        using Parent_t:: areEqual;
        using Parent_t:: length;
        
        // ----- Ring-wise operators
        using Parent_t:: operator==;
        using Parent_t:: operator!=;
        using Parent_t:: operator=;

        // ----- Initialisation
        using Parent_t:: init ;

        using Parent_t:: assign ;
    
        using Parent_t:: convert;

        using Parent_t:: reduce ;
        
        using Parent_t:: mul;
        using Parent_t:: div;
        using Parent_t:: add;
        using Parent_t:: sub;
        using Parent_t:: neg;
        using Parent_t:: inv;

        using Parent_t:: mulin;
        using Parent_t:: divin;
        using Parent_t:: addin;
        using Parent_t:: subin;
        using Parent_t:: negin;
        using Parent_t:: invin;
        
        using Parent_t:: axpy  ;
        using Parent_t:: axpyin;

        using Parent_t:: axmy;
        using Parent_t:: axmyin;

        using Parent_t:: maxpy;
        using Parent_t:: maxpyin;

        // ----- Random generators
        using Parent_t:: NonZeroRandIter;
        using Parent_t:: random;
        using Parent_t:: nonzerorandom;

        // --- IO methods
        std::istream& read (std::istream& s);
		{ string s; return Parent_t::read(s >> str)>> s >> _irred >> s >> _exponent; }
        std::ostream& write(std::ostream& s) const
		{ return Parent_t::write(s<<"Local- ") << "irred: " << _irred << ", exponent: " << _exponent; }
        std::istream& read ;
        std::ostream& write;
        
	protected:
		intType _irred;
		uint32_t _exponent;
	};


	template <>
	class DotProductDomain<PIRModular<int32_t> > : public  VectorDomainBase<PIRModular<int32_t> > {

	public:
		typedef int32_t Element;
		DotProductDomain (const PIRModular<int32_t> &F) :
			VectorDomainBase<PIRModular<int32_t> > (F)
		{}


	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;

			uint64_t y = 0;
			uint64_t t;

			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
				t = ( (uint64_t) *i ) * ( (uint64_t) *j );
				y += t;

				if (y < t)
					y += field()._two_64;
			}

			y %= (uint64_t) field().characteristic();
			return res = Element(y);

		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;

			uint64_t y = 0;
			uint64_t t;

			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
				t = ( (uint64_t) *i_elt ) * ( (uint64_t) v2[*i_idx] );
				y += t;

				if (y < t)
					y += field()._two_64;
			}


			y %= (uint64_t) field().characteristic();

			return res = (Element)y;
		}
	};


	// Specialization of MVProductDomain for int32_t modular field
	template <>
	class MVProductDomain<PIRModular<int32_t> > {
	public:

		typedef int32_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIRModular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64_t) *k) * ((uint64_t) *j);

				*l += t;

				if (*l < t)
                                    *l += VD.faxpy().field()._two_64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field().characteristic();

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIRModular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
	{
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64_t) k->second) * ((uint64_t) *j);

				_tmp[k->first] += t;

				if (_tmp[k->first] < t)
					_tmp[k->first] += VD.faxpy().field()._two_64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field().characteristic();

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIRModular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64_t) k->second) * ((uint64_t) *j);

				_tmp[k->first] += t;

				if (_tmp[k->first] < t)
					_tmp[k->first] += VD.faxpy().field()._two_64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field().characteristic();

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<PIRModular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<PIRModular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::first_type::const_iterator k_idx;
		typename Matrix::Column::second_type::const_iterator k_elt;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
			{
				t = ((uint64_t) *k_elt) * ((uint64_t) *j);

				_tmp[*k_idx] += t;

				if (_tmp[*k_idx] < t)
					_tmp[*k_idx] += VD.faxpy().field()._two_64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field().characteristic();

		return w;
	}


}


#endif //__LINBOX_pir_modular_int32_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
