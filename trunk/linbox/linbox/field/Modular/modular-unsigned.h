/* linbox/field/modular.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) 2011 LinBox
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * LargeGivaro::Modular is now replace by a class Givaro::Modular parameterized on the element
 * type. So, the old LargeGivaro::Modular is equivalent to Givaro::Modular<integer>. All other
 * interface details are exactly the same.
 *
 * Renamed from large-modular.h to modular.h
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_field_modular_unsigned_H
#define __LINBOX_field_modular_unsigned_H

namespace LinBox { /*  uint8_t */

	/*! Specialization of FieldAXPY for uint8_t modular field */

	template <>
	class FieldAXPY<Givaro::Modular<uint8_t> > {
	public:

		typedef uint8_t Element;
		typedef uint64_t Abnormal;
		typedef Givaro::Modular<uint8_t> Field;

		FieldAXPY (const Field &F) :
			_k (((uint64_t) -1LL) / ((F.characteristic() - 1) * (F.characteristic() - 1))),
			_field (&F),
			_y (0),
			i (_k)
		{
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_k (faxpy._k),
			_field (faxpy._field),
			_y (0),
			i (_k)
		{}

		FieldAXPY<Givaro::Modular<uint8_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			_k = faxpy._k;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint32_t t = (uint32_t) a * (uint32_t) x;

			if (!i--) {
				i = int(_k);
				return _y = _y % (uint32_t) field().characteristic() + t;
			}
			else
				return _y += t;
		}

		inline uint64_t& accumulate (const Element &t)
		{

			if (!i--) {
				i = int(_k);
				return _y = _y % (uint32_t) field().characteristic() + t;
			}
			else
				return _y += t;
		}

		inline Element &get (Element &y) const
		{
			const_cast<FieldAXPY<Field>*>(this)->_y %= (uint32_t) field().characteristic();
			if ((int32_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
			y = (uint8_t) _y;
			const_cast<FieldAXPY<Field>*>(this)->i = int(_k);
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			i = int(_k);
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

		inline const Field & field() const { return *_field; }
		
	public:
		uint64_t _k; 
		
	private:
		const Field *_field;
		uint64_t _y;
		int i;
	};

	//! Specialization of DotProductDomain for unsigned short modular field

	template <>
	class DotProductDomain<Givaro::Modular<uint8_t> > : public  VectorDomainBase<Givaro::Modular<uint8_t> > {
	public:

		typedef uint8_t Element;

		DotProductDomain(){}
		DotProductDomain (const Givaro::Modular<uint8_t> &F) :
			VectorDomainBase<Givaro::Modular<uint8_t> > (F)
		{}

		using VectorDomainBase<Givaro::Modular<uint8_t> >::field;
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        
	};

	//! Specialization of MVProductDomain for uint8_t modular field

	template <>
	class MVProductDomain<Givaro::Modular<uint8_t> > {
	public:

		typedef uint8_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Givaro::Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint32_t> _tmp;
	};

}

namespace LinBox { /*  uint16_t */

	/*! Specialization of FieldAXPY for uint16_t modular field */
	template <>
	class FieldAXPY<Givaro::Modular<uint16_t> > {
	public:

		typedef uint16_t Element;
		typedef Givaro::Modular<uint16_t> Field;

		FieldAXPY (const Field &F) :
			_k (((uint64_t) -1LL) / ((F.characteristic() - 1) * (F.characteristic() - 1))),
			_field (&F),
			_y (0),
			i (_k)
		{}
		
		FieldAXPY (const FieldAXPY &faxpy) :
			_k (faxpy._k), _field (faxpy._field), _y (0), i (_k)
		{}

		FieldAXPY<Givaro::Modular<uint16_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			_k = faxpy._k;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) ((long long) a * (long long) x);

			if (!i--) {
				i = (int)_k;
				return _y = _y % (uint64_t) field().characteristic() + t;
			}
			else
				return _y += t;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			if (!i--) {
				i = (int)_k;
				return _y = _y % (uint64_t) field().characteristic() + t;
			}
			else
				return _y += t;
		}

		inline Element &get (Element &y) const
		{
			const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
			if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
			y = (uint16_t) _y;
			const_cast<FieldAXPY<Field>*>(this)->i = int(_k);
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			i = (int)_k;
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

		inline const Field & field() const {return *_field;}
		
	public:
		uint64_t _k;
		
	private:
		const Field *_field;
		uint64_t _y;
		int i;
	};

	//! Specialization of DotProductDomain for unsigned short modular field

	template <>
	class DotProductDomain<Givaro::Modular<uint16_t> > : public VectorDomainBase<Givaro::Modular<uint16_t> > {
	public:

		typedef uint16_t Element;

		DotProductDomain () {}
		DotProductDomain (const Givaro::Modular<uint16_t> &F) :
			VectorDomainBase<Givaro::Modular<uint16_t> > (F)
		{}
		using VectorDomainBase<Givaro::Modular<uint16_t> >::field;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        
	};

	//! Specialization of MVProductDomain for uint16_t modular field

	template <>
	class MVProductDomain<Givaro::Modular<uint16_t> > {
	public:

		typedef uint16_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Givaro::Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

}

#include <givaro/modular-uint32.h>

namespace LinBox { /*  uint32_t */

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Givaro::Modular<Element> >;

	/*! Specialization of FieldAXPY for unsigned short modular field */

	template <>
	class FieldAXPY<Givaro::Modular<uint32_t> > {
	public:

		typedef uint32_t Element;
		typedef Givaro::Modular<uint32_t> Field;

		FieldAXPY (const Field &F) :
			_field (&F), _y(0)
		{
			_two_64 = (uint64_t(1) << 32) % uint64_t(F.characteristic());
			_two_64 = (_two_64 * _two_64) % uint64_t(F.characteristic());
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_two_64 (faxpy._two_64), _field (faxpy._field), _y (0)
		{}

		FieldAXPY<Givaro::Modular<uint32_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			_two_64 = faxpy._two_64;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;

			if (_y < t)
				return _y += _two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			_y += t;

			if (_y < t)
				return _y += _two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate_special (const Element &t)
		{
			return _y += t;
		}

		inline Element &get (Element &y) const
		{
			const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
			//if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
			return y = (uint32_t) _y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

		inline const Field & field() const { return *_field; }
		
	public:
	
		uint64_t _two_64;
		
	private:

		const Field *_field;
		uint64_t _y;
	};

	//! Specialization of DotProductDomain for uint32_t modular field

	template <>
	class DotProductDomain<Givaro::Modular<uint32_t> > : public VectorDomainBase<Givaro::Modular<uint32_t> > {
	public:

		typedef uint32_t Element;

		DotProductDomain () {}
		DotProductDomain (const Givaro::Modular<uint32_t> &F) :
			VectorDomainBase<Givaro::Modular<uint32_t> > (F)
		{}
		using VectorDomainBase<Givaro::Modular<uint32_t> >::field;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        
	};

	//! Specialization of MVProductDomain for uint32_t modular field

	template <>
	class MVProductDomain<Givaro::Modular<uint32_t> > {
	public:

		typedef uint32_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Givaro::Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

}

#include <givaro/modular-uint64.h>

namespace LinBox { /*  uint64_t */

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Givaro::Modular<Element> >;

	/*! Specialization of FieldAXPY for unsigned short modular field */

	template <>
	class FieldAXPY<Givaro::Modular<uint64_t> > {
	public:

		typedef uint64_t Element;
		typedef Givaro::Modular<uint64_t> Field;

		FieldAXPY (const Field &F) :
			_field (&F), _y(0)
		{
			_two_64 = (uint64_t(1) << 32) % uint64_t(F.characteristic());
			_two_64 = (_two_64 * _two_64) % uint64_t(F.characteristic());
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_two_64 (faxpy._two_64), _field (faxpy._field), _y (0)
		{}

		FieldAXPY<Givaro::Modular<uint64_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;

			if (_y < t)
				return _y += _two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			_y += t;

			if (_y < t)
				return _y += _two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate_special (const Element &t)
		{
			return _y += t;
		}

		inline Element &get (Element &y) const
		{
			const_cast<FieldAXPY<Field>*>(this)->_y %= (uint64_t) field().characteristic();
			//if ((int64_t) _y < 0) const_cast<FieldAXPY<Field>*>(this)->_y += field().characteristic();
			return y = (uint64_t) _y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

		inline const Field & field() const { return *_field; }
		
	public:
	
		uint64_t _two_64;
		
	private:

		const Field *_field;
		uint64_t _y;
	};

	//! Specialization of DotProductDomain for uint64_t modular field

	template <>
	class DotProductDomain<Givaro::Modular<uint64_t> > : public VectorDomainBase<Givaro::Modular<uint64_t> > {
	public:

		typedef uint64_t Element;

		DotProductDomain () {}
		DotProductDomain (const Givaro::Modular<uint64_t> &F) :
			VectorDomainBase<Givaro::Modular<uint64_t> > (F)
		{}
		using VectorDomainBase<Givaro::Modular<uint64_t> >::field;

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
        
	};

	//! Specialization of MVProductDomain for uint64_t modular field

	template <>
	class MVProductDomain<Givaro::Modular<uint64_t> > {
	public:

		typedef uint64_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Givaro::Modular<uint64_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint64_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint64_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint64_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Givaro::Modular<uint64_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

}

#include "linbox/field/Modular/modular.inl"

#endif // __LINBOX_field_modular_unsigned_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
