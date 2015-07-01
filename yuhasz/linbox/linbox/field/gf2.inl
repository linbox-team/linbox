/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/gf2.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FIELD_GF2_INL
#define __FIELD_GF2_INL

#include <iostream>
#include <time.h>

#include "linbox/field/gf2.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/vector/stream.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/matrix/matrix-domain.h"

namespace LinBox 
{ 

// Specialization of canonical vector types

template <>
class RawVector<bool>
{
    public:
	typedef BitVector Dense;
	typedef std::vector<uint32> Sparse;
};

// Specialization of DotProductDomain for GF2

template <>
class DotProductDomain<GF2> : private virtual VectorDomainBase<GF2>
{
    public:

	typedef bool Element;

	DotProductDomain (const GF2 &F)
		: VectorDomainBase<GF2> (F)
	{}

    protected:
	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

	template <class Vector1, class Vector2>
	inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;

	template <class Vector1, class Vector2>
	inline BitVector::reference dotSpecializedDD (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const;

	template <class Vector1, class Vector2>
	inline BitVector::reference dotSpecializedDSP (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const;
};

// Specialization of vector domain

template <>
class VectorDomain<GF2> : private virtual VectorDomainBase<GF2>, private DotProductDomain<GF2>
{
    public:
	typedef bool Element;

	VectorDomain (const VectorDomain &VD)
		: DotProductDomain<GF2> (VD._F), VectorDomainBase<GF2> (VD._F)
	{}

	VectorDomain &operator = (const VectorDomain &VD) { return *this; }

	const GF2 &field () const { return _F; }
    
	template <class Vector>
	inline std::ostream &write (std::ostream &os, const Vector &x) const
		{ return writeSpecialized (os, x, VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector>
	inline std::istream &read (std::istream &is, Vector &x) const
		{ return readSpecialized (is, x, VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
		{ return copySpecialized (res, v,
					  VectorTraits<Vector1>::VectorCategory (),
					  VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
		{ return copySpecialized (res, v, i, len,
					  VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		{ return areEqualSpecialized (v1, v2,
					      VectorTraits<Vector1>::VectorCategory (),
					      VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector>
	inline bool isZero (const Vector &v) const
		{ return isZeroSpecialized (v, VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dotSpecialized (res, v1, v2,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline BitVector::reference dot (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const
		{ return dotSpecialized (res, v1, v2,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dot (res, v1, v2); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory (),
					 VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   VectorTraits<Vector1>::VectorCategory (),
					   VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 VectorTraits<Vector1>::VectorCategory (),
					 VectorTraits<Vector2>::VectorCategory (),
					 VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   VectorTraits<Vector1>::VectorCategory (),
					   VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		{ copy (res, x); return res; }

	template <class Vector>
	inline Vector &negin (Vector &y) const
		{ return y; }

	template <class Vector1, class Vector2>
	inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element a) const
		{ return mulSpecialized (res, x, a, VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector>
	inline Vector &mulin (Vector &x, const Element a) const
		{ return mulinSpecialized (x, a, VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &axpy (Vector1 &res, const Element a, const Vector2 &x, const Vector3 &y) const
		{ if (a) add (res, x, y); else this->copy (res, y); return res; }

	template <class Vector1, class Vector2>
	inline Vector1 &axpyin (Vector1 &y, const Element a, const Vector2 &x) const
		{ if (a) addin (y, x); return y; }

	VectorDomain (const GF2 &F)
		: DotProductDomain<GF2> (F), VectorDomainBase<GF2> (F)
	{}

    protected:

	// Specialized function implementations
	template <class Vector, class Trait>
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::DenseZeroOneVectorTag<Trait>) const;
	template <class Vector, class Trait>
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::SparseZeroOneVectorTag<Trait>) const;

	template <class Vector, class Trait>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::DenseZeroOneVectorTag<Trait>) const;
	template <class Vector, class Trait>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::SparseZeroOneVectorTag<Trait>) const;

	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag<Trait1>,
				  VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return v1 == v2; }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag<Trait1>,
				  VectorCategories::SparseZeroOneVectorTag<Trait2>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag<Trait1>,
					 VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return areEqual (v2, v1); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::SparseZeroOneVectorTag<Trait1>,
				  VectorCategories::SparseZeroOneVectorTag<Trait2>) const
		{ return v1 == v2; }

	template <class Vector, class Trait>
	bool isZeroSpecialized (const Vector &v, VectorCategories::DenseZeroOneVectorTag<Trait>) const;
	template <class Vector, class Trait>
	inline bool isZeroSpecialized (const Vector &v,
				       VectorCategories::SparseZeroOneVectorTag<Trait>) const
		{ return v.empty (); }

	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::DenseZeroOneVectorTag<Trait1>,
					 VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ std::copy (v.begin (), v.end (), res.begin ()); return res; }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::SparseZeroOneVectorTag<Trait1>,
				  VectorCategories::DenseZeroOneVectorTag<Trait2>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::DenseZeroOneVectorTag<Trait1>,
				  VectorCategories::SparseZeroOneVectorTag<Trait2>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::SparseZeroOneVectorTag<Trait1>,
					 VectorCategories::SparseZeroOneVectorTag<Trait2>) const
		{ res = v; return res; }

	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag<Trait1>,
					VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag<Trait1>,
					VectorCategories::SparseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::SparseZeroOneVectorTag<Trait1>,
					VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
				 VectorCategories::SparseZeroOneVectorTag<Trait1>,
				 VectorCategories::SparseZeroOneVectorTag<Trait2>) const;

	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag<Trait1>,
					VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag<Trait1>,
					VectorCategories::SparseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::SparseZeroOneVectorTag<Trait1>,
					VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1); }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
				 VectorCategories::SparseZeroOneVectorTag<Trait1>,
				 VectorCategories::SparseZeroOneVectorTag<Trait2>) const;

	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag<Trait1>,
				 VectorCategories::DenseZeroOneVectorTag<Trait2>,
				 VectorCategories::DenseZeroOneVectorTag<Trait3>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag<Trait1>,
				 VectorCategories::DenseZeroOneVectorTag<Trait2>,
				 VectorCategories::SparseZeroOneVectorTag<Trait3>) const
		{ copy (res, y); addin (res, x); }
	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::SparseZeroOneVectorTag<Trait1>,
				 VectorCategories::SparseZeroOneVectorTag<Trait2>,
				 VectorCategories::SparseZeroOneVectorTag<Trait3>) const;

	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag<Trait1>,
				   VectorCategories::DenseZeroOneVectorTag<Trait2>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag<Trait1>,
				   VectorCategories::SparseZeroOneVectorTag<Trait2>) const;
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag<Trait1>,
				   VectorCategories::DenseZeroOneVectorTag<Trait2>) const
		{ Vector1 xp, res; copy (xp, x); add (res, y, xp); copy (y, res); return y; }
	template <class Vector1, class Trait1, class Vector2, class Trait2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag<Trait1>,
				   VectorCategories::SparseZeroOneVectorTag<Trait2>) const
		{ Vector1 res; add (res, y, x); this->copy (y, res); return y; }

	template <class Vector1, class Vector2, class Trait>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::DenseZeroOneVectorTag<Trait> tag) const
		{ if (a) this->copy (res, x); else std::fill (res.wordBegin (), res.wordEnd (), 0); return res; }
	template <class Vector1, class Vector2, class Trait>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::SparseZeroOneVectorTag<Trait> tag) const
		{ if (a) this->copy (res, x); else res.clear (); return res; }

	template <class Vector, class Trait>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::DenseZeroOneVectorTag<Trait>) const
		{ if (!a) std::fill (x.wordBegin (), x.wordEnd (), 0); return x; }

	template <class Vector, class Trait>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::SparseZeroOneVectorTag<Trait> tag) const
		{ if (!a) x.clear (); return x; }

	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag<Trait1>,
					VectorCategories::GenericVectorTag<Trait2>,
					VectorCategories::GenericVectorTag<Trait3>) const
	{
		typename LinBox::Vector<GF2>::Sparse v;
		typename LinBox::Vector<GF2>::Sparse w;
		typename LinBox::Vector<GF2>::Sparse u;

		copy (v, x);
		copy (w, y);
		add (u, w, v);
		copy (res, u);

		return u;
	}

	template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
	inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag<Trait1>,
					VectorCategories::GenericVectorTag<Trait2>,
					VectorCategories::GenericVectorTag<Trait3>) const
	{
		typename LinBox::Vector<GF2>::Sparse v;
		typename LinBox::Vector<GF2>::Sparse w;
		typename LinBox::Vector<GF2>::Sparse u;

		copy (v, x);
		copy (w, y);
		sub (u, w, v);
		copy (res, u);

		return u;
	}
};

// Specialization of RandomDenseStream

class RandomDenseStreamGF2 : public VectorStream<BitVector>
{
    public:
	typedef BitVector Vector;

	RandomDenseStreamGF2 (const GF2 &F, uint32 seed, size_t n, size_t m = 0)
		: _MT (seed), _n (n), _m (m), _j (0)
	{}

	Vector &get (Vector &v) 
	{
		Vector::word_iterator i;

		if (_m > 0 && _j++ >= _m)
			return v;

		for (i = v.wordBegin (); i != v.wordEnd (); i++)
			*i = _MT.randomInt ();

		return v;
	}

	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

    private:
	MersenneTwister _MT;
	size_t          _n;
	size_t          _m;
	size_t          _j;
};

// Specialization of RandomSparseStream

template <class _Vector = Vector<GF2>::Sparse>
class RandomSparseStreamGF2 : public VectorStream<_Vector>
{
    public:
	typedef GF2 Field;
	typedef _Vector Vector;

	RandomSparseStreamGF2 (const GF2 &F, uint32 seed, double p, size_t n, size_t m = 0)
		: _MT (seed), _n (n), _m (m), _j (0)
	{ setP (p); }

	Vector &get (Vector &v);

	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

	void setP (double p)
	{
		linbox_check ((p >= 0.0) && (p <= 1.0)); 
		_p = p;
		_1_log_1mp   = 1 / log (1 - _p);
	}

    private:
	MersenneTwister _MT;
	size_t _n;
	double _p;
	double _1_log_1mp;
	size_t _m;
	size_t _j;
};

template <class _Vector>
_Vector &RandomSparseStreamGF2<_Vector>::get (_Vector &v)
{
	size_t i = (size_t) -1;
	double val;
	int skip;

	if (_m > 0 && _j++ >= _m)
		return v;

	v.clear ();

	while (1) {
		val = (double) _MT.randomDouble ();
		skip = (int) (ceil (log (val) * _1_log_1mp));

		if (skip <= 0)
			i++;
		else
			i += skip;

		if (i >= _n) break;

		v.push_back (i);
	}

	return v;
}

template <class Vector1, class Vector2>
inline bool &DotProductDomain<GF2>::dotSpecializedDD
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	linbox_check (v1.size () == v2.size ());

	uint32 t = 0;
	uint32 mask;
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();

	while (i != v1.wordEnd () - 1)
		t ^= *i++ & *j++;

	mask = (1 << (v1.size () & 31)) - 1;
	if (mask == 0) mask = 0xffffffff;

	t ^= *i & *j & mask;

	t ^= (t >> 16);
	t ^= (t >> 8);
	t ^= (t >> 4);
	t ^= (t >> 2);
	t ^= (t >> 1);

	return res = bool (t & 1);
}

template <class Vector1, class Vector2>
inline bool &DotProductDomain<GF2>::dotSpecializedDSP
	(bool          &res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	typename Vector2::const_iterator i;

	res = 0;

	for (i = v2.begin (); i != v2.end (); ++i)
		res ^= v1[*i];

	return res;
}

template <class Vector1, class Vector2>
inline BitVector::reference DotProductDomain<GF2>::dotSpecializedDD
	(BitVector::reference res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	linbox_check (v1.size () == v2.size ());

	uint32 t = 0;
	uint32 mask;
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();

	while (i != v1.wordEnd () - 1)
		t ^= *i++ & *j++;

	mask = (1 << (v1.size () & 31)) - 1;
	if (mask == 0) mask = 0xffffffff;

	t ^= *i & *j & mask;

	t ^= (t >> 16);
	t ^= (t >> 8);
	t ^= (t >> 4);
	t ^= (t >> 2);
	t ^= (t >> 1);

	return res = bool (t & 1);
}

template <class Vector1, class Vector2>
inline BitVector::reference DotProductDomain<GF2>::dotSpecializedDSP
	(BitVector::reference res,
	 const Vector1 &v1,
	 const Vector2 &v2) const
{
	typename Vector2::const_iterator i;

	res = 0;

	for (i = v2.begin (); i != v2.end (); ++i)
		res ^= v1[*i];

	return res;
}

template <class Vector, class Trait>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::DenseZeroOneVectorTag<Trait>) const
{
	typename Vector::const_iterator i;

	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i)
		os << *i << ' ';

	os << ']';

	return os;
}

template <class Vector, class Trait>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::SparseZeroOneVectorTag<Trait>) const
{
	typename Vector::const_iterator i;
	size_t idx = 0;

	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i) {
		while (++idx <= *i)
			os << 0 << ' ';

		os << 1 << ' ';
	}

	os << ']';

	return os;
}

template <class Vector, class Trait>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
						  VectorCategories::DenseZeroOneVectorTag<Trait>) const
{
	typename Vector::iterator i;
	char c;

	while (!isdigit (is >> c));

	is.unget (c);

	for (i = x.begin (); i != x.end (); ++i)
		is >> *i;

	return is;
}

template <class Vector, class Trait>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &os, const Vector &x,
						  VectorCategories::SparseZeroOneVectorTag<Trait>) const
{
	char c;
	size_t idx;

	while (!isdigit (is >> c));

	is.unget (c);
	x.clear ();

	while (1) {
		is >> c;

		if (!isdigit (c) && c != ' ') break;
		is.unget (c);
		is >> idx;
		x.push_back (idx);
	}

	return is;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag<Trait1>,
					     VectorCategories::SparseZeroOneVectorTag<Trait2>) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();
	size_t idx = 0;

	for (; j != v2.end (); ++j, ++i, ++idx) {
		while (idx < *j) {
			if (*i) return false;
			++idx;
			++i;
		}

		if (!*i) return false;
	}

	for (; i != v1.end (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector, class Trait>
bool VectorDomain<GF2>::isZeroSpecialized (const Vector &v,
					   VectorCategories::DenseZeroOneVectorTag<Trait>) const
{
	typename Vector::const_word_iterator i;

	for (i = v.wordBegin (); i != v.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::SparseZeroOneVectorTag<Trait1>,
					     VectorCategories::DenseZeroOneVectorTag<Trait2>) const
{
	typename Vector2::const_iterator i;
	size_t idx = 0;

	res.clear ();

	for (i = v.begin (); i != v.end (); ++i, ++idx)
		if (*i) res.push_back (idx);

	return res;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::DenseZeroOneVectorTag<Trait1>,
					     VectorCategories::SparseZeroOneVectorTag<Trait2>) const
{
	typename Vector2::const_iterator i;

	std::fill (res.wordBegin (), res.wordEnd (), 0);

	for (i = v.begin (); i != v.end (); ++i)
		res[*i] = 1;

	return res;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
bool &VectorDomain<GF2>::dotSpecialized (bool &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag<Trait1>,
					 VectorCategories::SparseZeroOneVectorTag<Trait2>) const
{
	typename Vector1::const_iterator i = v1.begin ();
	typename Vector2::const_iterator j = v2.begin ();
	res = false;

	while (i != v1.end () || j != v2.end ()) {
		while (i != v1.end () && (j == v2.end () || *i < *j)) { res = !res; ++i; }
		while (j != v2.end () && (i == v1.end () || *j < *i)) { res = !res; ++j; }
		if (i != v1.end () && j != v2.end () && *i == *j) { ++i; ++j; }
	}

	return res;
}

template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::DenseZeroOneVectorTag<Trait1>,
					    VectorCategories::DenseZeroOneVectorTag<Trait2>,
					    VectorCategories::DenseZeroOneVectorTag<Trait3>) const
{
	linbox_check (res.size () == y.size ());
	linbox_check (res.size () == x.size ());

	typename Vector1::word_iterator i = res.wordBegin ();
	typename Vector2::const_word_iterator j = y.wordBegin ();
	typename Vector3::const_word_iterator k = x.wordBegin ();

	for (; i != res.wordEnd (); ++i)
		*i = *j++ ^ *k++;

	return res;
}

template <class Vector1, class Trait1, class Vector2, class Trait2, class Vector3, class Trait3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::SparseZeroOneVectorTag<Trait1>,
					    VectorCategories::SparseZeroOneVectorTag<Trait2>,
					    VectorCategories::SparseZeroOneVectorTag<Trait3>) const
{
	typename Vector2::const_iterator i = y.begin ();
	typename Vector3::const_iterator j = x.begin ();

	res.clear ();

	while (i != y.end () || j != x.end ()) {
		while (i != y.end () && (j == x.end () || *i < *j)) { res.push_back (*i); ++i; }
		while (j != x.end () && (i == y.end () || *j < *i)) { res.push_back (*j); ++j; }
		if (i != y.end () && j != x.end () && *i == *j) { ++i; ++j; }
	}

	return res;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag<Trait1>,
					      VectorCategories::DenseZeroOneVectorTag<Trait2>) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::word_iterator i = y.wordBegin ();
	typename Vector2::const_word_iterator j = x.wordBegin ();

	for (; i != y.wordEnd (); ++i, ++j)
		*i ^= *j;

	return y;
}

template <class Vector1, class Trait1, class Vector2, class Trait2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag<Trait1>,
					      VectorCategories::SparseZeroOneVectorTag<Trait2>) const
{
	typename Vector2::const_iterator i;

	for (i = x.begin (); i != x.end (); ++i)
		y[*i] = !y[*i];

	return y;
}

// Specialization of MatrixVectorDomain for GF2
template <>
class MatrixVectorDomain<GF2>
{
    protected:
	MatrixVectorDomain (const GF2 &F) : _VD (F) {}

	template <class Vector1, class Matrix, class Vector2, class VectorTrait>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag<VectorTrait>) const;
	template <class Vector1, class Matrix, class Vector2, class VectorTrait>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::SparseZeroOneVectorTag<VectorTrait>) const;

	template <class Vector1, class VectorTrait1, class Matrix, class Vector2, class VectorTrait2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag<VectorTrait1>,
				    VectorCategories::DenseZeroOneVectorTag<VectorTrait2>) const;
	template <class Vector1, class VectorTrait1, class Matrix, class Vector2, class VectorTrait2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag<VectorTrait1>,
				    VectorCategories::SparseZeroOneVectorTag<VectorTrait2>) const;

	VectorDomain<GF2> _VD;
};

template <class Vector1, class Matrix, class Vector2, class VectorTrait>
Vector1 &MatrixVectorDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag<VectorTrait>) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector1::iterator j = w.begin ();

	for (; j != w.end (); ++j, ++i)
		_VD.dot (*j, v, *i);

	return w;
}

template <class Vector1, class Matrix, class Vector2, class VectorTrait>
Vector1 &MatrixVectorDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::SparseZeroOneVectorTag<VectorTrait>) const
{
	typename Matrix::ConstRowIterator i = A.rowBegin ();
	GF2::Element t;
	unsigned int idx = 0;

	w.clear ();

	for (; i != A.rowEnd (); ++i, ++idx) {
		_VD.dot (t, v, *i);

		if (t)
			w.push_back (t);
	}

	return w;
}

template <class Vector1, class VectorTrait1, class Matrix, class Vector2, class VectorTrait2>
Vector1 &MatrixVectorDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag<VectorTrait1>,
						     VectorCategories::DenseZeroOneVectorTag<VectorTrait2>) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstColIterator i = A.colBegin ();
	typename Vector2::const_iterator j = v.begin ();

	_VD.subin (w, w);

	for (; j != v.end (); ++j, ++i)
		_VD.axpyin (w, *j, *i);

	return w;
}

template <class Vector1, class VectorTrait1, class Matrix, class Vector2, class VectorTrait2>
Vector1 &MatrixVectorDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag<VectorTrait1>,
						     VectorCategories::SparseZeroOneVectorTag<VectorTrait2>) const
{
	linbox_check (A.rowdim () == w.size ());

	typename Vector2::const_iterator j = v.begin ();

	_VD.subin (w, w);

	for (; j != v.end (); ++j) {
		typename Matrix::ConstColIterator i = A.colBegin () + *j;
		_VD.axpyin (w, true, *i);
	}

	return w;
}

// Specialization of diagonal for GF2
template <>
class Diagonal<GF2, Vector<GF2>::Dense, VectorTraits<Vector<GF2>::Dense>::VectorCategory>
	: public BlackboxArchetype<Vector<GF2>::Dense>
{
    public:

	typedef GF2                       Field;
	typedef Vector<GF2>::Dense        Vector;
	typedef BlackboxArchetype<Vector> Blackbox;
	typedef bool                      Element;

	Diagonal (const Field &F, const BitVector &y)
		: _v (y) 
	{}

	Blackbox *clone() const
		{ return new Diagonal (*this); }

	Vector& apply (Vector& y, const Vector& x) const
	{
		linbox_check (y.size () == x.size ());
		linbox_check (y.size () == _v.size ());

		BitVector::word_iterator i = y.wordBegin ();
		BitVector::const_word_iterator j1 = x.wordBegin (), j2 = _v.wordBegin ();

		for (; i != y.wordEnd (); ++i, ++j1, ++j2)
			*i = *j1 & *j2;

		return y;
	}

	Vector& applyTranspose (Vector& y, const Vector& x) const { return apply (y, x); }
	size_t rowdim () const { return _v.size (); } 
	size_t coldim () const { return _v.size (); } 

    private:

	// Bit vector of elements
	BitVector _v;
    
}; // template <Field, Vector> class Diagonal<DenseVectorTag>

} // namespace LinBox

#endif // __FIELD_GF2_INL
