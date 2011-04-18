/* linbox/field/gf2.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen, Dumas, bds
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_field_gf2_INL
#define __LINBOX_field_gf2_INL

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

#include <cctype> //isdigit

template<typename Vector>
std::ostream& afficheVector (std::ostream& o, const Vector& C) {
          for(typename Vector::const_iterator refs =  C.begin();
                                refs != C.end() ;
                                      ++refs )
                          o << (*refs) << " " ;
            return o;
}


namespace LinBox 
{ 

// Specialization of canonical vector types

template <>
class RawVector<bool>
{
    public:
    typedef BitVector Dense;
    typedef std::vector<size_t> Sparse;
    typedef std::vector<size_t> SparseSeq;
    typedef std::vector<size_t> SparseMap;
    typedef std::vector<size_t> SparsePar;
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
		: VectorDomainBase<GF2> (VD._F), DotProductDomain<GF2> (VD._F)
	{}

	VectorDomain &operator = (const VectorDomain &) { return *this; }

	const GF2 &field () const { return _F; }
    
	template <class Vector>
	inline std::ostream &write (std::ostream &os, const Vector &x) const
		{ return writeSpecialized (os, x, typename VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector>
	inline std::istream &read (std::istream &is, Vector &x) const
		{ return readSpecialized (is, x, typename VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v) const
		{ return copySpecialized (res, v,
					  typename VectorTraits<Vector1>::VectorCategory (),
					  typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &copy (Vector1 &res, const Vector2 &v, size_t i, size_t len = 0) const
		{ return copySpecialized (res, v, i, len,
					  typename VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline bool areEqual (const Vector1 &v1, const Vector2 &v2) const
		{ return areEqualSpecialized (v1, v2,
					      typename VectorTraits<Vector1>::VectorCategory (),
					      typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector>
	inline bool isZero (const Vector &v) const
		{ return isZeroSpecialized (v, typename VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dot (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dotSpecialized (res, v1, v2,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline BitVector::reference dot (BitVector::reference res, const Vector1 &v1, const Vector2 &v2) const
		{ return dotSpecialized (res, v1, v2,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Element &dotprod (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{ return dot (res, v1, v2); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &add (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory (),
					 typename VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &addin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   typename VectorTraits<Vector1>::VectorCategory (),
					   typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &sub (Vector1 &res, const Vector2 &y, const Vector3 &x) const
		{ return addSpecialized (res, y, x,
					 typename VectorTraits<Vector1>::VectorCategory (),
					 typename VectorTraits<Vector2>::VectorCategory (),
					 typename VectorTraits<Vector3>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &subin (Vector1 &y, const Vector2 &x) const
		{ return addinSpecialized (y, x,
					   typename VectorTraits<Vector1>::VectorCategory (),
					   typename VectorTraits<Vector2>::VectorCategory ()); }

	template <class Vector1, class Vector2>
	inline Vector1 &neg (Vector1 &res, const Vector2 &x) const
		{ copy (res, x); return res; }

	template <class Vector>
	inline Vector &negin (Vector &y) const
		{ return y; }

	template <class Vector1, class Vector2>
	inline Vector1 &mul (Vector1 &res, const Vector2 &x, const Element a) const
		{ return mulSpecialized (res, x, a, typename VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector>
	inline Vector &mulin (Vector &x, const Element a) const
		{ return mulinSpecialized (x, a, typename VectorTraits<Vector>::VectorCategory ()); }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &axpy (Vector1 &res, const Element a, const Vector2 &x, const Vector3 &y) const
		{ if (a) add (res, x, y); else this->copy (res, y); return res; }

	template <class Vector1, class Vector2>
	inline Vector1 &axpyin (Vector1 &y, const Element a, const Vector2 &x) const
		{ if (a) addin (y, x); return y; }

	VectorDomain (const GF2 &F)
		: VectorDomainBase<GF2> (F), DotProductDomain<GF2> (F)
	{}


	// Specialized function implementations
	template <class Vector> 
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	std::ostream &writeSpecialized (std::ostream &os, const Vector &x,
					VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	std::istream &readSpecialized (std::istream &is, const Vector &x,
				       VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const;
	
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	inline bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ return areEqual (v2, v1); }
	template <class Vector1, class Vector2>
	bool areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
				  VectorCategories::SparseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
    

	template <class Vector>
	bool isZeroSpecialized (const Vector &v, VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector>
	inline bool isZeroSpecialized (const Vector &v,
				       VectorCategories::SparseZeroOneVectorTag) const
		{ return v.empty (); }

	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::DenseZeroOneVectorTag,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ std::copy (v.wordBegin (), v.wordEnd (), res.wordBegin ()); return res; }

	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len, VectorCategories::DenseZeroOneVectorTag) const
	{
		std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
		return res;
	}
	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v, size_t i, size_t len, VectorCategories::DenseVectorTag) const
	{
		std::copy (v.begin (), (len == 0) ? v.end () : v.begin () + len, res.begin () + i);
		return res;
	}



	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::SparseZeroOneVectorTag,
				  VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
				  VectorCategories::DenseZeroOneVectorTag,
				  VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	inline Vector1 &copySpecialized (Vector1 &res, const Vector2 &v,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const
		{ res = v; return res; }

	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::SparseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2); }
	template <class Vector1, class Vector2>
	inline Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::SparseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1); }
	template <class Vector1, class Vector2>
	Element &dotSpecialized (Element &res, const Vector1 &v1, const Vector2 &v2,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDD (res, v1, v2); }
	template <class Vector1, class Vector2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::DenseZeroOneVectorTag,
					VectorCategories::SparseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v1, v2); }
	template <class Vector1, class Vector2>
	inline BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
					VectorCategories::SparseZeroOneVectorTag,
					VectorCategories::DenseZeroOneVectorTag) const
		{ return DotProductDomain<GF2>::dotSpecializedDSP (res, v2, v1); }
	template <class Vector1, class Vector2>
	BitVector::reference dotSpecialized (BitVector::reference res, const Vector1 &v1, const Vector2 &v2,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::DenseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const
		{ copy (res, y); addin (res, x); }
	template <class Vector1, class Vector2, class Vector3>
	Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag,
				 VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag,
				   VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::DenseZeroOneVectorTag,
				   VectorCategories::SparseZeroOneVectorTag) const;
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag,
				   VectorCategories::DenseZeroOneVectorTag) const
		{ Vector1 xp, res; copy (xp, x); add (res, y, xp); copy (y, res); return y; }
	template <class Vector1, class Vector2>
	Vector1 &addinSpecialized (Vector1 &y, const Vector2 &x,
				   VectorCategories::SparseZeroOneVectorTag,
				   VectorCategories::SparseZeroOneVectorTag) const
		{ Vector1 res; add (res, y, x); this->copy (y, res); return y; }

	template <class Vector1, class Vector2>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::DenseZeroOneVectorTag ) const
		{ if (a) this->copy (res, x); else std::fill (res.wordBegin (), res.wordEnd (), 0); return res; }
	template <class Vector1, class Vector2>
	Vector1 &mulSpecialized (Vector1 &res, const Vector2 &x, const Element a,
				 VectorCategories::SparseZeroOneVectorTag ) const
		{ if (a) this->copy (res, x); else res.clear (); return res; }

	template <class Vector>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::DenseZeroOneVectorTag) const
		{ if (!a) std::fill (x.wordBegin (), x.wordEnd (), 0); return x; }

	template <class Vector>
	inline Vector &mulinSpecialized (Vector &x, const Element a,
					 VectorCategories::SparseZeroOneVectorTag ) const
		{ if (!a) x.clear (); return x; }

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag) const
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

	template <class Vector1, class Vector2, class Vector3>
	inline Vector1 &subSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag,
					VectorCategories::GenericVectorTag) const
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
template<size_t bitsize> struct MTrandomInt {
    template<typename M32Twister>
    unsigned __LINBOX_INT32 operator() (M32Twister& MT) const {
        return MT.randomInt();
    }
};    

template<> struct MTrandomInt<64> {
    template<typename M32Twister>
    unsigned __LINBOX_INT64 operator() (M32Twister& MT) const {
        unsigned __LINBOX_INT64 tmp = MT.randomInt();
        tmp <<=32;
        return tmp += MT.randomInt();
    }
};

class RandomDenseStreamGF2 : public VectorStream<BitVector>
{
    public:
	typedef BitVector Vector;

	RandomDenseStreamGF2 (const GF2 &, uint32 seed, size_t n, size_t m = 0)
		: MT (seed), _n (n), _m (m), _j (0)
	{}

	Vector &get (Vector &v) 
	{
		Vector::word_iterator i;

		if (_m > 0 && _j++ >= _m)
			return v;

		for (i = v.wordBegin (); i != v.wordEnd (); i++)
			*i = MTrandomInt<__LINBOX_BITSOF_LONG>()(MT);
                
                const size_t zeroing = __LINBOX_BITSOF_LONG - (v.size() % __LINBOX_BITSOF_LONG);
                *(v.wordRbegin()) <<= zeroing;
                *(v.wordRbegin()) >>= zeroing;
		return v;
	}

	size_t size () const { return _m; }
	size_t pos () const { return _j; }
	size_t dim () const { return _n; }
	operator bool () const { return _m == 0 || _j < _m; }
	void reset () { _j = 0; }

    private:
	MersenneTwister MT;
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

	RandomSparseStreamGF2 (const GF2 &, uint32 seed, double p, size_t n, size_t m = 0)
		: MT (seed), _n (n), _m (m), _j (0)
	{ setP (p); }

    	RandomSparseStreamGF2 (const GF2 &F, const GF2RandIter& r, double p, size_t n, size_t m = 0)
		: MT (r.getMT()), _n (n), _m (m), _j (0)
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
	MersenneTwister MT;
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
		val = (double) MT.randomDouble ();
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

	unsigned long t = 0;
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();

	while (i != v1.wordEnd () - 1)
		t ^= *i++ & *j++;
        
        const size_t zeroing = __LINBOX_BITSOF_LONG - (v1.size() % __LINBOX_BITSOF_LONG);
        unsigned long lastdot = *i & *j;
        lastdot <<= zeroing;
        lastdot >>= zeroing;
        
        t ^= lastdot;
        return res = __LINBOX_PARITY(t);
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
    bool tmp;
    return res = dotSpecializedDD(tmp, v1, v2);
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

template <class Vector>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::DenseZeroOneVectorTag) const
{
	

// TO BE REMOVED
	os << "writeSpec DenseZO, of size " << x.size() << ' ';

	os << "[ ";

	for (typename Vector::const_iterator i = x.begin (); i != x.end (); ++i)
		os << *i << ' ';

	os << ']';

	os << "( ";

	for (typename Vector::const_word_iterator i = x.wordBegin (); i != x.wordEnd (); ++i)
		os << *i << ' ';

	os << ')';

	return os;
}

template <class Vector>
std::ostream &VectorDomain<GF2>::writeSpecialized (std::ostream &os, const Vector &x,
						   VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector::const_iterator i;
	size_t idx = 0;

// TO BE REMOVED
	os << "writeSpec SparseZO, of size " << x.size() << ' ';
	os << "[ ";

	for (i = x.begin (); i != x.end (); ++i) {
		while (++idx <= *i)
			os << 0 << ' ';

		os << 1 << ' ';
	}
	os << ']';

	return os;
}

template <class Vector>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
						  VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector::iterator i;
	char c;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();

	for (i = x.begin (); i != x.end (); ++i)
		is >> *i;

	return is;
}

template <class Vector>
std::istream &VectorDomain<GF2>::readSpecialized (std::istream &is, const Vector &x,
						  VectorCategories::SparseZeroOneVectorTag) const
{
	char c;
	size_t idx;

	do { is >> c ; } while (!std::isdigit (c));

	is.unget ();
	x.clear ();

	while (1) {
		is >> c;

		if (!std::isdigit (c) && c != ' ') break;
		is.unget ();
		is >> idx;
		x.push_back (idx);
	}

	return is;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
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

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector1::const_word_iterator i = v1.wordBegin ();
	typename Vector2::const_word_iterator j = v2.wordBegin ();
	for (; j != v2.wordEnd (); ++j, ++i)
		if (*i != *j) return false;
	return true;
}

template <class Vector1, class Vector2>
bool VectorDomain<GF2>::areEqualSpecialized (const Vector1 &v1, const Vector2 &v2,
					     VectorCategories::SparseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{ return v1 == v2;}

template <class Vector>
bool VectorDomain<GF2>::isZeroSpecialized (const Vector &v,
					   VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector::const_word_iterator i;

	for (i = v.wordBegin (); i != v.wordEnd (); ++i)
		if (*i) return false;

	return true;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::SparseZeroOneVectorTag,
					     VectorCategories::DenseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i;
	size_t idx = 0;

	res.clear ();

	for (i = v.begin (); i != v.end (); ++i, ++idx)
		if (*i) res.push_back (idx);

	return res;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::copySpecialized (Vector1 &res, const Vector2 &v,
					     VectorCategories::DenseZeroOneVectorTag,
					     VectorCategories::SparseZeroOneVectorTag) const
{
    	size_t sparsesize = *(v.rbegin());
    	if (sparsesize > res.size()) res.resize( *(v.rbegin()) );
	std::fill (res.wordBegin (), res.wordEnd (), 0);

	for (typename Vector2::const_iterator i = v.begin (); 
             i != v.end (); 
             ++i)
        	res[*i] = true;
	return res;
}

template <class Vector1, class Vector2>
bool &VectorDomain<GF2>::dotSpecialized (bool &res, const Vector1 &v1, const Vector2 &v2,
					 VectorCategories::SparseZeroOneVectorTag,
					 VectorCategories::SparseZeroOneVectorTag) const
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

template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::DenseZeroOneVectorTag,
					    VectorCategories::DenseZeroOneVectorTag) const
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

template <class Vector1, class Vector2, class Vector3>
Vector1 &VectorDomain<GF2>::addSpecialized (Vector1 &res, const Vector2 &y, const Vector3 &x,
					    VectorCategories::SparseZeroOneVectorTag,
					    VectorCategories::SparseZeroOneVectorTag,
					    VectorCategories::SparseZeroOneVectorTag) const
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

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (y.size () == x.size ());

	typename Vector1::word_iterator i = y.wordBegin ();
	typename Vector2::const_word_iterator j = x.wordBegin ();

	for (; i != y.wordEnd (); ++i, ++j)
		*i ^= *j;

	return y;
}

template <class Vector1, class Vector2>
Vector1 &VectorDomain<GF2>::addinSpecialized (Vector1 &y, const Vector2 &x,
					      VectorCategories::DenseZeroOneVectorTag,
					      VectorCategories::SparseZeroOneVectorTag) const
{
	typename Vector2::const_iterator i;

	for (i = x.begin (); i != x.end (); ++i)
		y[*i] = !y[*i];

	return y;
}

// Specialization of MatrixDomain for GF2
template <>
class MatrixDomain<GF2>
{
    public:
	MatrixDomain (const GF2 &F) : _VD (F) {}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &vectorMul (Vector1 &w, const Matrix &A, const Vector2 &v) const
        { return mulSpecialized (w, A, v, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::RowMatrixTag) const
		{ return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::ColMatrixTag) const
		{ return mulColSpecialized (w, A, v,
					    typename VectorTraits<Vector1>::VectorCategory (),
					    typename VectorTraits<Vector2>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::RowColMatrixTag) const
		{ return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::SparseZeroOneVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag,
				    VectorCategories::DenseZeroOneVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseZeroOneVectorTag,
				    VectorCategories::SparseZeroOneVectorTag) const;

	VectorDomain<GF2> _VD;
};

template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag) const
{
	linbox_check (A.coldim () == v.size ());
	linbox_check (A.rowdim () == w.size ());

	typename Matrix::ConstRowIterator i = A.rowBegin ();
	typename Vector1::iterator j = w.begin ();

	for (; j != w.end (); ++j, ++i)
		_VD.dot (*j, v, *i);

	return w;
}

template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<GF2>::mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::SparseZeroOneVectorTag) const
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

template <class Vector1, class Matrix, class Vector2 >
Vector1 &MatrixDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag,
						     VectorCategories::DenseZeroOneVectorTag) const
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

template <class Vector1, class Matrix, class Vector2>
Vector1 &MatrixDomain<GF2>::mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
						     VectorCategories::DenseZeroOneVectorTag,
						     VectorCategories::SparseZeroOneVectorTag) const
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
class Diagonal<GF2, VectorTraits<Vector<GF2>::Dense>::VectorCategory>
	: public BlackboxArchetype
{
    public:

	typedef GF2                       Field;
	typedef Vector<GF2>::Dense        Vector;
	typedef BlackboxArchetype         Blackbox;
	typedef bool                      Element;

	Diagonal (const Field &, const BitVector &y)
		: _v (y) 
	{}

        /// The field.	
        const Field& field() const {return *(new GF2());}

	Blackbox *clone() const
		{ return new Diagonal (*this); }


	template <class OutVector, class InVector>
        OutVector& apply (OutVector& y, const InVector& x) const
        {
            linbox_check (y.size () == x.size ());
            linbox_check (y.size () == _v.size ());
            typename InVector::const_iterator j1 = x.begin();
            typename OutVector::iterator i = y.begin();
            BitVector::const_iterator j2 = _v.begin();
            for (; i != y.end (); ++i, ++j1, ++j2)
                *i = *j1 & *j2;
            return y;
        }       

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

	template <class OutVector, class InVector>
        OutVector& applyTranspose (OutVector& y, const InVector& x) const 
        { return apply (y, x); }
    
	size_t rowdim () const { return _v.size (); } 
	size_t coldim () const { return _v.size (); } 

                /** Get an entry and store it in the given value
                 * @param x Element in which to store result
                 * @param i Row index
                 * @param j Column index
                 * @return Reference to x
                 */
        Element &getEntry (Element &x, size_t i, size_t j) const {
                return (i==j?x=this->_v[i]:x=false);
        }

    private:

	// Bit vector of elements
	BitVector _v;
    
}; // template <Field, Vector> class Diagonal<DenseVectorTag>

} // namespace LinBox

#include "linbox/switch/cekstv.h"
namespace LinBox 
{ 
// Specialization of Butterfly switch object 
template <>
class CekstvSwitch<GF2>
{
    public:
    	typedef GF2 Field;
	/// Typedef
	typedef Field::Element Element;
	typedef CekstvSwitch<Field> Self_t;
	typedef CekstvSwitchFactory<Field> Factory;

	/** Constructor from a field and a field element.
	 * @param F field in which arithmetic is done
	 * @param switches vector of switches
	 */
	CekstvSwitch (const Field::Element &a)
		: _a (a) 
	{}

	~CekstvSwitch () {}

	bool apply (const Field &F, Element &x, Element &y) const {
            F.axpyin (x, _a, y);
            F.addin (y, x);
            return true;
        }
    
	bool applyTranspose (const Field &F, Element &x, Element &y) const {
            F.addin (x, y);
            F.axpyin (y, _a, x);
            return true;
        }

	bool apply (const Field &F, std::_Bit_reference x, std::_Bit_reference y) const {
            F.axpyin (x, _a, y);
            F.addin (y, x);
            return true;
        }
    
	bool applyTranspose (const Field &F, std::_Bit_reference x, std::_Bit_reference y) const {
            F.addin (x, y);
            F.axpyin (y, _a, x);
            return true;
        }

        template<typename _Tp1>
        struct rebind
        { 
            typedef CekstvSwitch<_Tp1> other;

                // special rebind operator() with two fields, 
                // indeed local field is not stored in the switch
            void operator() (other *& Ap, const Self_t& A, const _Tp1& T, const Field& F) {
                typename _Tp1::Element u;
                Hom<Field, _Tp1>(F,T).image(u, A._a);
                Ap = new other(u);
            }
        };
    

   private:

	// Parameter of this 2x2 block
	Field::Element _a;
};


} // namespace LinBox

#endif // __LINBOX_field_gf2_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
