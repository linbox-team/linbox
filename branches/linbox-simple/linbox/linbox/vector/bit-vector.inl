/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/vector/bit-vector.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BIT_VECTOR_INL
#define __BIT_VECTOR_INL

#include <stdexcept>

#include "linbox/vector/vector-traits.h"
#include "linbox/vector/bit-vector.h"

namespace std 
{
	template <>
	struct iterator_traits<LinBox::BitVector::iterator>
	{
		typedef random_access_iterator_tag iterator_category;
		typedef LinBox::BitVector::reference reference;
		typedef bool *pointer;
		typedef bool value_type;
		typedef long difference_type;
	};

	template <>
	struct iterator_traits<LinBox::BitVector::const_iterator>
	{
		typedef random_access_iterator_tag iterator_category;
		typedef LinBox::BitVector::const_reference reference;
		typedef const bool *pointer;
		typedef bool value_type;
		typedef long difference_type;
	};
}

namespace LinBox
{

class BitVector::reference
{
    public:

	reference (std::vector<unsigned long>::iterator word, uint8 position)
		: _word (word), _pos (position) {}

	~reference () {}

	reference &operator = (reference &a) 
		{ return *this = (bool) a; }

	reference &operator = (bool v) 
		{ *_word = v ? (*_word | (1 << _pos)) : (*_word & ~(1 << _pos)); return *this; }

	reference &operator &= (reference &a) 
		{ *_word &= ~(1 << _pos) | (a.get_bit () << (_pos - a._pos)); return *this; }

	reference &operator &= (bool v) 
		{ *_word &= ~(1 << _pos) | (v & (1 << _pos)); return *this; }

	reference &operator |= (reference &a) 
		{ *_word |= a.get_bit () << (_pos - a._pos); return *this; }

	reference &operator |= (bool v) 
		{ *_word |= v & (1 << _pos); return *this; }

	reference &operator ^= (reference &a) 
		{ *_word ^= a.get_bit () << (_pos - a._pos); return *this; }

	reference &operator ^= (bool v) 
		{ *_word ^= v & (1 << _pos); return *this; }

	operator bool (void) const
		{ return (*_word >> _pos) & 1; return *this; }

    private:
	friend class iterator;
	friend class const_iterator;

	unsigned long neg_mask_word (void) { return *_word & ~(1 << _pos); }
	unsigned long get_bit ()           { return *_word & (1 << _pos); }

	std::vector<unsigned long>::iterator _word;
	uint8                         _pos;
};

std::istream &operator >> (std::istream &is, BitVector::reference &a) 
	{ bool v; is >> v; a = v; return is; }

std::ostream &operator << (std::ostream &os, BitVector::reference &a) 
	{ os << bool (a); return os; }

class BitVector::const_reference
{
    public:

	const_reference (std::vector<unsigned long>::const_iterator word, uint8 position)
		: _word (word), _pos (position) {}

	~const_reference () {}

	operator bool (void) const
		{ return (*_word >> _pos) & 1; }

    private:
	friend class const_iterator;

	std::vector<unsigned long>::const_iterator _word;
	uint8                               _pos;
};

std::ostream &operator << (std::ostream &os, BitVector::const_reference &a) 
	{ os << bool (a); return os; }

class BitVector::iterator : public std::iterator <std::random_access_iterator_tag, bool>
{
    public:

	typedef std::iterator_traits<iterator>::iterator_category iterator_category;
	typedef std::iterator_traits<iterator>::reference reference;
	typedef std::iterator_traits<iterator>::pointer pointer;
	typedef std::iterator_traits<iterator>::value_type value_type;
	typedef std::iterator_traits<iterator>::difference_type difference_type;

	iterator () : _ref (std::vector<unsigned long>::iterator (), 0) {}
	iterator (std::vector<unsigned long>::iterator word, uint8 position) : _ref (word, position) {}
	iterator (const iterator &i) : _ref (i._ref._word, i._ref._pos) {}

	iterator &operator = (const iterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	iterator &operator ++ () 
	{
		if (++_ref._pos > 31) {
			++_ref._word;
			_ref._pos = 0;
		}

		return *this;
	}

	iterator operator ++ (int) 
	{
		iterator tmp (*this);
		++*this;
		return tmp;
	}

	iterator operator + (difference_type i) const
	{
		std::vector<unsigned long>::iterator new_word = _ref._word + (i >> 5);
		uint8 new_pos = _ref._pos + (i & 0x1f);

		new_word += new_pos >> 5;
		new_pos &= 0x1f;

		return iterator (new_word, new_pos);
	}

	iterator &operator += (difference_type i) 
	{
		_ref._word += i >> 5;
		_ref._pos  += i & 0x1f;
		_ref._word += _ref._pos >> 5;
		_ref._pos  &= 0x1f;
		return *this;
	}

	iterator &operator -- () 
	{
		if (--_ref._pos > 31) {
			--_ref._word;
			_ref._pos = 31;
		}

		return *this;
	}

	iterator operator -- (int) 
	{
		iterator tmp (*this);
		--*this;
		return tmp;
	}

	iterator operator - (difference_type i) const
		{ return *this + -i; }

	iterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (iterator &i) const 
		{ return (_ref._word - i._ref._word) * 32 + (_ref._pos - i._ref._pos); }

	reference operator [] (int i) 
		{ return *(*this + i); }

	reference operator * () 
		{ return _ref; }

	bool operator == (const iterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator != (const iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

    private:
	friend class const_iterator;

	reference _ref;
};
 
class BitVector::const_iterator : public std::iterator <std::random_access_iterator_tag, bool>
{
    public:

	typedef std::iterator_traits<const_iterator>::iterator_category iterator_category;
	typedef std::iterator_traits<const_iterator>::reference reference;
	typedef std::iterator_traits<const_iterator>::pointer pointer;
	typedef std::iterator_traits<const_iterator>::value_type value_type;
	typedef std::iterator_traits<const_iterator>::difference_type difference_type;

	const_iterator () : _ref (std::vector<unsigned long>::const_iterator (), 0) {}
	const_iterator (std::vector<unsigned long>::const_iterator word, uint8 position) : _ref (word, position) {}
	const_iterator (const const_iterator &i) : _ref (i._ref._word, i._ref._pos) {}

	const_iterator &operator = (const const_iterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	const_iterator &operator = (const iterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	const_iterator &operator ++ () 
	{
		if (++_ref._pos > 31) {
			++_ref._word;
			_ref._pos = 0;
		}

		return *this;
	}

	const_iterator operator ++ (int) 
	{
		const_iterator tmp (*this);
		++*this;
		return tmp;
	}

	const_iterator operator + (int i) const
	{
		std::vector<unsigned long>::const_iterator new_word = _ref._word + (i >> 5);
		uint8 new_pos = _ref._pos + (i & 0x1f);

		new_word += new_pos >> 5;
		new_pos &= 0x1f;

		return const_iterator (new_word, new_pos);
	}

	const_iterator &operator += (int i) 
	{
		_ref._word += i >> 5;
		_ref._pos  += i & 0x1f;
		_ref._word += _ref._pos >> 5;
		_ref._pos  &= 0x1f;
		return *this;
	}

	const_iterator &operator -- () 
	{
		if (--_ref._pos > 31) {
			--_ref._word;
			_ref._pos = 31;
		}

		return *this;
	}

	const_iterator operator -- (int) 
	{
		const_iterator tmp (*this);
		--*this;
		return tmp;
	}

	const_iterator operator - (difference_type i) const 
		{ return *this + -i; }

	const_iterator &operator -= (difference_type i) 
		{ return *this += -i; }

	difference_type operator - (const_iterator &i) const 
		{ return (_ref._word - i._ref._word) * 32 + (_ref._pos - i._ref._pos); }

	reference operator [] (difference_type i) const
		{ return *(*this + i); }

	reference operator * () const
		{ return _ref; }

	bool operator == (const const_iterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator == (const iterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator != (const const_iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

	bool operator != (const iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

    private:

	const_reference _ref;
};

inline BitVector::iterator BitVector::begin (void)
	{ return iterator (_v.begin (), 0); }

inline BitVector::const_iterator BitVector::begin (void) const
	{ return const_iterator (_v.begin (), 0); }

inline BitVector::iterator BitVector::end (void)
{
	if (_size & 0x1f == 0)
		return iterator (_v.end (), 0);
	else
		return iterator (_v.end () - 1, _size & 0x1f);
}

inline BitVector::const_iterator BitVector::end (void) const
{
	if (_size & 0x1f == 0)
		return const_iterator (_v.end (), 0);
	else
		return const_iterator (_v.end () - 1, _size & 0x1f);
}

inline BitVector::reverse_iterator BitVector::rbegin (void)
	{ return reverse_iterator (end () - 1); }

inline BitVector::const_reverse_iterator BitVector::rbegin (void) const
	{ return const_reverse_iterator (end () - 1); }

inline BitVector::reverse_iterator BitVector::rend (void)
	{ return reverse_iterator (begin () - 1); }

inline BitVector::const_reverse_iterator BitVector::rend (void) const
	{ return const_reverse_iterator (begin () - 1); }

inline BitVector::reference BitVector::operator[] (BitVector::size_type n)
	{ return *(begin () + n); }

inline BitVector::const_reference BitVector::operator[] (BitVector::size_type n) const
	{ return *(begin () + n); }

BitVector::reference BitVector::at (BitVector::size_type n)
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

BitVector::const_reference BitVector::at (BitVector::size_type n) const
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

inline BitVector::reference BitVector::front (void)
	{ return reference (_v.begin (), 0); }

inline BitVector::const_reference BitVector::front (void) const
	{ return const_reference (_v.begin (), 0); }

inline BitVector::reference BitVector::back (void)
{
	if (_size & 0x1f == 0)
		return reference (_v.end (), 0);
	else
		return reference (_v.end () - 1, _size & 0x1f);
}

inline BitVector::const_reference BitVector::back (void) const
{
	if (_size & 0x1f == 0)
		return const_reference (_v.end (), 0);
	else
		return const_reference (_v.end () - 1, _size & 0x1f);
}

template<class Container>
BitVector &BitVector::operator = (const Container &v)
{
	typename Container::const_iterator i;
	typename Container::const_iterator i_end = v.begin () + (v.size () >> 5);
	std::vector<unsigned long>::iterator j;
	unsigned int idx;

	_v.resize ((v.size () >> 5) + ((v.size () & 0x1F) ? 1 : 0));

	for (j = _v.begin (); i != i_end; ++j) {
		*j = 0;
		for (idx = 0; idx < 32; ++idx, ++i) {
			*j <<= 1;
			*j |= *i & 1;
		}
	}

	if (v.size () & 0x1f) {
		*j = 0;

		for (idx = 0; idx < v.size () & 0x1f; ++idx) {
			*j <<= 1;
			*j |= *i & 1;
		}
	}

	_size = v.size ();

	return *this;
}

void BitVector::resize (BitVector::size_type new_size, bool val)
	{ _v.resize ((new_size >> 5) + ((new_size & 0x1F) ? 1 : 0), val ? 0xffffffff : 0); _size = new_size; }

bool BitVector::operator == (const BitVector &v) const
{
	const_word_iterator i, j;
	unsigned long mask;

	if (_size != v._size) return false;

	for (i = wordBegin (), j = v.wordBegin (); i != wordEnd () - 1; ++i, ++j)
		if (*i != *j) return false;

	mask = (1 << (_size & (8 * sizeof (unsigned long) - 1))) - 1;
	if (mask == 0) mask = (unsigned long) -1;

	if ((*i & mask) == (*j & mask))
		return true;
	else
		return false;
}

/* 
namespace VectorWrapper 
{
	template <class Field, class Vector, class Trait>
	inline BitVector::reference refSpecialized
		(Vector &v, size_t i, VectorCategories::DenseZeroOneVectorTag<Trait>)
	{ return v[i]; }

	template <class Field, class Vector, class Trait>
	inline BitVector::const_reference constRefSpecialized
		(const Vector &v, size_t i, VectorCategories::DenseZeroOneVectorTag<Trait>)
	{ return v[i]; }

}
*/

} // namespace LinBox

#endif // __BIT_VECTOR_INL
