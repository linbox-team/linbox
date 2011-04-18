/* linbox/vector/bit-vector.inl
 * Copyright (C) 2003 Bradford Hovinen
 *
 * -------------------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_bit_vector_INL
#define __LINBOX_bit_vector_INL

#include <stdexcept>
#include <vector>


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
		{ 
                    *_word = v ? (*_word | (1UL << _pos)) : (*_word & ~(1UL << _pos));
                    return *this;
}

	reference &operator &= (reference &a) 
		{ *_word &= ~(1UL << _pos) | (a.get_bit () << (_pos - a._pos)); return *this; }

	reference &operator &= (bool v) 
		{ *_word &= ~(1UL << _pos) | (v & (1UL << _pos)); return *this; }

	reference &operator |= (reference &a) 
		{ *_word |= a.get_bit () << (_pos - a._pos); return *this; }

	reference &operator |= (bool v) 
		{ *_word |= v & (1UL << _pos); return *this; }

	reference &operator ^= (reference &a) 
		{ *_word ^= a.get_bit () << (_pos - a._pos); return *this; }

	reference &operator ^= (bool v) 
		{ *_word ^= v & (1UL << _pos); return *this; }

	operator bool (void) const
		{ return (*_word >> _pos) & 1UL; }

    private:
	friend class iterator;
	friend class const_iterator;
	friend class const_reference;

	unsigned long neg_mask_word (void) { return *_word & ~(1UL << _pos); }
	unsigned long get_bit ()           { return *_word & (1UL << _pos); }

	std::vector<unsigned long>::iterator _word;
	uint8                         _pos;
};

inline std::istream &operator >> (std::istream &is, BitVector::reference &a) 
	{ bool v; is >> v; a = v; return is; }

inline std::ostream &operator << (std::ostream &os, BitVector::reference &a) 
	{ os << bool (a); return os; }

class BitVector::const_reference
{
    public:

	const_reference (BitVector::reference r)
		: _word (r._word), _pos (r._pos) {}

	const_reference (std::vector<unsigned long>::iterator word, uint8 position)
		: _word (word), _pos (position) {}

	const_reference (std::vector<unsigned long>::const_iterator word, uint8 position)
		: _word (word), _pos (position) {}

	~const_reference () {}

	operator bool (void) const
		{ return (*_word >> _pos) & 1UL; }

    private:
	friend class const_iterator;

	std::vector<unsigned long>::const_iterator _word;
	uint8                               _pos;
};

inline std::ostream &operator << (std::ostream &os, BitVector::const_reference &a) 
	{ os << bool (a); return os; }

// class BitVector::iterator : public std::iterator <std::random_access_iterator_tag, bool>
class BitVector::iterator : public std::_Bit_iterator
{
    public:

	typedef std::iterator_traits<iterator>::iterator_category iterator_category;
	typedef std::iterator_traits<iterator>::reference reference;
	typedef std::iterator_traits<iterator>::pointer pointer;
	typedef std::iterator_traits<iterator>::value_type value_type;
	typedef std::iterator_traits<iterator>::difference_type difference_type;

	iterator () : _ref (std::vector<unsigned long>::iterator (), 0UL) {}
	iterator (std::vector<unsigned long>::iterator word, uint8 position) : _ref (word, position) {}
	iterator (const iterator &i) : std::_Bit_iterator(),_ref (i._ref._word, i._ref._pos) {}

	iterator &operator = (const iterator &i) {
		_ref._word = i._ref._word;
		_ref._pos = i._ref._pos;
		return *this;
	}

	iterator &operator ++ () 
	{
		if (++_ref._pos > __LINBOX_BITSOF_LONG_MUN) {
			++_ref._word;
			_ref._pos = 0UL;
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
		std::vector<unsigned long>::iterator new_word = _ref._word + (i >> __LINBOX_LOGOF_SIZE);
		uint8 new_pos = _ref._pos + (i & __LINBOX_POS_ALL_ONES);

		new_word += new_pos >> __LINBOX_LOGOF_SIZE;
		new_pos &= __LINBOX_POS_ALL_ONES;

#ifndef __CUDACC__
		return iterator (new_word, new_pos);
#else
		std::cout << "nvcc not happy here" << std::endl;
		exit(-4);
#endif
	}

	iterator &operator += (difference_type i) 
	{
		_ref._word += i >> __LINBOX_LOGOF_SIZE;
		_ref._pos  += i & __LINBOX_POS_ALL_ONES;
		_ref._word += _ref._pos >> __LINBOX_LOGOF_SIZE;
		_ref._pos  &= __LINBOX_POS_ALL_ONES;
		return *this;
	}

	iterator &operator -- () 
	{
		if (--_ref._pos > __LINBOX_BITSOF_LONG_MUN) {
			--_ref._word;
			_ref._pos = __LINBOX_BITSOF_LONG_MUN;
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
		{ return (_ref._word - i._ref._word) * __LINBOX_BITSOF_LONG + (_ref._pos - i._ref._pos); }

	reference operator [] (long i) 
		{ return *(*this + i); }

	reference operator * () 
		{ return _ref; }

	const_reference operator * () const 
		{ return _ref; }

	bool operator == (const iterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator != (const iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

    protected:
	friend class const_iterator;

	reference _ref;
};
 
class BitVector::const_iterator : public std::iterator <std::random_access_iterator_tag, bool>
{
    public:

	typedef std::iterator_traits<const_iterator>::iterator_category iterator_category;
	typedef std::iterator_traits<const_iterator>::reference const_reference;
	typedef std::iterator_traits<const_iterator>::pointer pointer;
	typedef std::iterator_traits<const_iterator>::value_type value_type;
	typedef std::iterator_traits<const_iterator>::difference_type difference_type;
	typedef BitVector::iterator iterator;

	const_iterator () : _ref (std::vector<unsigned long>::const_iterator (), 0UL) {}
	const_iterator (std::vector<unsigned long>::const_iterator word, uint8 position) : _ref (word, position) {}
	const_iterator (const const_iterator &i) : _ref (i._ref._word, i._ref._pos) {}

	const_iterator &operator = (const const_iterator &i) {
		this->_ref._word = i._ref._word;
		this->_ref._pos = i._ref._pos;
		return *this;
	}

	const_iterator &operator = (const iterator &i) {
		this->_ref._word = i._ref._word;
		this->_ref._pos = i._ref._pos;
		return *this;
	}

	const_iterator &operator ++ () 
	{
		if (++_ref._pos > __LINBOX_BITSOF_LONG_MUN) {
			++_ref._word;
			_ref._pos = 0UL;
		}

		return *this;
	}

	const_iterator operator ++ (int) 
	{
		const_iterator tmp (*this);
		++*this;
		return tmp;
	}

	const_iterator operator + (long i) const
	{
		std::vector<unsigned long>::const_iterator new_word = _ref._word + (i >> __LINBOX_LOGOF_SIZE);
		uint8 new_pos = _ref._pos + (i & __LINBOX_POS_ALL_ONES);

		new_word += new_pos >> __LINBOX_LOGOF_SIZE;
		new_pos &= __LINBOX_POS_ALL_ONES;

		return const_iterator (new_word, new_pos);
	}

	const_iterator &operator += (long i) 
	{
		_ref._word += i >> __LINBOX_LOGOF_SIZE;
		_ref._pos  += i & __LINBOX_POS_ALL_ONES;
		_ref._word += _ref._pos >> __LINBOX_LOGOF_SIZE;
		_ref._pos  &= __LINBOX_POS_ALL_ONES;
		return *this;
	}

	const_iterator &operator -- () 
	{
		if (--_ref._pos > __LINBOX_BITSOF_LONG_MUN) {
			--_ref._word;
			_ref._pos = __LINBOX_BITSOF_LONG_MUN;
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
		{ return (_ref._word - i._ref._word) * __LINBOX_BITSOF_LONG + (_ref._pos - i._ref._pos); }

	const_reference operator [] (difference_type i) const
		{ return *(*this + i); }

	const_reference operator * () const
		{ return _ref; }

	bool operator == (const const_iterator &c) const 
		{ return (_ref._word == c._ref._word) && (_ref._pos == c._ref._pos); }

	bool operator == (const BitVector::iterator &c) const 
		{ return (this->_ref._word == c._ref._word) && (this->_ref._pos == c._ref._pos); }

	bool operator != (const const_iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

	bool operator != (const BitVector::iterator &c) const 
		{ return (_ref._word != c._ref._word) || (_ref._pos != c._ref._pos); }

    protected:

	const_reference _ref;
};

inline BitVector::iterator BitVector::begin (void)
	{ return iterator (_v.begin (), 0UL); }

inline BitVector::const_iterator BitVector::begin (void) const
	{ return const_iterator (_v.begin (), 0UL); }

inline BitVector::iterator BitVector::end (void)
{
	if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return iterator (_v.end (), 0UL);
	else
		return iterator (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

inline BitVector::const_iterator BitVector::end (void) const
{
	if ((_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return const_iterator (_v.end (), 0UL);
	else
		return const_iterator (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

inline BitVector::reverse_iterator BitVector::rbegin (void)
	{ return reverse_iterator (end () - 1UL); }

inline BitVector::const_reverse_iterator BitVector::rbegin (void) const
	{ return const_reverse_iterator (end () - 1UL); }

inline BitVector::reverse_iterator BitVector::rend (void)
	{ return reverse_iterator (begin () - 1UL); }

inline BitVector::const_reverse_iterator BitVector::rend (void) const
	{ return const_reverse_iterator (begin () - 1UL); }

inline BitVector::reference BitVector::operator[] (BitVector::size_type n)
	{ return *(begin () + n); }

inline BitVector::const_reference BitVector::operator[] (BitVector::size_type n) const
	{ return *(begin () + n); }

inline BitVector::reference BitVector::at (BitVector::size_type n)
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

inline BitVector::const_reference BitVector::at (BitVector::size_type n) const
{
	if (n >= _size)
		throw std::out_of_range ("LinBox::BitVector");
	else
		return (*this)[n];
}

inline BitVector::reference BitVector::front (void)
	{ return reference (_v.begin (), 0UL); }

inline BitVector::const_reference BitVector::front (void) const
	{ return const_reference (_v.begin (), 0UL); }

inline BitVector::reference BitVector::back (void)
{
	if ( (_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return reference (_v.end (), 0UL);
	else
		return reference (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

inline BitVector::const_reference BitVector::back (void) const
{
	if ( (_size & __LINBOX_POS_ALL_ONES) == 0UL)
		return const_reference (_v.end (), 0UL);
	else
		return const_reference (_v.end () - 1UL, _size & __LINBOX_POS_ALL_ONES);
}

template<class Container>
inline BitVector &BitVector::operator = (const Container &v)
{
	typename Container::const_iterator i;
	typename Container::const_iterator i_end = v.begin () + (v.size () >> __LINBOX_LOGOF_SIZE);
	std::vector<unsigned long>::iterator j;
	unsigned long idx;

	_v.resize ((v.size () >> __LINBOX_LOGOF_SIZE) + ((v.size () & __LINBOX_POS_ALL_ONES) ? 1UL : 0UL));

	for (j = _v.begin (); i != i_end; ++j) {
		*j = 0UL;
		for (idx = 0UL; idx < __LINBOX_BITSOF_LONG; ++idx, ++i) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	if (v.size () & __LINBOX_POS_ALL_ONES) {
		*j = 0UL;

		for (idx = 0UL; idx < (v.size () & __LINBOX_POS_ALL_ONES); ++idx) {
			*j <<= 1UL;
			*j |= *i & 1UL;
		}
	}

	_size = v.size ();

	return *this;
}

inline void BitVector::resize (BitVector::size_type new_size, bool val)
	{_v.resize ((new_size >> __LINBOX_LOGOF_SIZE) + ((new_size & __LINBOX_POS_ALL_ONES) ? 1UL : 0UL), val ? __LINBOX_ALL_ONES : 0UL); _size = new_size; }

inline bool BitVector::operator == (const BitVector &v) const
{
	const_word_iterator i, j;
	unsigned long mask;

	if (_size != v._size) return false;

	for (i = wordBegin (), j = v.wordBegin (); i != wordEnd () - 1UL; ++i, ++j)
		if (*i != *j) return false;

	mask = (1UL << (_size & (8 * sizeof (unsigned long) - 1UL))) - 1UL;
	if (mask == 0UL) mask = (unsigned long) -1UL;

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

#endif // __LINBOX_bit_vector_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
