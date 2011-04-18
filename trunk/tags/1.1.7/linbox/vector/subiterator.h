/* linbox/vector/subiterator.h
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@acm.org>
 * Mods by -bds
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_subiterator_H
#define __LINBOX_subiterator_H

#include <iterator>
#include <vector>

// namespace in which all LinBox code resides
namespace LinBox
{

/** \brief Subvector iterator class provides striding iterators.
\ingroup vector

 *  A Subiterator steps by a fixed stride thru the underlying container.
 *  Subiter<Iterator> requires that Iterator be a random access iterator class 
 *  and then itself provides the full functionality of a random access iterator
 *  class.  See STL documentation for that functionality.  
 *  Documented here is only the constructor from (1) an iterator of an 
 *  underlying container and (2) a stride amount. 
 */
template <typename Iterator>
class Subiterator
{
    public:
	// Types
    	
	typedef typename std::iterator_traits<Iterator>::iterator_category	iterator_category;
	typedef typename std::iterator_traits<Iterator>::value_type		value_type;
	typedef typename std::iterator_traits<Iterator>::difference_type	difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer		pointer;
	typedef typename std::iterator_traits<Iterator>::reference		reference;

	// Constructors

	Subiterator () {}

	/** Subiterator p (pp, 3) provides an iterator which initially  has
	 *  the same reference, but for which increments and offsets step by
	 *  the amount stride rather than 1.  
	 *  Thus p+k is equivalent to pp+(stride*k).
	 * 
	 *  Striding iterators are easily positioned beyond the bounds of the 
	 *  underlying container.  It is up to the user to dereference the 
	 *  iterator only when it has a valid reference.
	 */
	Subiterator (const Iterator &iter, const difference_type& stride = 1)
		: _iter (iter), _stride (stride) {}

	template<class Iterator2>
	Subiterator (const Subiterator<Iterator2>& iter)
		: _iter (iter._iter), _stride (iter._stride) {}

	
	template<class Iterator2>
	Subiterator& operator = (const Subiterator<Iterator2>& sub)
	{
		_iter=sub._iter;
		_stride=sub._stride;
		return *this;
	}

	// Access operations
	reference operator * () const 
		{ return *_iter; }

	Iterator operator -> () const 
		{ return _iter; }

	reference operator [] (difference_type n) const 
		{ return _iter[n * _stride]; }

	// Iteration operations
    	
	Subiterator& operator ++ () 
		{ _iter += _stride; return *this; }
    	
	Subiterator operator ++ (int) 
		{ Subiterator tmp = *this; _iter += _stride; return tmp; }
    	
	Subiterator& operator -- () 
		{ _iter -= _stride; return *this; }
    	
	Subiterator operator -- (int) 
		{ Subiterator tmp = *this; _iter -= _stride; return tmp; }

	Subiterator operator + (difference_type n) const 
		{ return Subiterator (_iter + (n * _stride), _stride); }
    	
	Subiterator& operator += (difference_type n) 
		{ _iter += (n * _stride); return *this; }
    	
	Subiterator operator - (difference_type n) const 
		{ return Subiterator (_iter - (n * _stride), _stride); }
    	
	difference_type operator - (const Subiterator& x) const 
		{ return (_iter - x._iter)/_stride; }
    	
	Subiterator& operator -= (difference_type n) 
		{ _iter -= (n * _stride); return *this; }

	// Comparison operations

	bool operator == (const Subiterator& i) const 
		{ return ( (this->_stride == i._stride) && (this->_iter == i._iter) ); }
    	
	bool operator != (const Subiterator& i) const 
		{ return !(*this == i); }
    	
	bool operator < (const Subiterator& i) const 
		{ return (this->_iter < i._iter); }
    	
	bool operator > (const Subiterator& i) const 
		{ return (this->_iter > i._iter); }
    	
	bool operator <= (const Subiterator& i) const 
		{ return (this->_iter <= i._iter); }
    	
	bool operator >= (const Subiterator& i) const 
		{ return (this->_iter >= i._iter); }

	void swap (Subiterator& x)
		{ std::swap (_iter, x._iter); std::swap (_stride, x._stride); }
	
    protected:

	Iterator	_iter;		// wrapped iterator
	difference_type	_stride;	// length between iterations

}; // template <class Iterator> class Subiterator

} // namespace LinBox

#endif // __LINBOX_subiterator_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
