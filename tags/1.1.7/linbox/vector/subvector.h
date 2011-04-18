/* linbox/vector/subvector.h
 * Copyright (C) 2002 William J. Turner
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 * Written by William J. Turner <wjturner@acm.org>
 * Mods by -bds 
 * Maintained by: -bds <saunders@udel.edu> 
 * (where there is missing or buggy function, please contact me rather than workaround)
 */

#ifndef __LINBOX_subvector_H
#define __LINBOX_subvector_H

#include <linbox/vector/subiterator.h>
#include <iterator>
#include <linbox/vector/vector-traits.h>
#include <stdexcept>

namespace LinBox 	{
//wrapper Iterator to get a const Iterator

/** \brief Dense subvector 
\ingroup vector

 * This class provides a statically sized subvector of a 
 * random access container (such as std::vector, deque).
 * It does not work on sparse linbox vectors.
 * It implements all of the types and methods of a std::vector
 * except for those that invalidate iterators, i.e.,
 * those (potentially) involving vector resizing, such as
 * push_back(), insert(), resize().
 */
template <typename Iterator, typename ConstIterator = Iterator> 
class Subvector //: public Vector // for types
{
    public:
	// Types
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	// should allocator_type even be offered?
	//include <memory>
	//typedef allocator<value_type>	allocator_type;
	typedef typename VectorCategories::DenseVectorTag VectorCategory;
    
	typedef size_t                                              size_type;
	typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
	typedef typename std::iterator_traits<Iterator>::pointer	    pointer;
	typedef typename std::iterator_traits<Iterator>::reference	    reference;
	typedef const reference	                                    const_reference; 
	typedef Iterator                                            iterator;
	//typedef typename ConstIteratorType<Iterator>::const_iterator         const_iterator;

	typedef ConstIterator                                            const_iterator;

	typedef std::reverse_iterator<iterator>	                    reverse_iterator;
	typedef std::reverse_iterator<const_iterator>               const_reverse_iterator;
   	
	Subvector() : _begin(Iterator ()), _end(Iterator ()) {}
   	
	template<class Vector>
	Subvector (Vector& v, size_type start, size_type stride, size_type length)
	: _begin (iterator (v.begin() + start, stride)),
	  _end   (iterator (v.begin() + start + (stride * length), stride)) {}
	
	Subvector(iterator begin, iterator end) : _begin(begin), _end(end) {}
	
	Subvector(iterator begin, size_type length) : _begin(begin), _end(begin + length) {}
	
	//copy constructor
	Subvector(const Subvector& x) : _begin(x._begin), _end(x._end) {}
	
	~Subvector() {}
	
	// Iterators
	
	iterator               begin  (void)       { return _begin; }
	const_iterator         begin  (void) const { return _begin; }
	iterator               end    (void)       { return _end; }
	const_iterator         end    (void) const { return _end; }
	
	reverse_iterator       rbegin (void)       { return reverse_iterator (_end); }
	const_reverse_iterator rbegin (void) const { return reverse_iterator (_end); }
	reverse_iterator       rend   (void)       { return reverse_iterator (_begin); }
	const_reverse_iterator rend   (void) const { return reverse_iterator (_begin); }
	
	// Element access
	
	reference       operator[] (size_type n)       { return _begin[n]; }
	const_reference operator[] (size_type n) const { return _begin[n]; }
	
	// the method "at" does appear to be implemented 
	// in the gnu implementation of the STL
	reference at(size_type n)  // validity is relative to valid _begin, _end
	{   
		iterator p = _begin + n;
		if ( _begin <= p && p < _end ) 
			return *p;
		else 
			throw std::out_of_range("out of range"); //out of range error message.
	}
	
	const_reference at(size_type n) const 
	{
		const_iterator p = _begin + n;
		if ( _begin <= p && p < _end)
			return *p;
		else 
			throw std::out_of_range("out of range"); //out of range error message
	}
	
	reference       front (void)       { return *_begin; }
	const_reference front (void) const { return *_begin; }
	reference       back  (void)       { return *( _end - 1 ); }
	const_reference back  (void) const { return *( _end - 1 ); }
	
	template<class Container>
	/** assign the elements of Container one by one to *this.
	 *  Container must be at least as long as this.
	 */
	Subvector& operator= (const Container& x)
	{
		typename Container::const_iterator q = x.begin ();
		
		for (iterator p = _begin (); p != _end (); ++p, ++q)
			*p = *q;
		
		return *this;
	}
	
	template<class Iterator2, class ConstIterator2>
	Subvector& operator = (const Subvector<Iterator2, ConstIterator2>& sub)
	{
		_begin=sub.begin();
		_end=sub.end();
		return *this;
	}
	
	//		template <class In> void assign(In first, In last);
	//		void assign(size_type n, const T& val);
	
	// Stack operations:  
	// 	not implemented because they invalidate iterators
	
	// List operations:
	// 	not implemented because they invalidate iterators
	
	// Capacity
	// 	resize, reserve: not implemented because they 
	// 		invalidate iterators
	
	//copy assignment
	Subvector& operator=(const Subvector& sub)
	{ _begin = sub._begin; _end = sub._end; return *this; }
	
	size_type size      (void) const { return _end - _begin; }
	bool      empty     (void) const { return _end == _begin; }
	size_type max_size  (void) const { return _end - _begin; }
	//size_type capacity(void) const { return _end - _begin; }
	
	// Swap

	void swap (Subvector& x)
	{ std::swap (_begin,x._begin); std::swap(_end, x._end); }
	
	
    protected:
	
	iterator _begin; // a subiterator of wrapped vector
	iterator _end;	 // a subiterator of wrapped vector
	
}; // template <class Vector> class Subvector
  
// Vector traits for Subvector wrapper
template <typename Iterator, typename ConstIterator> 
struct VectorTraits<Subvector<Iterator, ConstIterator> >
{ 
	typedef VectorCategories::DenseVectorTag VectorCategory; 
};
  
  /*     Equality and unequality operators may be desirable, both for raw vectors of elements
	 and for vectors over a field with non-unique element representations.   In the latter
	 case a vector domain can provide it thru it's areEqual predicate.  The following
	 thought makes the valid point that raw vector equality could easily be misused.  
	 I remain agnostic on the subject.   -bds

 	 These and also < type operator comparisons are inappropriate for use
	 by linbox programmers, since we are interested in comparing vectors
	 of field elements not vectors of underlying representation type.
	 Not wishing to use these functions, we don't bother to implement them.
	 template<class Iterator>
	 bool operator==(const Subvector<Iterator>& sub1, const Subvector<Iterator>& sub2) const;
	 
	 template<class Iterator>
	 bool operator!=(const Subvector<Iterator>& sub1, const Subvector<Iterator>& sub2) const;
  */

} // namespace LinBox


namespace std {

	template<class TP>
       		void swap (TP&, TP&);


	template<class Iterator, class ConstIterator>
	void swap ( LinBox::Subvector<Iterator, ConstIterator>& x, 
		    LinBox::Subvector<Iterator, ConstIterator>& y ) {

		x. swap (y);

		
	}
}
	
#endif //__LINBOX_subvector_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
