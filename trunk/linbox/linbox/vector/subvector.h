/* linbox/vector/subvector.h
 * Copyright (C) 2002 William J. Turner
 * See COPYING for license information.
 *
 * Written by William J. Turner <wjturner@acm.org>
 * Mods by -bds 
 * Maintainer: -bds 
 * (where there is missing or buggy function, please contact me rather than workaround)
 */

#include "subiterator.h"

// namespace in which all LinBox code resides
namespace LinBox
{

	/** Dense subvector class
	 * This class provides a statically sized subvector of a 
	 * random access container (such as std::vector, deque).
	 * It does not work on sparse linbox vectors.
	 * It implements all of the types and methods of a std::vector
	 * except for those that invalidate iterators, i.e.,
	 * those (potentially) involving vector resizing, such as
	 * push_back(), insert(), resize().
	 */
//	template <class Vector>
//	template <typename Iterator>//= typename Vector::iterator > 
	template <class Vector, typename Iterator > 
	class Subvector //: public Vector // for types
	{
	public:
		// Types

		typedef typename Vector::value_type	value_type;
		// should allocator_type even be offered?
		typedef typename Vector::allocator_type	allocator_type;
		typedef typename Vector::size_type	size_type;
		typedef typename Vector::difference_type	difference_type;

		typedef Iterator iterator;
		typedef const iterator 	const_iterator;
		typedef std::reverse_iterator<iterator>	reverse_iterator;
		typedef std::reverse_iterator<const_iterator>	const_reverse_iterator;

		typedef typename Vector::pointer	pointer;
		typedef typename Vector::const_pointer	const_pointer;
		typedef typename Vector::reference	reference;
		typedef typename Vector::const_reference	const_reference;

		// Constructors.   ... which should be explicit?

		Subvector()
		: _begin(0), _end(0) {}

		Subvector(Vector& v, size_type start, size_type stride, size_type length)
		: _begin(iterator (v.begin() + start, stride) ),
		 _end(iterator (v.begin() + start + (stride*length), stride) )
		{}
		
		Subvector(iterator begin, iterator end)
		: _begin(begin), _end(end) {}
		
		Subvector(iterator begin, size_type length)
		: _begin(begin), _end(begin + length) {}
		
		Subvector(const Subvector& x) 
		: _begin(x._begin), _end(x._end) {}

		~Subvector() {}

		// Iterators

		iterator begin(void)			{ return _begin; }
		const_iterator begin(void) const	{ return _begin; }
		iterator end(void)			{ return _end; }
		const_iterator end(void) const		{ return _end; }

		reverse_iterator rbegin(void)	{return reverse_iterator(_end);}
		const_reverse_iterator rbegin(void) const {return reverse_iterator(_end);}
		reverse_iterator rend(void) 	{return reverse_iterator(_begin);}
		const_reverse_iterator rend(void) const {return reverse_iterator(_begin);}

		// Element access

		reference operator[] (size_type n) 	{ return _begin[n]; }
		const_reference operator[] (size_type n) const { return _begin[n]; }

		// the method "at" does appear to be implemented 
		// in the gnu implementation of the STL
		reference at(size_type n)  // validity is relative to valid _begin, _end
		{   iterator p = _begin + n;
		    if ( _begin <= p && p < _end ) return *p;
		    else /*fixme: throw error; */ return *p;
		}

		const reference at(size_type n) const 
		{   const_iterator p = _begin + n;
		    if ( _begin <= p && p < _end ) return *p;
		    else /*fixme: throw error; */ return *p;
		}

		reference front(void) { return *_begin; }
		const_reference front(void) const { return *_begin; }
		reference back(void) { return *( _end - 1 ); }
		const_reference back(void) const { return *( _end - 1 ); }

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

//		template <class In> void assign(In first, In last);
//		void assign(size_type n, const T& val);

		// Stack operations:  
		// 	not implemented because they invalidate iterators

		// List operations:
		// 	not implemented because they invalidate iterators

		// Capacity
		// 	resize, reserve: not implemented because they 
		// 		invalidate iterators

		size_type size(void) const { return _end - _begin; }
		bool empty(void) const { return _end == _begin; }
		size_type max_size(void) const { return _end - _begin; }
		//size_type capacity(void) const { return _end - _begin; }

		// Swap
		void swap(Subvector& x)	// does this invalidate iterators?
		{ _begin.swap(x._begin); swap(_end, x._end); }
		
	protected:

		iterator _begin; // a subiterator of wrapped vector
		iterator _end;	 // a subiterator of wrapped vector

	}; // template <class Vector> class Subvector

	// Vector traits for Subvector wrapper
	template <class Vector, typename Iterator> 
	struct VectorTraits< Subvector<Vector, Iterator> >
       	{ 
		typedef typename VectorTraits<Vector>::VectorCategory VectorCategory; 
	};


} // namespace LinBox

