/* linbox/vector/subvector.h
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@acm.org>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "subiterator.h"

// namespace in which all LinBox code resides
namespace LinBox
{

	/** Dense subvector class
	 * This class accesses a subvector of a dense vector.
	 * It does not work on any sparse vectors.
	 * It implements all of the methods of an STL vector and its
	 * iterators except for those that invalidate iterators.
	 */
	template <class Vector> class Subvector
	{
	public:
		// Types

		typedef typename Vector::value_type	value_type;
		typedef typename Vector::allocator_type	allocator_type;
		typedef typename Vector::size_type	size_type;
		typedef typename Vector::difference_type	difference_type;

		typedef Subiterator<typename Vector::iterator>	iterator;
		typedef Subiterator<typename Vector::const_iterator>	const_iterator;
		typedef std::reverse_iterator<iterator>	reverse_iterator;
		typedef std::reverse_iterator<const_iterator>	const_reverse_iterator;

		typedef typename Vector::pointer	pointer;
		typedef typename Vector::const_pointer	const_pointer;
		typedef typename Vector::reference	reference;
		typedef typename Vector::const_reference	const_reference;

		// Iterators

		iterator begin(void) { return iterator(_v.begin() + _start, _stride); }
		const_iterator begin(void) const { return const_iterator(_v.begin() + _start, _stride); }
		iterator end(void) { return (this->begin() + _length); }
		const_iterator end(void) const { return (this->begin() + _length); }

#if 1
		reverse_iterator rbegin(void) { return reverse_iterator( this->end() ); }
		const_reverse_iterator rbegin(void) const { return reverse_iterator( this->end() ); }
		reverse_iterator rend(void) { return reverse_iterator( this->begin() ); }
		const_reverse_iterator rend(void) const { return reverse_iterator( this->begin() ); }
#endif

		// Element access

		reference operator[] (size_type n) { return _v[_start + (n * _stride)]; }
		const_reference operator[] (size_type n) const { return _v[_start + (n * _stride)]; }

#if 0	// the method "at" does appear to be implemented 
	// in the gnu implementation of the STL
		
		reference at(size_type n) { return _v.at(_start + (n * _stride)); };
		const reference at(size_type n) const { return _v.at(_start + (n * _stride)); };
#endif

		reference front(void) { return *( this->begin() ); }
		const_reference front(void) const { return *( this->begin() ); }
		reference back(void) { return *( this->end() - 1 ); }
		const_reference back(void) const { return *( this->end() - 1 ); }

		// Constructors, etc

		Subvector(){}

		Subvector(Vector& v, size_type start, size_type stride, size_type length)
			: _v(v), _start(start), _stride(stride), _length(length) {}

#if 0
/* I think we can't have any of these constructors because we want to avoid
constructors that copy the underlying vector.  -bds
*/
//		explicit Subvector(const A& = A());
//		explicit Subvector(size_type n, const T& val = T(), const A& = A());
		template <class In> Subvector(In first, In last, const allocator_type& A = allocator_type())
			: _v(first, last, A), _start(0), _stride(1), _length(1)
			{ _length = _v.size(); }
		
#endif
		// allow copy construction for subvectors of the same vector.
		Subvector(const Subvector& x) 
			: _v(x._v), _start(x._start), _stride(x._stride), _length(x._length) {}

		~Subvector() {}

                template<class Container>
		/** assign the elements of Container one by one to *this.
		 *  Container must be at least as long as this.
		 */
		Vector& operator=(const Container& x)
		{
			typename Container::const_iterator q = x.begin();
			for( iterator p = begin(); p != end(); ++p, ++q )
				*p = *q;	
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

		size_type size(void) const { return _length; }
		bool empty(void) const { return size() == 0; }
		size_type max_size(void) const { return _length; }
		size_type capacity(void) const { return _length; }

		// Swap

#if 0
		void swap(Subvector&);	// does this invalidate iterators?
#endif
		
//	private:

		Vector& _v;		// wrapped vector
		size_type _start;	// starting position
		size_type _stride;	// length between iterations
		size_type _length;	// length of subvector

	}; // template <class Vector> class Subvector

	// Helper functions
	// Vector traits for Subvector wrapper
	template <class Vector> struct VectorTraits< Subvector<Vector> >
       	{ 
		typedef typename VectorTraits<Vector>::VectorCategory VectorCategory; 
	};


} // namespace LinBox

