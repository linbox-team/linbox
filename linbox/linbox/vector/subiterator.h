/* linbox/vector/subiterator.h
 * Copyright (C) 2002 William J. Turner
 *
 * Written by William J. Turner <wjturner@acm.org>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

// namespace in which all LinBox code resides
namespace LinBox
{
	/** Subvector iterator class
	 */
	template <class Iterator> class Subiterator
	{
	public:
		// Types
		
		typedef iterator_traits<Iterator>::iterator_category	iterator_category;
		typedef iterator_traits<Iterator>::value_type		value_type;
		typedef iterator_traits<Iterator>::difference_type	difference_type;
		typedef iterator_traits<Iterator>::pointer		pointer;
		typedef iterator_traits<Iterator>::reference		reference;

		// Constructors

		Subiterator(const Iterator& iter, const difference_type& stride) 
			: _iter(iter), _stride(stride) {}

		// Access operations

		reference operator*() const { return *_iter; }
		pointer operator->() const { return &(operator*()); }
		reference operator[](difference_type n) const { return _iter[n * _stride]; }

		// Iteration operations
		
		Subiterator& operator++() 
		{
		       _iter += _stride; 
		       return *this; 
		}
		
		Subiterator operator++(int) 
		{ 
			Subiterator tmp = *this; 
			_iter += _stride; 
			return tmp; 
		}
		
		Subiterator& operator--() 
		{ 
			_iter -= _stride; 
			return *this; 
		}
		
		Subiterator operator--(int) 
		{ 
			Subiterator tmp = *this; 
			_iter -= _stride; 
			return tmp; 
		}

		Subiterator operator+(difference_type n) const 
		{ 
			return Subiterator(_iter + (n * _stride), _stride); 
		}
		
		Subiterator& operator+=(difference_type n) 
		{ 
			_iter += (n * _stride); 
			return *this; 
		}
		
		Subiterator operator-(difference_type n) const 
		{ 
			return Subiterator(_iter - (n * _stride), _stride); 
		}
		
		Subiterator& operator-=(difference_type n) 
		{ 
			_iter -= (n * _stride); 
			return *this; 
		}

		// Comparison operations

		bool operator==(const Subiterator& i) const 
		{ 
			return ( (this->_stride == i._stride) 
					&& (this->_iter == i._iter) ); 
		}
		
		bool operator!=(const Subiterator& i) const 
		{ 
			return !(*this == i); 
		}
		
		bool operator<(const Subiterator& i) const 
		{ 
			return (this->_iter < i._iter); 
		}
		
		bool operator>(const Subiterator& i) const 
		{ 
			return (this->_iter > i._iter); 
		}
		
		bool operator<=(const Subiterator& i) const 
		{ 
			return (this->_iter <= i._iter); 
		}
		
		bool operator>=(const Subiterator& i) const 
		{ 
			return (this->_iter >= i._iter); 
		}

	private:

		Iterator	_iter;		// wrapped iterator
		difference_type	_stride;	// length between iterations

	}; // template <class Iterator> class Subiterator

} // namespace LinBox

