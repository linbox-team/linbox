/* Copyright (C) 2010 LinBox
 * Written by William J Turner 
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_randiter_unparametric_H
#define __LINBOX_randiter_unparametric_H

#include <ctime>
#include <vector>

// Namespace in which all LinBox library code resides
namespace LinBox
{
	// forward declarations
	template <class K> class UnparametricField;
		
	/** Unparameterized random field element generator template.
	 * Implements LinBox random field element generator common object interface 
	 * for unparameterized fields.
	 * Used to generate efficient field classes for unparameterized fields.
	 * Constructs LinBox unparameterized random field element generators from 
	 * field types K.
	 * In particular, constructs LinBox random field element generators for
	 * unparameterized fields from field types that
	 * adhere to the operations for double, for
	 * example UnparametricRandIter< float >.
	 * Can be used as a pattern to write a particular
	 * field interface, such as, UnparametricRandIter< SaclibQ > as
	 * a template specialization.
	 * This implementation uses the standard C++ random number generator.  Thus,
	 * only one random field element generator can be used at a time since 
	 * creating a new one will re-seed the built-in generator and affect all 
	 * current LinBox generators.
	 * @param  K unparameterized field class
	 */
	template <class K> class UnparametricRandIter
	{
	public:
		
		/** @name Common Object Interface.
		 * These methods are required of all LinBox random field element generators.
		 */
		//@{
	 
		/** Field element type.
		 * The field element must contain a default constructor, 
		 * a copy constructor, a destructor, and an assignment operator.
		 */
		typedef K Element;    

		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min(size, F.cardinality(c)).
		 * A sampling size of zero means to sample from the entire field.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * This implementation sets the sampling size to be no more than the
		 * cardinality of the field.
		 * @param F LinBox field archetype object in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = 0)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		UnparametricRandIter(
				const UnparametricField<K>& F, 
				const integer& size = 0, 
				const integer& seed = 0
				)
			: _size(size), _seed(seed)
		{
			if (_seed == integer(0)) _seed = integer(time(NULL));
			
			integer cardinality; F.cardinality(cardinality);
			if ( (cardinality != integer(-1)) && (_size > cardinality) )
				_size = cardinality;

#ifdef TRACE
			cout << "created random generator with size " << _size 
	   << " and seed " << _seed << endl;
#endif // TRACE
			
			// Seed random number generator
			srand(static_cast<long>(_seed));

		} // UnparametricRandIter(const UnparametricField<K>&, const integer&, const integer&)

		/** Copy constructor.
		 * Constructs UnparametricRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * In this implementation, this means copying the random field element
		 * generator to which R._randIter_ptr points.
		 * @param  R UnparametricRandIter object.
		 */
		UnparametricRandIter(const UnparametricRandIter& R)
			: _size(R._size), _seed(R._seed) {}

		/** Destructor.
		 * This destructs the random field element generator object.
		 * In this implementation, this destroys the generator by deleting 
		 * the random generator object to which _randIter_ptr points.
		 */
		~UnparametricRandIter(void) {}
		
		/** Assignment operator.
		 * Assigns UnparametricRandIter object R to generator.
		 * In this implementation, this means copying the generator to
		 * which R._randIter_ptr points.
		 * @param  R UnparametricRandIter object.
		 */
		UnparametricRandIter& operator=(const UnparametricRandIter& R)
		{
			if (this != &R) // guard against self-assignment
			{
				_size = R._size;
				_seed = R._seed;
			}

			return *this;
		}
 
		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * @return random field element
		 */
		Element& random (Element& x) const
		{
			// Create new random elements
			if (_size == 0)
				return x = rand();
			else
				return x = static_cast<integer>((double(rand())/RAND_MAX)*double(_size));
		
		} // element& operator() (void)

		//@} Common Object Iterface
	 
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all 
		 * \Ref{LinBox Random field element generators}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/// Default constructor
		UnparametricRandIter(void) : _size(0), _seed(0) { time(NULL); }
		
		//@}

	private:

		/// Sampling size
		integer _size;
		
		/// Seed
		integer _seed;

	}; // template <class K> class UnparametricRandIter

} // namespace LinBox

#endif // __LINBOX_randiter_unparametric_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
