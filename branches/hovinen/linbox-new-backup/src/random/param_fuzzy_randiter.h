/* File: src/library/objects/random/param_fuzzy_randiter.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _PARAM_FUZZY_RANDITER_
#define _PARAM_FUZZY_RANDITER_

#include <iostream>
#include <vector>
#include "LinBox/integer.h"
#include "LinBox/param_fuzzy.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  /** Random field base element generator.
   * This encapsulated class is a generator of random field base elements for 
   * the encapsulating field.
   * It is required to contain constructors from a field object and
   * two integers.  The first integer being a cardinality of a set to 
   * draw the random elements from, and the second being a seed for the 
   * random number generator.
   * It is also required to contain a copy constructor, a destructor, and
   * an operator() which acts on a reference to a field base element.  In this 
   * operator(), the random element is placed into the input field base element 
   * and also returned as a reference.
   */
  class param_fuzzy_randIter
  {
  public:

    /// Element type
    typedef double element;

    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of zero means to use some arbitrary seed for the generator.
     * Purely virtual.
     * @param F LinBox field archetype object in which to do arithmetic
     * @param size constant integer reference of sample size from which to 
     *             sample (default = 0)
     * @param seed constant integer reference from which to seed random number
     *             generator (default = 0)
     */
    param_fuzzy_randIter(const param_fuzzy& F, 
			     const integer& size = 0, 
			     const integer& seed = 0)
      : _F(F), _size(size), _seed(seed)
    { 
      if (_size == 0) _size = F.cardinality();
      if (_seed == 0) _seed = time(NULL);    
    } // param_fuzzy_randIter(const param_fuzzy&, ...)

    /** Copy constructor.
     * Constructs param_fuzzy_randIter object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * @param  R param_fuzzy_randIter object.
     */
    param_fuzzy_randIter(const param_fuzzy_randIter& R) 
      : _F(R._F), _size(R._size), _seed(R._seed) {}

    /** Destructor.
     * This destructs the random field element generator object.
     */
    ~param_fuzzy_randIter() {}
    
    /** Assignment operator.
     * Assigns param_fuzzy_randIter object R to generator.
     * @param  R param_fuzzy_randIter object.
     */
    param_fuzzy_randIter& operator=(const param_fuzzy_randIter& R)
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
     * Required by abstract base class.
     * @return reference to random field element
     */
    element& operator() (void) 
    {
      // Create new random elements
      if (_size == 0)
	return *(new element(rand()));
      else
	return *(new element(static_cast<long>((double(rand())/RAND_MAX)*_size)));
    } // element& operator() (void)

  private:

    /// Field in which arithmetic is done
    param_fuzzy _F;

    /// Sampling size
    integer _size;
    
    /// Seed
    integer _seed;

  }; // class param_fuzzy_randIter : public param_fuzzy_randIter

} // namespace LinBox 

#endif // _PARAM_FUZZY_RANDITER_
