/* File: src/library/objects/random/param_modular_randiter.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _PARAM_MODULAR_RANDITER_
#define _PARAM_MODULAR_RANDITER_

#include <iostream>
#include <vector>
#include "LinBox/integer.h"
#include "LinBox/param_modular.h"

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
  class param_modular_randIter
  {
  public:

    /// Element type
    typedef integer element;

    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of zero means to use some arbitrary seed for the generator.
     * Purely virtual.
     * @param F LinBox field archetype object in which to do arithmetic
     * @param size constant integer reference of sample size from which to 
     *             sample (default = modulus of field)
     * @param seed constant integer reference from which to seed random number
     *             generator (default = 0)
     */
    param_modular_randIter(const param_modular& F, 
			   const integer& size = 0, 
			   const integer& seed = 0)
      : _F(F), _size(size), _seed(seed)
    { 
      if (_seed == 0) _seed = time(NULL);    

      integer cardinality = F.cardinality();
      if ( (_size == 0) 
	   || ( (cardinality != integer(-1)) && (_size > cardinality) ) )
	_size = cardinality;
    } // param_modular_randIter(const param_modular&, ...)

    /** Copy constructor.
     * Constructs param_modular_randIter object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * @param  R param_modular_randIter object.
     */
    param_modular_randIter(const param_modular_randIter& R) 
      : _F(R._F), _size(R._size), _seed(R._seed) {}

    /** Destructor.
     * This destructs the random field element generator object.
     */
    ~param_modular_randIter() {}
    
    /** Assignment operator.
     * Assigns param_modular_randIter object R to generator.
     * @param  R param_modular_randIter object.
     */
    param_modular_randIter& operator=(const param_modular_randIter& R)
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
      long temp_long;
      temp_long = static_cast<long>((double(rand())/RAND_MAX)*_size);
      temp_long %= _F.cardinality();
      if (temp_long < 0) temp_long += _F.cardinality();
      return *(new element(temp_long));
    } // element& operator() (void)

  private:

    /// Field in which arithmetic is done
    param_modular _F;

    /// Sampling size
    integer _size;
    
    /// Seed
    integer _seed;

  }; // class param_modular_randIter

} // namespace LinBox 

#endif // _PARAM_MODULAR_RANDITER_
