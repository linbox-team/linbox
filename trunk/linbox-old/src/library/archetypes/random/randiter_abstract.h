/* File: src/library/archetypes/field/field_abstract.h
 * Author: William Turner for the LinBox group
 */

#ifndef _RANDITER_ABSTRACT_
#define _RANDITER_ABSTRACT_
#include "LinBox/integer.h"

#include <iostream>

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
  // forward declarations
  class Field_abstract;
  class Element_abstract;

  /** Random field element generator.
   * This encapsulated class is a generator of random field elements for 
   * the encapsulating field.
   * It is required to contain constructors from a field object and
   * two integers.  The first integer being a cardinality of a set to 
   * draw the random elements from, and the second being a seed for the 
   * random number generator.
   * It is also required to contain a copy constructor, a destructor, and
   * an operator() which acts on a reference to a field element.  In this 
   * operator(), the random element is placed into the input field element 
   * and also returned as a reference.
   */
  class RandIter_abstract
  {
  public:
    
    /** Virtual constructor from field, sampling size, and seed.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of size_t(-1) means to use some arbitrary seed for the generator.
     * Purely virtual.
     * @param F LinBox field archetype object in which to do arithmetic
     * @param size unsigned integer of sample size from which to sample
     *             (default = 0)
     * @param seed unsigned integer from which to seed random number generator
     *             (default = -1)
     */
    virtual RandIter_abstract* construct(const Field_abstract& F, 
					 integer size = 0, 
					 size_t seed = size_t(-1)) const = 0;

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * Purely virtual.
     * @return pointer to new RandIter_abstract object in dynamic memory.
     */
    virtual RandIter_abstract* clone(void) const = 0;

    /** Assignment operator.
     * Purely virtual.
     * @param  x constant reference to RandIter_abstract object
     * @return reference to self
     */
    virtual RandIter_abstract& operator=(const RandIter_abstract& x) = 0;

    /** Destructor.
     */
    virtual ~RandIter_abstract(void) {}

    /** Random field element creator.
     * Purely virtual.
     * @return reference to Element_abstract object
     */
    virtual Element_abstract& operator() (void) = 0;

  protected:

    /** Default constructor
     * Required by derived classes, but protected because this class should
     * never be constructed by itself.
     */
    RandIter_abstract(void) {}

  }; // class RandIter_abstract

} // namespace LinBox

#endif // _RANDITER_ABSTRACT_
