/* File: src/library/archetypes/field/randiter_envelope.h
 * Author: William Turner for the LinBox group
 */

#ifndef _RANDITER_ENVELOPE_
#define _RANDITER_ENVELOPE_

#include <iostream>
#include "LinBox/field_envelope.h"
#include "LinBox/element_envelope.h"
#include "LinBox/randiter_abstract.h"

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
  template <class Field>
  class RandIter_envelope : public RandIter_abstract
  {
  public:

    /// Element type
    typedef Element_envelope<Field> element;

    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of size_t(-1) means to use some arbitrary seed for the generator.
     * For both of the last two parameters, a value of -1
     * @param F LinBox field envelope object in which to do arithmetic
     * @param size unsigned integer of sample size from which to sample
     *             (default = 0)
     * @param seed unsigned integer from which to seed random number generator
     *             (default = -1)
     */
    RandIter_envelope(const Field_envelope<Field>& F, 
		      size_t size = 0, 
		      size_t seed = size_t(-1))
      : _randIter(F._field, size, seed) {}

    /** Constructor from random field element generator to be wrapped
     * @param R random field element generator object to be wrapped
     */
    RandIter_envelope(const typename Field::randIter& R) : _randIter(R) {}

    /** Copy constructor.
     * Constructs RandIter_envelope object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * @param  R RandIter_envelope object.
     */
    RandIter_envelope(const RandIter_envelope& R) 
    { _randIter = R._randIter; }

    /** Destructor.
     * Required by abstract base class.
     * This destructs the random field element generator object.
     */
    ~RandIter_envelope() {}
    
    /** Assignment operator.
     * Assigns RandIter_envelope object R to generator.
     * Required by abstract base class.
     * @param  R RandIter_envelope object.
     */
    RandIter_abstract& operator=(const RandIter_abstract& R)
    {
      if (this != &R) // guard against self-assignment
	_randIter = static_cast<const RandIter_envelope&>(R)._randIter;

      return *this;
    }
 
    /** Virtual constructor from field, sampling size, and seed.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of size_t(-1) means to use some arbitrary seed for the generator.
     * Required by abstract base class.
     * @param F LinBox field abstract object in which to do arithmetic
     * @param size unsigned integer of sample size from which to sample
     *             (default = 0)
     * @param seed unsigned integer from which to seed random number generator
     *             (default = -1)
     */
    RandIter_abstract* construct(const Field_abstract& F, 
				 size_t size = 0, 
				 size_t seed = size_t(-1)) const
    { 
      return new 
	RandIter_envelope(static_cast<const Field_envelope<Field>&>(F)._field,
			  size,
			  seed);
    } // RandIter_abstract* construct(const Field_abstract&, size_t, size_t)

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * Required by abstract base class.
     * @return pointer to new RandIter_abstract object in dynamic memory.
     */
    RandIter_abstract* clone(void) const
    { return new RandIter_envelope(*this); }

    /** Random field element creator.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * Required by abstract base class.
     * @return reference to random field element
     */
    Element_abstract& operator() (void) 
    { return *(new Element_envelope<Field>(_randIter())); }

  private:

    typename Field::randIter _randIter;

  }; // class RandIter_envelope

} // namespace LinBox

#endif // _RANDITER_ENVELOPE_

