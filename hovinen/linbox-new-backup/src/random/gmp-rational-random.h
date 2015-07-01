/* File: src/library/archetypes/random/ranom_archetype.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _GMP_RANDITER_ARCHETYPE_
#define _GMP_RANDITER_ARCHETYPE_

#include "LinBox/gmp-rational-number.h"
#include "LinBox/gmp-rational-field.h"

extern "C" {
#    include <sys/time.h>
#    include <stdlib.h>
}

// Namespace in which all LinBox code resides
namespace LinBox
{
  // forward declarations

  /** Random field element generator archetype.
   * Archetype for the random field element generator
   * common object interface to \Ref{LinBox}.
   *
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
  class GMP_Rational_Random
  {
  public:
    
    /** @name Common Object Interface.
     * These methods are required of all \Ref{LinBox} field element generators.
     */
    //@{
    
    /// Element type
    typedef GMP_Rational_Number element;
    
    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of zero means to use some arbitrary seed for the generator.
     * @param F LinBox field archetype object in which to do arithmetic
     * @param size constant integer reference of sample size from which to 
     *             sample (default = 0)
     * @param seed constant integer reference from which to seed random number 
     *             generator (default = 0)
     *
     * If size > 0, then set is {1,...,size}
     *
     * FIXME: We should really use a rational number for a seed here, but I'm
     * not entirely sure.
     */
    GMP_Rational_Random(const GMP_Rational_Field &F,
			const integer& size = 0,
			const integer& seed = 0)
	    : _F (F), _size (size), _seed (seed)
	 {
	      if (seed == 0)
		   _seed = time (NULL);
	 }

    /** Copy constructor.
     * Constructs GMP_Rational_Random object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * In this implementation, this means copying the random field element
     * generator to which R._randIter_ptr points.
     * @param  R GMP_Rational_Random object.
     */
    GMP_Rational_Random(const GMP_Rational_Random& R)
	 : _F (R._F), _size (R._size), _seed (R._seed) {}

    /** Destructor.
     * This destructs the random field element generator object.
     *
     * Vacuous
     */
    ~GMP_Rational_Random() 
    {}
    
    /** Assignment operator.
     * Assigns GMP_Rational_Random object R to generator.
     * In this implementation, this means copying the generator to
     * which R._randIter_ptr points.
     * @param  R GMP_Rational_Random object.
     */
    GMP_Rational_Random& operator=(const GMP_Rational_Random& R)
    {
      if (this != &R) // guard against self-assignment
      {
	   _F = R._F;
	   _seed = R._seed;
	   _size = R._size;
      }
      return *this;
    }
 
    /** Random field element creator.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return reference to random field element
     *
     * FIXME: How do we do this without causing serious memory leaks?
     */
    element& operator() (void) 
    {
	 element *new_elem;
	 unsigned int s;
	 int value;

	 if (_size == 0) {
	      new_elem = new element;
	      s = _seed;

	      mpz_set_si (mpq_numref (new_elem->rep), value);

	      do {
		   value = rand_r (&s);
	      } while (value == 0);

	      _seed = s;
	      mpz_set_si (mpq_denref (new_elem->rep), value);

	      return *new_elem;
	 }
	 else {
	      unsigned int s;
	      integer num, den;

	      s = _seed;
	      num = rand_r (&s);

	      if (_size > 0) {
		   num %= _size;
		   den = 1L;
	      } else {
		   den = rand_r (&s);
	      }

	      _seed = s;

	      new_elem = new element (num, den);
	      return *new_elem;
	 }
    }

    //@} Common Object Iterface
    
  private:

    /** Field in which we are doing arithmetic
     */
    GMP_Rational_Field _F;

    /** Size of set
     */
    integer _size;

    /** Seed for random number generation
     */
    integer _seed;
     
  }; // class GMP_Rational_Random
 
} // namespace LinBox

#endif // _GMP_RANDITER_ARCHETYPE_
