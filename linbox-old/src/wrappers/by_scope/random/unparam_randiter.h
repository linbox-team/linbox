/* File: src/wrappers/by_scope/random/unparam_randiter.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _UNPARAM_RANDITER_WRAPPER_
#define _UNPARAM_RANDITER_WRAPPER_

#include <vector>

// Namespace in which all LinBox library code resides
namespace LinBox
{
  // forward declarations
  template <class K> class unparam_field;
    
  /** Unparameterized random field element generator template.
   * Implements LinBox random field element generator common object interface 
   * for unparameterized fields.
   * Used to generate efficient field classes for unparameterized fields.
   * Constructs LinBox unparameterized random field element generators from 
   * field types K.
   * In particular, constructs LinBox random field element generators for
   * unparameterized fields from field types that
   * adhere to the operations for double, for
   * example unparam_randIter< float >.
   * Can be used as a pattern to write a particular
   * field interface, such as, unparam_randIter< SaclibQ > as
   * a template specialization.
   * @param  K unparameterized field class
   */
  template <class K> class unparam_randIter
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
    typedef K element;    

    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
     * than 1/min(size, F.cardinality()).
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
    unparam_randIter(const unparam_field<K>& F, 
		     const integer& size = 0, 
		     const integer& seed = 0)
      : _size(size), _seed(seed), _loops(0)
    { 
      _randIter = _random.begin();
      
      if (_seed == integer(0)) _seed = time(NULL);
      
      integer cardinality = F.cardinality();
      if ( (cardinality != integer(-1)) && (_size > cardinality) )
	_size = cardinality;

    } // unparam_randIter(const unparam_field<K>&, const integer&, const integer&)

    /** Copy constructor.
     * Constructs unparam_randIter object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * In this implementation, this means copying the random field element
     * generator to which R._randIter_ptr points.
     * @param  R unparam_randIter object.
     */
    unparam_randIter(const unparam_randIter& R)
      : _size(R._size), _seed(R._seed), _random(R._random), _loops(R._loops)
    { _randIter = _random.begin() + (R._randIter - R._random.begin()); }

    /** Destructor.
     * This destructs the random field element generator object.
     * In this implementation, this destroys the generator by deleting 
     * the random generator object to which _randIter_ptr points.
     */
    ~unparam_randIter(void) {}
    
    /** Assignment operator.
     * Assigns unparam_randIter object R to generator.
     * In this implementation, this means copying the generator to
     * which R._randIter_ptr points.
     * @param  R unparam_randIter object.
     */
    unparam_randIter& operator=(const unparam_randIter& R)
    {
      if (this != &R) // guard against self-assignment
      {
	_size = R._size;
	_seed = R._seed;
	_random = R._random;
	_loops = R._loops;
      }

      _randIter = _random.begin() + (R._randIter - R._random.begin());

      return *this;
    }
 
    /** Random field element creator.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */
    element& operator() (void)
    {
      // If at end of vector, lengthen it
      if (_randIter == _random.end())
      {
	// Create new random vector
	_random = std::vector<K>(100, K());
	
	// Seed random number generator
	srand(_seed + _loops);

	// Create new random elements
	if (_size == 0)
	  for (_randIter = _random.begin(); 
	       _randIter != _random.end(); 
	       _randIter++)
	    *_randIter = rand();
	else
	  for (_randIter = _random.begin(); 
	       _randIter != _random.end(); 
	       _randIter++)
	    *_randIter = static_cast<long>((double(rand())/RAND_MAX)*_size);

	// Reset iterator, and update _loops
	_randIter = _random.begin();
	_loops++;
	
      } // if (_randIter == _random.end())

      return *(new K(*_randIter++));
      
    } // element& operator() (void)

    //@} Common Object Iterface
   
    /** @name Implementation-Specific Methods.
     * These methods are not required of all 
     * \Ref{LinBox Random field element generators}
     * and are included only for this implementation of the archetype.
     */
    //@{

    /// Default constructor
    unparam_randIter(void) : _size(0), _seed(0) 
    { 
      _randIter = _random.begin();
      if (_seed == integer(0)) _seed = time(NULL);    
    } // unparam_randIter(void)
    
    //@}

  private:

    /// Sampling size
    integer _size;
    
    /// Seed
    integer _seed;

    /// STL vector of random field elements
    std::vector<K> _random;

    /// STL vector iterator pointing to next random field element
    std::vector<K>::iterator _randIter;

    /// Number of times vector has been looped over; used to seed rand
    integer _loops;

  }; // template <class K> class unparam_randIter

} // namespace LinBox

#endif // _UNPARAM_RANDITER_WRAPPER_
