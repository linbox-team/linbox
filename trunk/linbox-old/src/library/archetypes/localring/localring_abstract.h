/* File: src/library/archetypes/localring/localring.h
 * Author: William Turner for the LinBox group
 */

#ifndef _LOCALRING_ABSTRACT_
#define _LOCALRING_ABSTRACT_

#include <iostream>
#include "LinBox/element_abstract.h"
#include "LinBox/randiter_abstract.h"
#include "LinBox/integer.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

  /** Abstract localring base class.
   * Found in the file \URL{src/library/archetypes/localring/localring_abstract.h}.
   * Abstract base class used to implement the localring archetype to minimize
   * code bloat.  All public member functions of this class are purely
   * virtual and must be implemented by all derived classes.
   *
   * This modifies the Field_abstract class by adding the divides and isUnit
   * member functions and modifying the spec of div members.
   *
   * If a template is instantiated on the localring archetype, we can change the
   * localring it is using by changing the derived class of this class.  This allows
   * us to change the localring used in a template without having to reinstantiate
   * it.  This minimizes code bloat, but it also introduces indirection through
   * the use of pointers and virtual functions which is inefficient.
   */
  class Localring_abstract
  {
  public:

    /// Element type.
    typedef Element_abstract element;

    /// Random iterator generator type.
    typedef RandIter_abstract randIter;
 
    /** @name Object Management
     * There are no public constructors for this class.
     * It should only be used in tandem with \Ref{Localring_archetype}.
     */
    //@{

    /** Destructor.
     * Required because of virtual member functions.
     * Virtual.
     */
    virtual ~Localring_abstract(void) {}

    /** Virtual copy constructor.
     * Required because constructors cannot be virtual.
     * Passes construction on to derived classes.
     * Purely virtual.
     * This function is not part of the common object interface.
     * @return pointer to new object in dynamic memory.
     */
    virtual Localring_abstract* clone() const = 0;

    /** Assignment operator.
     * Purely virtual.
     * @return reference to self
     * @param F constant reference to Localring_abstract object
     */
    virtual Localring_abstract& operator= (const Localring_abstract& F) = 0;

    /** Initialization of Localring element from an integer.
     * Behaves like C++ allocator construct.
     * This function assumes the output localring element x has already been
     * constructed, but that it is not already initialized.
     * Purely virtual.
     * @return reference to localring element.
     * @param x localring element to contain output (reference returned).
     * @param y integer.
     */
    virtual element& init(element& x, const integer& y) const = 0;
 
   /** Conversion of localring element to an integer.
     * This function assumes the output localring element x has already been
     * constructed, but that it is not already initialized.
     * Purely virtual.
     * @return reference to integer.
     * @param x reference to interger to contain output (reference returned).
     * @param y constant localring element.
     */
    virtual integer& convert(integer& x, const element& y) const = 0;
 
    /** Assignment of one localring element to another.
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& assign(element& x, const element& y) const = 0;

    /** Cardinality.
     * Return integer representing cardinality of the domain.
     * Returns a non-negative integer for all domains with finite
     * cardinality, and returns -1 to signify a domain of infinite
     * cardinality.
     * Purely virtual.
     * @return integer representing cardinality of the domain
     */
    virtual const integer& cardinality(void) const = 0;
 
    /** Characteristic.
     * Return integer representing characteristic of the domain.
     * Returns a positive integer to all domains with finite characteristic,
     * and returns 0 to signify a domain of infinite characteristic.
     * Purely virtual.
     * @return integer representing characteristic of the domain.
     */
    virtual const integer& characteristic(void) const = 0;

    //@} Object Management

    /** @name Arithmetic Operations
     * x <- y op z; x <- op y
     * These operations require all elements, including x, to be initialized
     * before the operation is called.  Uninitialized localring elements will
     * give undefined results.
     */
    //@{

    /** Equality of two elements.
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return boolean true if equal, false if not.
     * @param  x localring element
     * @param  y localring element
     */
    virtual bool areEqual(const element& x, const element& y) const = 0;

    /** Addition.
     * x = y + z
     * This function assumes all the localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    virtual element& add(element& x,
  			 const element& y, const element& z) const = 0;
 
    /** Subtraction.
     * x = y - z
     * This function assumes all the localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    virtual element& sub(element& x,
  			 const element& y, const element& z) const = 0;
 
    /** Multiplication.
     * x = y * z
     * This function assumes all the localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element.
     */
    virtual element& mul(element& x,
  			 const element& y, const element& z) const = 0;

    /** Divisibility
     * @return boolean true iff there exists ring element z such that xz = y.
     * @param  x localring element.
     * @param  y localring element.
     */
    bool divides(const element& x, const element& y) const = 0;

    /** Division.
     * x = y / z
     * Division is only defined for y,z such that divides(z, y).
     * This function assumes all the localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.  
     * Some x such that xz = y, and same x every time for given y,z.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     * @param  z localring element which divides y in the ring.
     */
    virtual element& div(element& x,
  			 const element& y, const element& z) const = 0;
 
    /** Additive Inverse (Negation).
     * x = - y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& neg(element& x, const element& y) const = 0;
 
    /** Multiplicative Inverse.
     * x = 1 / y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element which is a unit of the ring.
     */
    virtual element& inv(element& x, const element& y) const = 0;

    //@} Arithmetic Operations
 
    /** @name Inplace Arithmetic Operations
     * x <- x op y; x <- op x
     */
    //@{

    /** Zero equality.
     * Test if localring element is equal to zero.
     * This function assumes the localring element has already been
     * constructed and initialized.
     * Purely virtual.
     * @return boolean true if equals zero, false if not.
     * @param  x localring element.
     */
    virtual bool isZero(const element& x) const = 0;
 
    /** One equality.
     * Test if localring element is equal to one.
     * This function assumes the localring element has already been
     * constructed and initialized.
     * Purely virtual.
     * @return boolean true if equals one, false if not.
     * @param  x localring element.
     */
    virtual bool isOne(const element& x) const = 0;

    /** Inplace Addition.
     * x += y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& addin(element& x, const element& y) const = 0;
 
    /** Inplace Subtraction.
     * x -= y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& subin(element& x, const element& y) const = 0;
 
    /** Inplace Multiplication.
     * x *= y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& mulin(element& x, const element& y) const = 0;
 
    /** Inplace Division.
     * x /= y
     * This function assumes both localring elements have already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     * @param  y localring element.
     */
    virtual element& divin(element& x, const element& y) const = 0;
 
    /** Inplace Additive Inverse (Inplace Negation).
     * x = - x
     * This function assumes the localring element has already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     */
    virtual element& negin(element& x) const = 0;
 
    /** Inplace Multiplicative Inverse.
     * x = 1 / x
     * This function assumes the localring elementhas already been
     * constructed and initialized.
     * Purely virtual.
     * @return reference to x.
     * @param  x localring element (reference returned).
     */
    virtual element& invin(element& x) const = 0;

    //@} Inplace Arithmetic Operations

    /** @name Input/Output Operations */
    //@{

    /** Print localring.
     * Purely virtual.
     * @return output stream to which localring is written.
     * @param  os  output stream to which localring is written.
     */
    virtual ostream& write(ostream& os) const = 0;
 
    /** Read localring.
     * Purely virtual.
     * @return input stream from which localring is read.
     * @param  is  input stream from which localring is read.
     */
    virtual istream& read(istream& is) = 0;

    /** Print localring element.
     * This function assumes the localring element has already been
     * constructed and initialized.
     * Purely virtual.
     * @return output stream to which localring element is written.
     * @param  os  output stream to which localring element is written.
     * @param  x   localring element.
     */
    virtual ostream& write(ostream& os, const element& x) const = 0;
 
    /** Read localring element.
     * This function assumes the localring element has already been
     * constructed and initialized.
     * Purely virtual.
     * @return input stream from which localring element is read.
     * @param  is  input stream from which localring element is read.
     * @param  x   localring element.
     */
    virtual istream& read(istream& is, element& x) const = 0;

    //@}

  protected:

    /** Default Constructor.
     * Required by derived classes, but protected because this class should
     * never be constructed by itself.
     */
    Localring_abstract () {}
 
  private:

    /// Localring_archetype is friend.
    friend class Localring_archetype;


  }; // class Localring_abstract

} // namespace LinBox

#endif // _LOCALRING_ABSTRACT_
