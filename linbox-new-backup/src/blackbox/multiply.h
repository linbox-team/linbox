/* File: src/library/objects/blackbox/multiply.h
 * Authro: William J Turner for the LinBox group
 */

#ifndef _MULTIPLY_
#define _MULTIPLY_

#include "LinBox/blackbox_archetype.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Blackbox multiply matrix.
   * This is a class that multiplies two matrices by implementing an 
   * apply method that calls the apply methods of both of the consituent 
   * matrices.
   *
   * This class, like the Black Box archetype from which it is derived, 
   * is templatized by the vector type to which the matrix is applied.  
   * Both constituent matrices must also use this same vector type.
   * @param Vector \Ref{LinBox} dense or sparse vector of field elements
   */
  template <class Vector> class multiply : public Blackbox_archetype<Vector>
  {
  public:

    /** Constructor from two black box matrices.
     * This constructor creates a matrix that is a product of two black box
     * matrices: A*B.
     * @param A_ptr pointer to black box matrix A.
     * @param B_ptr pointer to black box matrix B.
     */
    multiply(Blackbox_archetype<Vector>* A_ptr, 
	     Blackbox_archetype<Vector>* B_ptr)
    {
      // create new copies of matrices in dynamic memory
      if ( (A_ptr != 0) && (B_ptr != 0) 
	   && (A_ptr->coldim() == B_ptr->rowdim()) )
      {
	_A_ptr = A_ptr->clone();
	_B_ptr = B_ptr->clone();
      } // if ( (A_ptr != 0) && (B_ptr != 0) && (...) )
      else
	cerr << "ERROR: Cannot construct multiplication matrix." << endl;
      
    } // multiply(A, B)

    /** Copy constructor.
     * Creates new black box objects in dynamic memory.
     * @param M constant reference to multiply black box matrix
     */
    multiply(const multiply<Vector>& M)
    {
      // create new copies of matrices in dynamic memory
      if ( (M._A_ptr != 0) && (M._B_ptr != 0) )
      {
	_A_ptr = M._A_ptr->clone();
	_B_ptr = M._B_ptr->clone();
      } // if ( (M._A_ptr != 0) && (M._B_ptr != 0) )
      else
	cerr << "ERROR: Cannot (copy) construct multiplication matrix." << endl;
    } // multiply(const multiply<Vector>& M)

    /// Destructor
    ~multiply(void)
    {
      if (_A_ptr != 0) delete _A_ptr;
      if (_B_ptr != 0) delete _B_ptr;
    } // ~multiply(void)

    /** Virtual constructor.
     * Required because constructors cannot be virtual.
     * Make a copy of the Blackbox_archetype object.
     * Required by abstract base class.
     * @return pointer to new blackbox object
     */
    Blackbox_archetype<Vector>* clone() const
    { return new multiply(*this); }

    /** Application of BlackBox matrix.
     * y= (A*B)*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    inline Vector& apply(const Vector& x) const
    {
      if ( (_A_ptr != 0) && (_B_ptr != 0) )
      {
	Vector* y_ptr = new Vector;
	*y_ptr = _B_ptr->apply(x);
	*y_ptr = _A_ptr->apply(*y_ptr);
	return *y_ptr;
      } // if ( (_A_ptr != 0) && (_B_ptr != 0) )
      else
	return *(new Vector);
    } // Vector& apply(const Vector& x) const

    /** Application of BlackBox matrix transpose.
     * y= transpose(A*B)*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    inline Vector& applyTranspose(const Vector& x) const
    {
      if ( (_A_ptr != 0) && (_B_ptr != 0) )
      {
	Vector* y_ptr = new Vector;
	*y_ptr = _A_ptr->applyTranspose(x);
	*y_ptr = _B_ptr->applyTranspose(*y_ptr);
	return *y_ptr;
      } // if ( (_A_ptr != 0) && (_B_ptr != 0) )
      else
	return *(new Vector);
    } // Vector& applyTranspose(const Vector& x) const

    /** Retreive row dimensions of BlackBox matrix.
     * This may be needed for applying preconditioners.
     * Required by abstract base class.
     * @return integer number of rows of black box matrix.
     */
    size_t rowdim(void) const
    {
      if (_A_ptr != 0) 
	return _A_ptr->rowdim();
      else 
	return 0;
    } // size_t row(void) const
    
    /** Retreive column dimensions of BlackBox matrix.
     * Required by abstract base class.
     * @return integer number of columns of black box matrix.
     */
    size_t coldim(void) const 
    {
      if (_B_ptr != 0) 
	return _B_ptr->coldim();
      else 
	return 0;
    } // size_t coldim(void) const
      

  private:

    // Pointers to A and B matrices
    Blackbox_archetype<Vector>* _A_ptr;
    Blackbox_archetype<Vector>* _B_ptr;

  }; // template <Vector> class multiply

} // namespace LinBox

#endif // _MULTIPLY_
