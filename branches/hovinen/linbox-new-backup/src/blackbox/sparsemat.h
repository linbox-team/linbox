/* File: sparsemat.h
 * Author: William J. Turner for the LinBox group
 */

#ifndef _SPARSEMAT_
#define _SPARSEMAT_

#include "LinBox/sparsemat_aux.h"
#include "LinBox/blackbox_archetype.h"
#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Blackbox sparse matrix template.
   * This is a class of sparse matrices templatized by the 
   * {@link Fields field} in
   * which the elements reside.  The matrix itself is stored as an
   * STL vector of \Ref{LinBox} sparse vectors of integers and field elements.
   * Each sparse vector corresponds to one row of the matrix, 
   * and each pair (j, a) in sparse vector i corresponds to the (i,j) 
   * entry of the matrix.
   *
   * The class is conforms to the 
   * {@link Archetypes archetype} for \Ref{Blackbox Matrices}.
   *
   * It has two base classes.  One, Blackbox_archetype, is an abstract base
   * class which ensures it adheres to the common object interface.  The
   * other, sparsemat, contains the mathematical functions and objects
   * necessary to implement the matrix.
   *
   * @param Field \Ref{LinBox} field
   * @param Row	  \Ref{LinBox} sparse vector implementation for rows of matrix
   * @param Vector \Ref{LinBox} dense or sparse vector of field elements
   */
  template <class Field, class Row, class Vector>
  class sparsemat
    : public Blackbox_archetype<Vector>, 
      public sparsemat_aux<Field, Row, Vector>
  {
  public:
 
    /** Constructor from sparsemat_aux<Field, Row, Vector>.
     * @param A constant reference to sparsemat object
     */
    sparsemat(const sparsemat_aux<Field, Row, Vector>& A) 
      : sparsemat_aux<Field, Row, Vector>(A) {}

    /** Constructor.
      * Note: the copy constructor and operator= will work as intended
      *       because of STL's container design
      * @param  F  the field of entries; passed so that a possible paramter 
      *            such as a modulus is known to the matrix.
      * @param  m  row dimension
      * @param  n  column dimension
      */
    sparsemat(const Field& F, size_t m, size_t n)
            : sparsemat_aux<Field, Row, Vector>(F, m, n) {}
    
    /** Virtual constructor.
     * Required because constructors cannot be virtual.
     * Make a copy of the Blackbox_archetype object.
     * Required by abstract base class.
     * @return pointer to new blackbox object
     */
    Blackbox_archetype<Vector>* clone() const 
    { return new sparsemat(*this); }

    /** Application of BlackBox matrix.
     * y= A*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    Vector& apply(const Vector& x) const
    { return sparsemat_aux<Field, Row, Vector>::apply(x); }

    /** Application of BlackBox matrix transpose.
     * y= transpose(A)*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    Vector& applyTranspose(const Vector& x) const
    { return sparsemat_aux<Field, Row, Vector>::applyTranspose(x); }

    /** Retreive row dimensions of BlackBox matrix.
     * This may be needed for applying preconditioners.
     * Required by abstract base class.
     * @return integer number of rows of black box matrix.
     */
    size_t rowdim(void) const 
    { return sparsemat_aux<Field, Row, Vector>::get_rowdim(); } 
    
    /** Retreive column dimensions of BlackBox matrix.
     * Required by abstract base class.
     * @return integer number of columns of black box matrix.
     */
    size_t coldim(void) const 
    { return sparsemat_aux<Field, Row, Vector>::get_coldim(); } 
  
  }; // sparsemat<Field> 

} // namespace LinBox

#endif // _SPARSEMAT_
