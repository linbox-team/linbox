/* File: sparsemat_map.h
 * Author: William J. Turner for the LinBox group
 */

#ifndef _SPARSEMAT_MAP_
#define _SPARSEMAT_MAP_

#include "LinBox/sparsemat_map_aux.h"
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
   * other, sparsemat_map, contains the mathematical functions and objects
   * necessary to implement the matrix.
   *
   * @param Field \Ref{LinBox} field
   * @param Row	  \Ref{LinBox} sparse vector implementation for rows of matrix
   * @param Vector \Ref{LinBox} dense or sparse vector of field elements
   */
  template <class Field, class Row, class Vector>
  class sparsemat_map
    : public Blackbox_archetype<Field, Vector>, 
      public sparsemat_map_aux<Field, Vector>
  {
  public:
 
    /** Constructor from sparsemat_map<Field>.
     * @param A constant reference to sparsemat_map object
     */
    sparsemat_map(const sparsemat_map_aux<Field, Vector>& A) 
      : sparsemat_map_aux<Field, Vector>(A) {}

    /** Virtual constructor.
     * Required because constructors cannot be virtual.
     * Make a copy of the Blackbox_archetype object.
     * Required by abstract base class.
     * @return pointer to new blackbox object
     */
    Blackbox_archetype<Field, Vector>* clone() const 
    { return new sparsemat_map(*this); }

    /** Application of BlackBox matrix.
     * y= A*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    Vector& apply(const Vector& x) const
    { return sparsemat_map_aux<Field, Vector>::apply(x); }

    /** Retreive row dimensions of BlackBox matrix.
     * This may be needed for applying preconditioners.
     * Required by abstract base class.
     * @return integer number of rows of black box matrix.
     */
    size_t rowdim(void) const 
    { return sparsemat_map_aux<Field, Vector>::get_rowdim(); } 
    
    /** Retreive column dimensions of BlackBox matrix.
     * Required by abstract base class.
     * @return integer number of columns of black box matrix.
     */
    size_t coldim(void) const 
    { return sparsemat_map_aux<Field, Vector>::get_coldim(); } 
  
  }; // sparsemat_map<Field> 

} // namespace LinBox

#endif // _SPARSEMAT_MAP_
