/* File: faxpy.h
 * Author: William J Turner for the LinBox group
 */
///////////  this needs adaptation to rings -bds ////////////

#ifndef _FAXPY_FIELD_OBJECT_
#define _FAXPY_FIELD_OBJECT_

// Namespace in which all LinBox library code resides
namespace LinBox
{
  /** Faxpy object.
   * This template class contains the faxpy object which is used for
   * the operation y = y + a*x or y += a*x.
   * This object is constructed from the field object F and a field 
   * element a which it stores and thus can use several times.
   * The use of an object instead of a static variable to store the element
   * a makes this method thread-safe.
   * @param Field \Ref{LinBox} {@link Fields field}
   */
  template< class Field > class faxpy 
  {
    
  public:

    /// Definition of element type
    typedef typename Field::element element;

    /** Constructor.
     * A faxpy object if constructed from a Field and a field element.
     * Copies of this objects are stored in the faxpy object.
     * @param F_init field F in which arithmetic is done
     * @param a_init field element a used for y += a*x operation
     */
    faxpy( const Field& F, const element& a ) : _F(F), _a(a) {}
   
    /** Apply method.
     * z = a*x + y.
     * @return reference to element z
     * @param z reference to element z
     * @param x constant reference to element x
     * @param y constant reference to element y
     */
    element& apply(element& z, const element& x, const element& y) const
    {
      _F.mul(z,_a,x);	// multiply a*x
      _F.addin(z,y);	// add y + a*x
      return z;
    }

    /** Inplace apply method.
     * y = a*x + y or y += a*x.
     * @return reference to element y
     * @param y reference to element y
     * @param x constant reference to element x
     */
    element& applyin(element& y, const element& x) const
    {
      element temp(x);
      _F.mul(temp,_a,x);	// multiply a*x
      _F.addin(y, temp);	// add y + a*x
      return y;
    }

    /** Assign method.
     * Stores new field element for arithmetic.
     * @return reference to self
     * @param a_init constant reference to element a
     */
    faxpy& assign(const element& a)
    {
      _a = a;
      return *this;
    }

  private:

    /// Field in which arithmetic is done
    Field _F;

    /// Field element for arithmetic
    element _a;

  }; // class faxpy

} // namespace LinBox

#endif // _FAXPY_FIELD_OBJECT_

