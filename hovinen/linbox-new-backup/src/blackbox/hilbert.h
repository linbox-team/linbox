/* File: src/library/objects/blackbox/hilbert.h
 * Authro: William J Turner for the LinBox group
 */

#ifndef _HILBERT_
#define _HILBERT_

#include "LinBox/blackbox_archetype.h"
#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Blackbox Hilbert matrix.
   * This is a class of n by n Hilbert matrices templatized by the 
   * {@link Fields field} in 
   * which the elements reside.  The class conforms to the 
   * {@link Archetypes archetype} for \Ref{BlackBox Matrices}.
   *
   * The matrix itself is not stored in memory.  Rather, its apply
   * methods use a vector of {@link Fields field} elements, which are 
   * used to "multiply" the matrix to a vector.
   * 
   * This class has three template parameters.  The first is the field in 
   * which the arithmetic is to be done.  The second is the type of 
   * \Ref{LinBox} vector to which to apply the matrix.  The 
   * third is chosen be defualt to be the \Ref{LinBox} vector trait
   * of the vector.  This class is then specialized for dense and sparse 
   * vectors.
   * 
   * The default class is not implemented.  It's functions should never
   * be called because partial template specialization should always be
   * done on the vector traits.
   * @param Field \Ref{LinBox} field
   * @param Vector \Ref{LinBox} dense or sparse vector of field elements
   * @param Trait  Marker whether to use dense or sparse LinBox vector 
   *               implementation.  This is chosen by a default parameter 
   *               and partial template specialization.
   */
  template <class Field, class Vector,
	    class Trait = vector_traits<Vector>::vector_category>
  class hilbert : public Blackbox_archetype<Vector>
  {
  public:

    /// Element definition
    typedef typename Field::element element;

    /** Constructor from integer and field.
     * @param n size_t integer number of rows and columns of matrix.
     */
    hilbert(Field F, size_t n);

    /** Virtual constructor.
     * Required because constructors cannot be virtual.
     * Make a copy of the Blackbox_archetype object.
     * Required by abstract base class.
     * @return pointer to new blackbox object
     */
    Blackbox_archetype<Vector>* clone() const;

    /** Application of BlackBox matrix.
     * y= A*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    Vector& apply(const Vector& x) const;

    /** Application of BlackBox matrix transpose.
     * y= transpose(A)*x.
     * Requires one vector conforming to the \Ref{LinBox}
     * vector {@link Archetypes archetype}.
     * Required by abstract base class.
     * Because the Hilbert matrix is symmetric, this is the same as calling 
     * the apply function.
     * @return reference to vector y containing output.
     * @param  x constant reference to vector to contain input
     */
    Vector& applyTranspose(const Vector& x) const;

   /** Retreive row dimensions of BlackBox matrix.
     * This may be needed for applying preconditioners.
     * Required by abstract base class.
     * @return integer number of rows of black box matrix.
     */
    size_t rowdim(void) const;
    
    /** Retreive column dimensions of BlackBox matrix.
     * Required by abstract base class.
     * @return integer number of columns of black box matrix.
     */
    size_t coldim(void) const;

  }; // template <Field, Vector> class hilbert
 
  // Specialization of hilbert for LinBox dense vectors
  template <class Field, class Vector>
  class hilbert<Field, Vector, vector_categories::dense_vector_tag>
    : public Blackbox_archetype<Vector>
  {
  public:

    typedef typename Field::element element;
    hilbert(Field F, size_t n);
    Blackbox_archetype<Vector>* clone() const 
      { return new hilbert(*this); }
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const { return apply(x); }
    size_t rowdim(void) const { return _n; } 
    size_t coldim(void) const { return _n; } 

  private:

    // Field for arithmetic
    Field _F;

    // Number of rows and columns of square matrix.
    size_t _n;

    // STL vector of field elements used in applying matrix.
    std::vector<element> _H;
    
  }; // template <Field, Vector> class hilbert<dense_vector_tag>
   
  // Specialization of hilbert for LinBox sparse sequence vectors
  template <class Field, class Vector>
  class hilbert<Field, Vector, vector_categories::sparse_sequence_vector_tag>
    : public Blackbox_archetype<Vector>
  {
  public:

    typedef typename Field::element element;
    hilbert(Field F, size_t n);
    Blackbox_archetype<Vector>* clone() const 
      { return new hilbert(*this); }
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const { return apply(x); }
    size_t rowdim(void) const { return _n; } 
    size_t coldim(void) const { return _n; } 

  private:

    // Field for arithmetic
    Field _F;

    // Number of rows and columns of square matrix.
    size_t _n;

    // STL vector of field elements used in applying matrix.
    std::vector<element> _H;
    
  }; // template <Field, Vector> class hilbert<sparse_sequence_vector_tag>

  // Specialization of hilbert for LinBox sparse associative vectors
  template <class Field, class Vector>
  class hilbert<Field, Vector, vector_categories::sparse_associative_vector_tag>
    : public Blackbox_archetype<Vector>
  {
  public:

    typedef typename Field::element element;
    hilbert(Field F, size_t n);
    Blackbox_archetype<Vector>* clone() const 
      { return new hilbert(*this); }
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const { return apply(x); }
    size_t rowdim(void) const { return _n; } 
    size_t coldim(void) const { return _n; } 

  private:

    // Field for arithmetic
    Field _F;

    // Number of rows and columns of square matrix.
    size_t _n;

    // STL vector of field elements used in applying matrix.
    std::vector<element> _H;
    
  }; // template <Field, Vector> class hilbert<sparse_associative_vector_tag>

  // Method implementations for dense vectors
 
  template <class Field, class Vector>
  inline hilbert<Field, Vector, vector_categories::dense_vector_tag>
  ::hilbert(Field F, size_t n) : _F(F), _n(n)
  {
    element one, temp;
    _F.init(one, 1);
    _F.init(temp, 0);

    _H = std::vector<element>(2*_n - 1, temp);

    std::vector<element>::iterator iter;

    for (iter = _H.begin(); iter!= _H.end(); iter++)
    {
      _F.addin(temp, one);
      _F.div(*iter, one, temp);
    }
 
  } // hilbert<dense_vctor_tag>::hilbert(Field, size_t)

  template <class Field, class Vector>
  inline Vector& hilbert<Field, Vector, vector_categories::dense_vector_tag>
  ::apply(const Vector& x) const
  {
    // Create zero vector to hold output
    element temp;
    _F.init(temp, 0);
    Vector* y_ptr = new Vector(_n, temp);
 
    if (_n != x.size())
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl
  	   << endl;
      return *y_ptr;
    }
 
    // Create iterators for input, output, and stored vectors
    std::vector<element>::const_iterator iter, start_iter;
    typename Vector::const_iterator x_iter;
    typename Vector::iterator y_iter;
 
    // Start at beginning of _H vector for first row
    start_iter = _H.begin();
 
    // Iterator over elements of output vector.
    // For each element, multiply row of matrix with input vector.
    // Each row of matrix starts one further in _H vector.
    for (y_iter = y_ptr->begin();
  	 y_iter != y_ptr->end();
  	 y_iter++, start_iter++)
    {
      // start matrix row at correct place
      iter = start_iter;
 
      // Multiply matrix row and input vector by iterating over both.
      for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)
      {
  	_F.mul(temp, *iter, *x_iter);
  	_F.addin(*y_iter, temp);
      } // for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)
 
    } // for (y_iter = y.begin(); y_iter != y.end(); y_iter++, start_iter++)
 
    return *y_ptr;
  } // Vector& hilbert<dense_vector_tag>::apply(const Vector&) const
  
  // Method implementations for sparse sequence vectors
 
  template <class Field, class Vector>
  inline hilbert<Field, Vector, vector_categories::sparse_sequence_vector_tag>
  ::hilbert(Field F, size_t n) : _F(F), _n(n)
  {
    element one, temp;
    _F.init(one, 1);
    _F.init(temp, 0);

    _H = std::vector<element>(2*_n - 1, temp);

    std::vector<element>::iterator iter;

    for (iter = _H.begin(); iter!= _H.end(); iter++)
    {
      _F.addin(temp, one);
      _F.div(*iter, one, temp);
    }
 
  } // hilbert<sparse_sequence_vector_tag>::hilbert(Field, size_t)

  template <class Field, class Vector>
  inline Vector& 
  hilbert<Field, Vector, vector_categories::sparse_sequence_vector_tag>
  ::apply(const Vector& x) const
  {
    // Create zero vector to hold output
    Vector* y_ptr = new Vector();
 
    if ( (!x.empty()) && (_n < x.back().first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl
  	   << endl;
      return *y_ptr;
    }

    // create field elements to be used in calculations
    element zero, entry, temp;
    _F.init(zero, 0);
    _F.init(entry, 0);
    _F.init(temp, 0);

    // Create iterators for input, output, and stored vectors
    std::vector<element>::const_iterator iter, start_iter;
    typename Vector::const_iterator x_iter;
 
    // Start at beginning of _H vector for first row
    start_iter = _H.begin();
 
    // Iterator over indices of output vector.
    // For each element, multiply row of matrix with input vector,
    // and insert non-zero elements into vector
    // Each row of matrix starts one further in _H vector.
    for (size_t i = 0; i < _n; i++, start_iter++)
    {
      entry = zero;
 
      // Multiply matrix row and input vector by iterating over both.
      for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)
      {
//  	_F.mul(temp, *(start_iter + x_iter->first), x_iter->second); 
	  // problems with deque
  	
	_F.mul(temp, *(start_iter + (*x_iter).first), (*x_iter).second);
  	_F.addin(entry, temp);
      } // for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)

      if (!_F.isZero(entry)) y_ptr->push_back(make_pair(i, entry));

    } // for (size_t i = 0; i < _n; i++, start_iter++)

    return *y_ptr;

  } // Vector& hilbert<sparse_sequence_vector_tag>::apply(const Vector&) const

  // Method implementations for sparse associative vectors
 
  template <class Field, class Vector>
  inline hilbert<Field, Vector, vector_categories::sparse_associative_vector_tag>
  ::hilbert(Field F, size_t n) : _F(F), _n(n)
  {
    element one, temp;
    _F.init(one, 1);
    _F.init(temp, 0);

    _H = std::vector<element>(2*_n - 1, temp);

    std::vector<element>::iterator iter;

    for (iter = _H.begin(); iter!= _H.end(); iter++)
    {
      _F.addin(temp, one);
      _F.div(*iter, one, temp);
    }
 
  } // hilbert<sparse_associative_vector_tag>::hilbert(Field, size_t)

  template <class Field, class Vector>
  inline Vector& hilbert<Field, 
                         Vector, 
			 vector_categories::sparse_associative_vector_tag>
  ::apply(const Vector& x) const
  {
    // Create zero vector to hold output
    Vector* y_ptr = new Vector();
 
    if ( (!x.empty()) && (_n < x.rbegin()->first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl
  	   << endl;
      return *y_ptr;
    }

    // create field elements to be used in calculations
    element zero, entry, temp;
    _F.init(zero, 0);
    _F.init(entry, 0);
    _F.init(temp, 0);

    // Create iterators for input, output, and stored vectors
    std::vector<element>::const_iterator iter, start_iter;
    typename Vector::const_iterator x_iter;
 
    // Start at beginning of _H vector for first row
    start_iter = _H.begin();
 
    // Iterator over indices of output vector.
    // For each element, multiply row of matrix with input vector,
    // and insert non-zero elements into vector
    // Each row of matrix starts one further in _H vector.
    for (size_t i = 0; i < _n; i++, start_iter++)
    {
      entry = zero;
 
      // Multiply matrix row and input vector by iterating over both.
      for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)
      {
  	_F.mul(temp, *(start_iter + x_iter->first), x_iter->second); 
  	_F.addin(entry, temp);
      } // for (x_iter = x.begin(); x_iter != x.end(); x_iter++, iter++)

      if (!_F.isZero(entry)) y_ptr->insert(y_ptr->end(), make_pair(i, entry));

    } // for (size_t i = 0; i < _n; i++, start_iter++)

    return *y_ptr;

  } // Vector& hilbert<sparse_associative_vector_tag>::apply(...) const

} // namespace LinBox

#endif // _HILBERT_
