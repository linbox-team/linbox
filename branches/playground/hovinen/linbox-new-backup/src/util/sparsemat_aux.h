/* File: sparsemat_aux.h
 * Author: William J. Turner for the LinBox group
 */

#ifndef _SPARSEMAT_AUX_
#define _SPARSEMAT_AUX_

#include <vector>    // STL vector
#include <utility>   // STL pair
#include <iostream>
#include <algorithm>

#include "LinBox/sparsemat_base.h"
#include "LinBox/faxpy.h"
#include "LinBox/vector_traits.h"

#include "Examples/vector_utility.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

  /** Auxillary sparse matrix template.
   * This class is derived from \Ref{sparsemat_base} which stores the actual
   * matrix implementation.  This class adds the functions gauss, linsolve, 
   * and apply to the sparse matrix through two new template parameters.  The
   * second of these is chosen be default to be the \Ref{LinBox} vector trait
   * of the vector.  This class is then specialized for dense and sparse 
   * vectors.  The separation of these methods from the rest of the class is 
   * done to avoid implementing all methods twice.  Now, only the methods
   * that depend on the vectors are implemented three times.
   *
   * The default class is not implemented.  It's functions should never
   * be called because partial template specialization should always be
   * done on the vector traits.
   *
   * This class is originally from the Programming Languages for Mathematicians
   * class taught by Erich Kaltofen at NCSU in Fall 1998.
   * @param Field LinBox field
   * @param Row   LinBox sparse vector type to use for rows of matrix
   * @param Vector LinBox dense or sparse vector
   * @param Trait  Marker whether to use dense or sparse LinBox vector 
   *               implementation.  This is chosen by a default parameter 
   *               and partial template specialization.
   */
  template <class Field, class Row, class Vector,
	    class Trait = vector_traits<Vector>::vector_category>
  class sparsemat_aux : public sparsemat_base<Field, Row>
  {
  public:

    /** Constructor.
      * Note: the copy constructor and operator= will work as intended
      *       because of STL's container design
      * @param  F  the field of entries; passed so that a possible paramter 
      * 	   such as a modulus is known to the matrix.
      * @param  m  row dimension
      * @param  n  column dimension
      */
    sparsemat_aux(const Field& F, size_t m, size_t n); 

    /** Constructor from sparsemat_base.
      * @param  A sparsemat_base object of same Field and Row types.
      */
    sparsemat_aux(const sparsemat_base<Field, Row>& B); 

    /** Destructor. */
    ~sparsemat_aux(); // {}

    /** Linear Solver.
     * Solves linear system Ax = b by performing Gaussian elimination
     * and back substitution.  If the system does not have a uinque
     * solution, an error is printed and a zero matrix is returned.
     * Templatized by LinBox vector class.
     * @param	   b must be encapsulated vector of size m.
     * @return     encapsualted vector of size n such that A * x = b
     */
    Vector& linsolve(Vector& b);

    /** Gaussian elimination.
     * Performs gaussian elimination of matrix A and corresponding
     * right side vector b.
     * Templatized by LinBox vector class.
     * @param  b corresponding right side vector b.  
     *         (Default is empty vector.)
     * @return boolean true if succesful, false if not
     * @see www.math.ncsu.edu/~kaltofen/courses/LinAlgebra/Maple/refpkg/Src/ref.mpl
     * 	  for a Maple row echelon form algorithm
     */
    bool gauss(Vector& b = Vector() );
 
    /** Application of BlackBox matrix.
     * y = A*x.
     * Templatized by LinBox vector class.
     * @return reference to output vector y
     * @param  x input vector
     */
    Vector& apply(const Vector& x) const;

    /** Application of BlackBox matrix transpose.
     * y = transpose(A)*x.
     * Templatized by LinBox vector class.
     * @return reference to output vector y
     * @param  x input vector
     */
    Vector& applyTranspose(const Vector& x) const;

  };// sparsemat_aux
	   
  // Specialization of sparsemat_aux for LinBox dense vectors
  template <class Field, class Row, class Vector>
  class sparsemat_aux<Field, Row, Vector, vector_categories::dense_vector_tag>
    : public sparsemat_base<Field, Row>
  {
  public:

    sparsemat_aux(const Field& F, size_t m, size_t n) 
      : sparsemat_base<Field, Row>(F, m, n) {}
    sparsemat_aux(const sparsemat_base<Field, Row>& B)
      : sparsemat_base<Field, Row>(B) {}
    ~sparsemat_aux() {}
    Vector& linsolve(Vector& b);
    bool gauss(Vector& b = Vector() );
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const;

  };// sparsemat_aux<dense_vector_tag>
	  
  // Specialization of sparsemat_aux for LinBox sparse sequence vectors
  template <class Field, class Row, class Vector>
  class sparsemat_aux<Field, 
                      Row, 
		      Vector, 
		      vector_categories::sparse_sequence_vector_tag>
    : public sparsemat_base<Field, Row>
  {
  public:

    sparsemat_aux(const Field& F, size_t m, size_t n) 
      : sparsemat_base<Field, Row>(F, m, n) {}
    sparsemat_aux(const sparsemat_base<Field, Row>& B)
      : sparsemat_base<Field, Row>(B) {}
    ~sparsemat_aux() {}
    Vector& linsolve(Vector& b);
    bool gauss(Vector& b = Vector() );
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const;

  private:
    // used in lower_bound as function object
    struct comp_w_index 
    {
      bool operator()
      (const pair< size_t, element >& entry, size_t col_in)

      { return entry.first < col_in; }
    }; // struct comp_w_index

  };// sparsemat_aux<sparse_sequence_vector_tag>
	  
  // Specialization of sparsemat_aux for LinBox sparse associative vectors
  template <class Field, class Row, class Vector>
  class sparsemat_aux<Field, 
                      Row, 
		      Vector, 
		      vector_categories::sparse_associative_vector_tag>
    : public sparsemat_base<Field, Row>
  {
  public:

    sparsemat_aux(const Field& F, size_t m, size_t n) 
      : sparsemat_base<Field, Row>(F, m, n) {}
    sparsemat_aux(const sparsemat_base<Field, Row>& B)
      : sparsemat_base<Field, Row>(B) {}
    ~sparsemat_aux() {}
    Vector& linsolve(Vector& b);
    bool gauss(Vector& b = Vector() );
    Vector& apply(const Vector& x) const;
    Vector& applyTranspose(const Vector& x) const;

  };// sparsemat_aux<sparse_associative_vector_tag>

  // Implementation of matrix methods for dense vectors

  template <class Field, class Row, class Vector>
  Vector& 
  sparsemat_aux<Field, Row, Vector, vector_categories::dense_vector_tag>
  ::linsolve(Vector& b)
  {
  // exceptions: if the system does not have a unique solution,
  //		 an error is printed and a zero matrix returned
    size_t bs = b.size();
    
    // Create vector to hold output
    typename Field::element zero;
    _F.init(zero,0);
    Vector* x_ptr = new Vector(_n,zero);
    
    // Check to see if dimensions are compatible for solving the system.
    if ( bs != _m )
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return *x_ptr;
    }

    // make copies so can restore originals
    std::vector< row_type > _AA(_A);
    Vector bb(b);

    // perform gaussian elimination
    if (!gauss(bb))
    {
      cerr << endl
	   << "ERROR:  Gaussian elimination not performed." << endl << endl;
      return *x_ptr;
    } // if (!gauss(bb))

    // Back substitution to construct solution

    size_t i, j, k;

    // find first pivot element from bottom by scanning rows
    for (i = _m - 1; i >= 0; i--)
      if ( (_A[i]).begin() != (_A[i]).end() )
      {
  	j = _A[i].begin()->first;
  	break;
      }

    // Check for inconsistancy
    for (k = i+1; k < _m; k++)
      if ( !_F.isZero(b[k]) )
      {
  	cerr << endl 
	     << "ERROR:  Found inconsistency in row " << k << "." << endl
	     << "	 No solution is possible." << endl << endl;
  	_A = _AA;
  	return *x_ptr;
      }

    if ( i < 0 ) // Zero matrix!
    {
      cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
		   << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }

    if ( j != _n-1 ) // Check to see if unique solution
    {
      cerr << endl << "ERROR:  No uniquie solution." << endl
	           << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }

    element value;
    _F.init(value, 0);
    row_iter iter;

    while ( i >= 0 )
    {
      iter = _A[i].begin();
      (*x_ptr)[(iter->first)] = bb[(iter->first)];
      iter++;

      for ( ; iter != _A[i].end(); iter++ ) // Subtract left values.
  	_F.subin((*x_ptr)[(_A[i].begin()->first)],
  		 _F.mul(value,
  			iter->second,
  			(*x_ptr)[(_A[iter->first].begin()->first)]));

      _F.divin((*x_ptr)[(_A[i].begin()->first)], _A[i].begin()->second);
  	// Divide by pivot

      if ( i > 0 )  // Find new pivot element.
      {
	i--; // move to previous row

  	k = _A[i].begin()->first;

  	if ( k < j-1 )  // Check to see if unique solution (redundant)
  	{
  	  cerr << endl << "ERROR:  No uniquie solution." << endl
	               << "	   No solution found." << endl
		       << "	   (Redundant Check)" << endl << endl;
  	  _A = _AA;
  	  return *x_ptr;
  	} // checking for unique solution

  	j = k;
      } // if ( i >= 0 )  // Find new pivot element.
      else
	break;

    } // while (i > 0)

    // retore original matrix A and vector b.
    _A = _AA;
    return *x_ptr;

  } // Vector& sparsemat_aux<dense_vector_tag>::linsolve(Vector&)

  template <class Field, class Row, class Vector>
  bool sparsemat_aux<Field, Row, Vector, vector_categories::dense_vector_tag>
  ::gauss(Vector& b)
  {
    size_t i, j, k, bs;

    bs = b.size();

    // Check to see if dimensions are compatible for solving the system.
    // Empty vectors are ignored
    if ( (bs != 0) && (bs != _m) ) 
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return false;
    }

    i = 0;  // current row
    j = 0;  // current column
    element value, temp;
    _F.init(value, 0);
    _F.init(temp, 0);

    while ( (i < _m - 1) && (j < _n) )
    {
#ifdef TRACE
      clog << "In row " << i << " and column " << j << endl;
#endif // TRACE
      
      // look for non zero element in column
      for (k = i; k < _m; k++) 
	if ( _A[k].begin()->first == j ) break;

      if (k < _m)  // found pivot
      {
#ifdef TRACE
	clog << "  found pivot in row " << k << " and column " 
	     << _A[k].begin()->first << endl;
#endif // TRACE
	
  	if (i != k)  // if current row has zero (i.e., not pivot), swap rows
  	{
#ifdef TRACE
	  clog << "  pivot is not in current row.  swapping rows " 
	       << i << " and " << k << endl;
#endif // TRACE
	  
  	  swaprow(i,k);
	  
	  if (bs != 0) // swap vector elements if not empty
	  {
	    temp = b[i];
	    b[i] = b[k];
	    b[k] = temp;
	  }
  	} // if

#ifdef TRACE
	clog << "  zeroing out rest of column starting at row " 
	     << i + 1 << endl;
#endif // TRACE

  	for (k = i+1; k < _m; k++)
	{
#ifdef TRACE
  	  clog << "    zeroing out row " << k 
	       << ", which has its first element in column " 
	       <<  _A[k].begin()->first << endl;
#endif // TRACE

	  if ( _A[k].begin()->first == j )
	  {
	    _F.negin(_F.div(value,
			    _A[k].begin()->second,
			    _A[i].begin()->second));

#ifdef TRACE
    	    clog << "      adding ";
	    _F.write(clog, value);
	    clog << " times row " << i << " to row " << k << endl;
#endif // TRACE

	    addrow(i, k, value);
	    if (bs != 0)
	    {
	      _F.mul(temp, value, b[i]);
	      _F.addin(b[k], temp);
	    }

#ifdef TRACE
	    clog << "      Matrix is now:" << endl;
	    write(clog);
	    clog << "      and vector is now:" << endl;
	    vector_utility<Field, Vector> VU(_F);
	    VU.write(clog, b);
#endif // TRACE

	  } // if ( _A[k].begin()->first == j )
	} // for (k = i+1; k < _m; k++)
	
	i++;
  	j++;
      } // found pivot
      else  // if there is no pivot, go to next column
  	j++;
    } // while

    return true;
  } //  bool sparsemat_aux<dense_vector_tag>::gauss(Vector&)

  template <class Field, class Row, class Vector>
  inline Vector&
  sparsemat_aux<Field, Row, Vector, vector_categories::dense_vector_tag>
  ::apply(const Vector& x) const
  {
    if (_n != x.size())
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    const_row_iter iter;
    typename Vector::iterator y_iter;

    element temp;
    _F.init(temp, 0);
 
    Vector* y_ptr = new Vector(_m, temp);  // Create output vector of zeros

    y_iter = y_ptr->begin();
    for (size_t i = 0; i < _m; i++, y_iter++)
      for (iter = _A[i].begin(); iter != _A[i].end(); iter++)
  	_F.addin(*y_iter, _F.mul(temp, (*iter).second, x[(*iter).first]));

    return *y_ptr;
 
  } // Vector& sparsemat_aux<dense_vector_tag>::apply(const Vector&) const

  template <class Field, class Row, class Vector>
  inline Vector&
  sparsemat_aux<Field, Row, Vector, vector_categories::dense_vector_tag>
  ::applyTranspose(const Vector& x) const
  {
    if (_m != x.size())
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    const_row_iter iter;
    typename Vector::iterator y_iter;

    element temp;
    _F.init(temp, 0);
 
    Vector* y_ptr = new Vector(_n, temp);  // Create output vector of zeros

    for (size_t i = 0; i < _m; i++, y_iter++)
      for (iter = _A[i].begin(); iter != _A[i].end(); iter++)
	_F.addin((*y_ptr)[(*iter).first], _F.mul(temp, (*iter).second, x[i]));
    
    return *y_ptr;
 
  } // Vector& sparsemat_aux<dense_vector_tag>::applyTranspose(const Vector&) const
  
  // Implementation of matrix methods for sparse sequence vectors

  template <class Field, class Row, class Vector>
  Vector& sparsemat_aux<Field, 
                        Row, 
			Vector, 
			vector_categories::sparse_sequence_vector_tag>
  ::linsolve(Vector& b)
  { 
  // exceptions: if the system does not have a unique solution,
  //		 an error is printed and a zero matrix returned
    size_t bs = b.size();
    
    // Create vector to hold output
    Vector* x_ptr = new Vector();
    
    // Check to see if dimensions are compatible for solving the system.
    if ( (bs != 0) && (_m < b.back().first) )
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return *x_ptr;
    } // if ( (bs != 0) && (_m < b.back().first) )

    // make copies so can restore originals
    std::vector< row_type > _AA(_A);
    Vector bb(b);

    // perform gaussian elimination
    if (!gauss(bb))
    {
      cerr << endl
	   << "ERROR:  Gaussian elimination not performed." << endl << endl;
      return *x_ptr;
    } // if (!gauss(bb))

    // Back substitution to construct solution

    size_t i, j, k;

    // find first pivot element from bottom by scanning rows
    for (i = _m - 1; i >= 0; i--)
      if ( (_A[i]).begin() != (_A[i]).end() )
      {
  	j = _A[i].begin()->first;
  	break;
      } // if ( (_A[i]).begin() != (_A[i]).end() )

    // Check for inconsistancy
    for (k = i+1; k < _m; k++)
      if ( (bb.back().first == k) && (!_F.isZero(bb.back().second)) )
      {
  	cerr << endl 
	     << "ERROR:  Found inconsistency in row " << k << "." << endl
	     << "	 No solution is possible." << endl << endl;
  	_A = _AA;
  	return *x_ptr;
      }

    if ( i < 0 ) // Zero matrix!
    {
      cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
		   << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }

    if ( j != _n-1 ) // Check to see if unique solution
    {
      cerr << endl << "ERROR:  No uniquie solution." << endl
	           << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }
   
    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
    row_iter iter;
    typename Vector::iterator b_iter, x_iter;

    while ( i >= 0 )
    {
      // Set iterators.
      iter = _A[i].begin();
      iter++;
      x_iter = x_ptr->begin();

      // Obtain original value to start subtractions
      b_iter = lower_bound(bb.begin(), bb.end(), i, comp_w_index());
      if ( (b_iter != bb.end()) && (i == b_iter->first) )
	value = b_iter->second;
      else
	value = zero;

      for ( ; iter != _A[i].end(); iter++ ) // Subtract left values.
      {
	k = _A[iter->first].begin()->first;
      
	x_iter = lower_bound(x_iter, x_ptr->end(), k, comp_w_index());

	if ( (x_iter != x_ptr->end()) && (k == x_iter->first) )
	  _F.subin(value, _F.mul(temp, iter->second, x_iter->second));
      } // for ( ; iter != _A[i].end(); iter++ ) // Subtract left values.

      _F.divin(value, _A[i].begin()->second); // Divide by pivot

      // Insert non-zero element in solution vector
      if (!_F.isZero(value))
	x_ptr->insert(x_ptr->begin(), make_pair(i, value));

      if ( i > 0 )  // Find new pivot element.
      {
        i--; // move to previous row
  	
	k = _A[i].begin()->first;

  	if ( k < j-1 )  // Check to see if unique solution (redundant)
  	{
  	  cerr << endl << "ERROR:  No uniquie solution." << endl
	               << "	   No solution found." << endl
		       << "	   (Redundant Check)" << endl << endl;
  	  _A = _AA;
  	  return *x_ptr;
  	} // if ( k < j-1 )  // Check to see if unique solution (redundant)

  	j = k;

      } // if ( i >= 0 )  // Find new pivot element.
      else
	break;

    } // while (i >= 0)

    // retore original matrix A and vector b.
    _A = _AA;

    return *x_ptr;

  } // Vector& sparsemat_aux<sparse_sequence_vector_tag>::linsolve(Vector&)

  template <class Field, class Row, class Vector>
  bool sparsemat_aux<Field, 
                     Row, 
		     Vector, 
		     vector_categories::sparse_sequence_vector_tag>
  ::gauss(Vector& b)
  {
    size_t i, j, k, bs;

    bs = b.size();

    // Check to see if dimensions are compatible for solving the system.
    // Empty vectors are ignored
    if ( (bs != 0) && (_m < b.back().first) )
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return false;
    } // if ( (bs != 0) && (m < b.back().first) )

    i = 0;  // current row
    j = 0;  // current column
    element value, temp;
    _F.init(value, 0);
    _F.init(temp, 0);
    typename Vector::iterator iter_i, iter_k;

    while ( (i < _m - 1) && (j < _n) )
    {
      // look for non zero element in column
      for (k = i; k < _m; k++) 
	if ( _A[k].begin()->first == j ) break;

      if (k < _m)  // found pivot
      {
  	if (i != k)  // if current row has zero, swap rows
  	{
  	  swaprow(i,k);
	  
	  if (bs != 0) // swap vector elements if not empty
	  {
	    // find where elements occur in vector
	    // element k will always be after element i.
	    iter_i = lower_bound(b.begin(), b.end(), i, comp_w_index());
	    iter_k = lower_bound(iter_i, b.end(), k, comp_w_index());

	    // perform swap
	    if ( (iter_i != b.end()) && (iter_k != b.end()) 
		 && (i == iter_i->first) && (k == iter_k->first) )
	    {
              temp = iter_i->second;
	      iter_i->second = iter_k->second;
	      iter_k->second = temp;
	    } // if ( (i == iter_i->first) && (k == iter_k->first) )
            else if ( ( (iter_i != b.end()) && (i == iter_i->first) )
		      && ( (iter_k == b.end()) || (k != iter_k->first) ) )
	    {
              b.insert(iter_k, make_pair(k, iter_i->second));
	      iter_i = lower_bound(b.begin(), b.end(), i, comp_w_index());
	      b.erase(iter_i);
	    } // else if ( (i == iter_i->first) && (k != iter_k->first) )
            else if ( ( (iter_k != b.end()) && (k == iter_k->first) )
		      && ( (iter_i == b.end()) || (i != iter_i->first) ) )
	    {
              iter_i = b.insert(iter_i, make_pair(i, iter_k->second));
	      iter_k = lower_bound(iter_i, b.end(), k, comp_w_index());
	      b.erase(iter_k);
	    } // else if ( (k == iter_k->first) && (i != iter_i->first) )
	  } // if (bs != 0) // swap vector elements if not empty
  	} // if (i != k)  // if current row has zero, swap rows

  	for (k = i+1; k < _m; k++)  // zero out rest of column
	  if ( _A[k].begin()->first == j )
	  {
	    _F.negin(_F.div(value,
			    _A[k].begin()->second,
			    _A[i].begin()->second));
	    addrow(i, k, value);
	    if (bs != 0)
	    {
	      // find where elements occur in vector
	      // element k will always be after element i.
	      iter_i = lower_bound(b.begin(), b.end(), i, comp_w_index());
	      iter_k = lower_bound(iter_i, b.end(), k, comp_w_index());
	      
  	      // perform swap
	      if ( (iter_i != b.end()) && (i == iter_i->first) )
	      {
		if ( (iter_k != b.end()) && (k == iter_k->first) )
		{
		  _F.addin(iter_k->second, _F.mul(temp, value, iter_i->second));
		} // if (k == iter_k->first)
		else
		  b.insert(iter_k, 
			   make_pair(k, 
				     _F.mul(temp, value, iter_i->second)));
	      } // if (i == iter_i->first)
	    } // if (bs != 0)
	  } // if ( _A[k].begin()->first == j ) and for (k = i+1; k < _m; k++)

	i++;
  	j++;
      } // if (k < _m)  // found pivot
      else  // if there is no pivot, go to next column
  	j++;
    } // while ( (i < _m - 1) && (j < _n) )

    return true;
  } //  bool sparsemat_aux<sparse_sequence_vector_tag>::gauss(Vector&)

  template <class Field, class Row, class Vector>
  inline Vector& sparsemat_aux<Field, 
                               Row, 
			       Vector, 
			       vector_categories::sparse_sequence_vector_tag>
  ::apply(const Vector& x) const
  {
    if ( (x.size() != 0) && (_n < x.back().first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    Vector* y_ptr = new Vector();  // Create output vector of zeros
    
    const_row_iter iter;
    typename Vector::const_iterator x_iter;

    size_t k;

    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
 
    for (size_t i = 0; i < _m; i++)
    {
      value = zero;
      x_iter = x.begin();

      for (iter = _A[i].begin(); iter != _A[i].end(); iter++)
      {
	k = (*iter).first;
      
	x_iter = lower_bound(x_iter, x.end(), k, comp_w_index());

	if ( (x_iter != x.end()) && (k == (*x_iter).first) )
	  _F.addin(value, _F.mul(temp, (*iter).second, (*x_iter).second));
      } // for (iter = _A[i].begin(); iter != _A[i].end(); iter++)

      // Insert non-zero element in solution vector
      if (!_F.isZero(value))
	y_ptr->push_back(make_pair(i, value));

    } // for (size_t i = 0; i < _m; i++)
    
    return *y_ptr;
  } // Vector& sparsemat_aux<sparse_sequence_vector_tag>::apply(const Vector&) const

  template <class Field, class Row, class Vector>
  inline Vector&
  sparsemat_aux<Field, 
                Row, 
		Vector, 
		vector_categories::sparse_sequence_vector_tag>
  ::applyTranspose(const Vector& x) const
  {
    if ( (x.size() != 0) && (_m < x.back().first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    Vector* y_ptr = new Vector();  // Create output vector of zeros
    
    const_row_iter iter;
    typename Vector::const_iterator x_iter;

    size_t k;

    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
 
    std::vector<element> y(_n, zero); // temporary vector for calculating output
    
    for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
    {
      k = (*x_iter).first;       // vector index
      value = (*x_iter).second;  // value in vector

      // apply vector to column k
      for (iter = _A[k].begin(); iter != _A[k].end(); iter++)
	_F.addin(y[(*iter).first], _F.mul(temp, (*iter).second, value));
             
    } // for (x_iter = x.begin(); x_iter != x.end(); x_iter++)

    // Convert temporary vector to sparse vector for output
    for (size_t i = 0; i < y.size(); i++)
      if (!_F.isZero(y[i])) y_ptr->push_back(make_pair(i, y[i]));
 
    return *y_ptr;
  } // Vector& sparsemat_aux<sparse_sequence_vector_tag>::applyTranspose(...) const

  // Implementation of matrix methods for sparse associative vectors

  template <class Field, class Row, class Vector>
  Vector& 
  sparsemat_aux<Field, 
                Row, 
		Vector, 
		vector_categories::sparse_associative_vector_tag>
  ::linsolve(Vector& b)
  { 
    // exceptions: if the system does not have a unique solution,
    //		 an error is printed and a zero matrix returned
    size_t bs = b.size();
    
    // Create vector to hold output
    Vector* x_ptr = new Vector();
    
    // Check to see if dimensions are compatible for solving the system.
    if ( (bs != 0) && (_m < b.rbegin()->first) )
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return *x_ptr;
    } // if ( (bs != 0) && (_m < b.rbegin()->first) )

    // make copies so can restore originals
    std::vector< row_type > _AA(_A);
    Vector bb(b);

    // perform gaussian elimination
    if (!gauss(bb))
    {
      cerr << endl
	   << "ERROR:  Gaussian elimination not performed." << endl << endl;
      return *x_ptr;
    } // if (!gauss(bb))

    // Back substitution to construct solution

    size_t i, j, k;

    // find first pivot element from bottom by scanning rows
    for (i = _m - 1; i >= 0; i--)
      if ( (_A[i]).begin() != (_A[i]).end() )
      {
  	j = _A[i].begin()->first;
  	break;
      } // if ( (_A[i]).begin() != (_A[i]).end() )

    // Check for inconsistancy
    for (k = i+1; k < _m; k++)
      if ( (bb.rbegin()->first == k) && (!_F.isZero(bb.rbegin()->second)) )
      {
  	cerr << endl 
	     << "ERROR:  Found inconsistency in row " << k << "." << endl
	     << "	 No solution is possible." << endl << endl;
  	_A = _AA;
  	return *x_ptr;
      }

    if ( i < 0 ) // Zero matrix!
    {
      cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
		   << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }

    if ( j != _n-1 ) // Check to see if unique solution
    {
      cerr << endl << "ERROR:  No uniquie solution." << endl
	           << "        No solution found." << endl << endl;
      _A = _AA;
      return *x_ptr;
    }
   
    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
    row_iter iter;
    typename Vector::iterator b_iter, x_iter;

    while ( i >= 0 )
    {
      // Set iterators.
      iter = _A[i].begin();
      iter++;
      x_iter = x_ptr->begin();

      // Obtain original value to start subtractions
      b_iter = bb.find(i);
      if (b_iter != bb.end())
	value = b_iter->second;
      else
	value = zero;

      for ( ; iter != _A[i].end(); iter++ ) // Subtract left values.
      {
	k = _A[iter->first].begin()->first;
      
	x_iter = x_ptr->find(k);

	if (x_iter != x_ptr->end())
	  _F.subin(value, _F.mul(temp, iter->second, x_iter->second));
      } // for ( ; iter != _A[i].end(); iter++ ) // Subtract left values.

      _F.divin(value, _A[i].begin()->second); // Divide by pivot

      // Insert non-zero element in solution vector
      if (!_F.isZero(value))
	x_ptr->insert(x_ptr->begin(), make_pair(i, value));

      if ( i > 0 )  // Find new pivot element.
      {
        i--; // move to previous row
  	
	k = _A[i].begin()->first;

  	if ( k < j-1 )  // Check to see if unique solution (redundant)
  	{
  	  cerr << endl << "ERROR:  No uniquie solution." << endl
	               << "	   No solution found." << endl
		       << "	   (Redundant Check)" << endl << endl;
  	  _A = _AA;
  	  return *x_ptr;
  	} // if ( k < j-1 )  // Check to see if unique solution (redundant)

  	j = k;

      } // if ( i >= 0 )  // Find new pivot element.
      else
	break;

    } // while (i >= 0)

    // retore original matrix A and vector b.
    _A = _AA;

    return *x_ptr;
  } // Vector& sparsemat_aux<sparse_associative_vector_tag>::linsolve(Vector&)

  template <class Field, class Row, class Vector>
  bool sparsemat_aux<Field, 
                     Row, 
		     Vector, 
		     vector_categories::sparse_associative_vector_tag>
  ::gauss(Vector& b)
  {
    size_t i, j, k, bs;

    bs = b.size();

    // Check to see if dimensions are compatible for solving the system.
    // Empty vectors are ignored
    if ( (bs != 0) && (_m < b.rbegin()->first) )
    {
      cerr << endl
	   << "ERROR:  Matrix and vector indices are incompatible." << endl
	   << "        No solution is possible." << endl << endl;
      return false;
    } // if ( (bs != 0) && (_m < b.rbegin()->first) )

    i = 0;  // current row
    j = 0;  // current column
    element value, temp;
    _F.init(value, 0);
    _F.init(temp, 0);
    typename Vector::iterator iter_i, iter_k;

    while ( (i < _m - 1) && (j < _n) )
    {
      // look for non zero element in column
      for (k = i; k < _m; k++) 
	if ( _A[k].begin()->first == j ) break;

      if (k < _m)  // found pivot
      {
  	if (i != k)  // if current row has zero, swap rows
  	{
  	  swaprow(i,k);
	  
	  if (bs != 0) // swap vector elements if not empty
	  {
	    // find where elements occur in vector
	    // element k will always be after element i.
	    iter_i = b.find(i);
	    iter_k = b.find(k);

	    // perform swap
	    if ( (iter_i != b.end()) && (iter_k != b.end()) )
	    {
              temp = iter_i->second;
	      iter_i->second = iter_k->second;
	      iter_k->second = temp;
	    } // if ( (iter_i != b.end()) && (iter_k != b.end()) )
            else if ( (iter_i != b.end()) && (iter_k == b.end()) )
	    {
              b.insert(make_pair(k, iter_i->second));
	      iter_i = b.find(i);
	      b.erase(iter_i);
	    } // else if ( (iter_i != b.end()) && (iter_k == b.end()) )
            else if ( (iter_k != b.end()) && (iter_i == b.end()) )
	    {
              b.insert(make_pair(i, iter_k->second));
	      iter_k = b.find(k);
	      b.erase(iter_k);
	    } // else if ( (iter_k != b.end()) && (iter_i == b.end()) )
	  } // if (bs != 0) // swap vector elements if not empty
  	} // if (i != k)  // if current row has zero, swap rows

  	for (k = i+1; k < _m; k++)  // zero out rest of column
	  if ( _A[k].begin()->first == j )
	  {
	    _F.negin(_F.div(value,
			    _A[k].begin()->second,
			    _A[i].begin()->second));
	    addrow(i, k, value);
	    if (bs != 0)
	    {
	      // find where elements occur in vector
	      // element k will always be after element i.
	      iter_i = b.find(i);
	      iter_k = b.find(k);
	      
  	      // perform addition
	      if (iter_i != b.end())
	      {
		if (iter_k != b.end())
		{
		  _F.addin(iter_k->second, _F.mul(temp, value, iter_i->second));
		} // if (iter_k != b.end())
		else
		  b.insert(make_pair(k, 
				     _F.mul(temp, value, iter_i->second)));
	      } // if (iter_i != b.end())
	    } // if (bs != 0)
	  } // if ( _A[k].begin()->first == j ) and for (k = i+1; k < _m; k++)

	i++;
  	j++;
      } // if (k < _m)  // found pivot
      else  // if there is no pivot, go to next column
  	j++;
    } // while ( (i < _m - 1) && (j < _n) )

    return true;
  } //  bool sparsemat_aux<sparse_associative_vector_tag>::gauss(Vector&)

  template <class Field, class Row, class Vector>
  inline Vector&
  sparsemat_aux<Field, 
                Row, 
		Vector, 
		vector_categories::sparse_associative_vector_tag>
  ::apply(const Vector& x) const
  {
    if ( (x.size() != 0) && (_n < x.rbegin()->first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    Vector* y_ptr = new Vector();  // Create output vector of zeros
    
    const_row_iter iter;
    typename Vector::const_iterator x_iter;

    size_t k;

    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
 
    for (size_t i = 0; i < _m; i++)
    {
      value = zero;
      x_iter = x.begin();

      for (iter = _A[i].begin(); iter != _A[i].end(); iter++)
      {
	k = (*iter).first;
      
	x_iter = x.find(k);

	if (x_iter != x.end())
	  _F.addin(value, _F.mul(temp, (*iter).second, (*x_iter).second));
      } // for (iter = _A[i].begin(); iter != _A[i].end(); iter++)

      // Insert non-zero element in solution vector
      if (!_F.isZero(value))
	y_ptr->insert(make_pair(i, value));

    } // for (size_t i = 0; i < _m; i++)
    
    return *y_ptr;
  } // Vector& sparsemat_aux<sparse_associative_vector_tag>::apply(...) const

  template <class Field, class Row, class Vector>
  inline Vector&
  sparsemat_aux<Field, 
                Row, 
		Vector, 
		vector_categories::sparse_associative_vector_tag>
  ::applyTranspose(const Vector& x) const
  {
    if ( (x.size() != 0) && (_m < x.rbegin()->first) )
    {
      cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
      return *(new Vector);
    }
 
    Vector* y_ptr = new Vector();  // Create output vector of zeros
    
    const_row_iter iter;
    typename Vector::const_iterator x_iter;

    size_t k;

    element zero;
    _F.init(zero, 0);
    element value(zero), temp(zero);
 
    std::vector<element> y(_n, zero); // temporary vector for calculating output
    
    for (x_iter = x.begin(); x_iter != x.end(); x_iter++)
    {
      k = (*x_iter).first;       // vector index
      value = (*x_iter).second;  // value in vector

      // apply vector to column k
      for (iter = _A[k].begin(); iter != _A[k].end(); iter++)
	_F.addin(y[(*iter).first], _F.mul(temp, (*iter).second, value));
             
    } // for (x_iter = x.begin(); x_iter != x.end(); x_iter++)

    // Convert temporary vector to sparse vector for output
    for (size_t i = 0; i < y.size(); i++)
      if (!_F.isZero(y[i])) y_ptr->insert(make_pair(i, y[i]));
 
    return *y_ptr;
  } // Vector& sparsemat_aux<sparse_associative_vector_tag>::applyTranspose(...) const

} // namespace LinBox

#endif // _SPARSEMAT_AUX_
