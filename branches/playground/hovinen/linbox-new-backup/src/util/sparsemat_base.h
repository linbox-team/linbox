/* File: sparsemat_base.h
 * Author: William J. Turner for the LinBox group
 */

#ifndef _SPARSEMAT_BASE_
#define _SPARSEMAT_BASE_

#include <vector>    // STL vector
#include <utility>   // STL pair
#include <iostream>
#include <algorithm>

#include "LinBox/faxpy.h"
#include "LinBox/vector_traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
  /** Auxillary sparse matrix base template.
   * This is a class of sparse matrices templatized by the field in
   * which the elements reside.  The matrix itself is stored as an
   * STL vector of \Ref{LinBox} sparse vectors of integers and field elements.
   * Each sparse vector corresponds to one row of the matrix, 
   * and each pair (j, a) in sparse vector i corresponds to the (i,j) 
   * entry of the matrix.
   *
   * It is templatized by the \Ref{LinBox} field in which the arithmetic
   * is done, the sparse vector type with which the rows of the matrix 
   * are implemented, and the vector category of the row implementation.
   * This third template parameter is defaulted to be the \Ref{LinBox} 
   * vector trait of the vector.  This class is then specialized for 
   * sequence and associative sparse vectors.
   *
   * This class does not contain any functions operating on vectors.
   * These functions must all be implemented in derived classes which must
   * be able to access the actual data structes stored in this class.  
   * This separation occurs to avoid implementing all methods twice.  
   * Now, only the methods that depend on the vectors are implemented twice.
   *
   * This class is originally from the Programming Languages for Mathematicians
   * class taught by Erich Kaltofen at NCSU in Fall 1998.
   * @param Field LinBox field
   * @param Row   LinBox sparse vector type to use for rows of matrix
   */
  template <class Field, 
            class Row, 
	    class Trait = vector_traits<Row>::vector_category>
  class sparsemat_base
  {
  public:

    /// Element type
    typedef typename Field::element element;
 
    /// Row type
    typedef Row row_type;

    /// Row iterator type
    typedef typename row_type::iterator row_iter;

    /// Constant row iterator
    typedef typename row_type::const_iterator const_row_iter;

    /** Constructor.
      * Note: the copy constructor and operator= will work as intended
      *       because of STL's container design
      * @param  F  the field of entries; passed so that a possible paramter 
      * 	   such as a modulus is known to the matrix.
      * @param  m  row dimension
      * @param  n  column dimension
      */
    sparsemat_base(const Field& F, size_t m, size_t n);

    /** Destructor. */
    ~sparsemat_base() {}

    /** Retreive row dimensions of Sparsemat matrix.
      * @return integer number of rows of sparsemat_base matrix.
      */
    size_t get_rowdim(void) const { return _m; }

    /** Retreive column dimensions of Sparsemat matrix.
      * @return integer number of columns of sparsemat_base matrix.
      */
    size_t get_coldim(void) const { return _n; }

    /** Retrieve matrix element.
      * If the indices are out of range, an error is printed and
      * the zero element is returned.
      * Unlike STL vector's operator[] this member does not
      * return a reference to the Field element in the matrix
      * it is therefore impossible to write
      * A[make_pair(3, 5)] = element(1);
      * the reason is that zero entries have no memory where
      * one could place the right side Field element
      * @return     a copy of the element in row i, column j
      * @param      ind pair of indices (i,j)
      */
    element operator[] (const pair<size_t, size_t>& ind) const;

    /** Insert matrix element.
      * sets A[ind.first, ind.second] = a.
      * For example, A.put_value(make_pair(3, 5), element(-4))
      * If the indices are out of range, an error is printed.
      * If the element attempting to be inserted is the zero element,
      * no cell will be inserted and any already existing cell for the
      * entry will be erased.
      * @param ind  pair of integers (i,j) for row i and column j
      * @param a    field element to insert in matrix
      */
    void put_value (const pair<size_t, size_t>& ind, const element& a);

    /** Print matrix.
      * Prints rows as lists.
      * Can be used in operator <<.
      * @param  os  output stream on which to print matrix.
      */
    ostream& write(ostream& os) const;

    /** Read matrix.
      * can be called in operator>>.
      * Works by reading integers for row index and column index
      * and then element to insert in matrix.
      * Stops reading upon non-positive row index or end of file.
      * @param  is  input stream from which to read matrix.
      */
    istream& read(istream& is);

    /** Exchange two rows.
      * Exchanges rows i and j of the matrix.
      * @param  i first row index
      * @param  j second row index
      */
    void swaprow(size_t i, size_t j);

    /** Adds multiple of one row to another.
      * Adds a*(row i) to (row j).
      * @param  i first row index
      * @param  j second row index
      * @param  a multiple of row i to add to row j
      */
    void addrow(size_t i, size_t j, const element& a);

  protected:

    /* Sparse matrix data structure.
     * _A[i], the i-th row, is a LinBox sparse vector
     * Each sparse vector entry is a column index/value pair.
     */
    std::vector< row_type > _A;

    // Field used for all arithmetic
    Field _F;

    // Row dimension from 1..m; the actual dimensions
    size_t _m;

    // Column dimension from 1..n; the actual dimensions
    size_t _n;

  };// end template class sparsemat_base

  // Specialization of sparsemat_base for sequence rows
  template <class Field, class Row>
  class sparsemat_base<Field, 
                       Row, 
		       vector_categories::sparse_sequence_vector_tag>
  {
  public:

    typedef typename Field::element element;
    typedef Row row_type;
    typedef typename row_type::iterator row_iter;
    typedef typename row_type::const_iterator const_row_iter;
    sparsemat_base(const Field& F, size_t m, size_t n);
    ~sparsemat_base() {}
    size_t get_rowdim(void) const { return _m; }
    size_t get_coldim(void) const { return _n; }
    element operator[] (const pair<size_t, size_t>& ind) const;
    void put_value (const pair<size_t, size_t>& ind, const element& a);
    ostream& write(ostream& os) const;
    istream& read(istream& is);
    void swaprow(size_t i, size_t j);
    void addrow(size_t i, size_t j, const element& a);

  protected:

    std::vector< row_type > _A;
    Field _F;
    size_t _m;
    size_t _n;

  private:
    // used in lower_bound as function object
    struct comp_w_col_index 
    {
      bool operator()
      (const pair< size_t, element >& entry, size_t col_in)

      { return entry.first < col_in; }
    }; // struct comp_w_col_index

/*
    // Friend declarations
    friend ostream& operator<< <> (ostream&, const sparsemat_base<Field, Row, vector_categories::sparse_sequence_vector_tag>&);
    friend istream& operator>> <> (istream&, sparsemat_base<Field, Row, vector_categories::sparse_sequence_vector_tag>&);
*/
  };// class sparsemat_base<sparse_sequence_vector_tag>

  // Specialization of sparsemat_base for associative rows
  template <class Field, class Row>
  class sparsemat_base<Field, 
                       Row, 
		       vector_categories::sparse_associative_vector_tag>
  {
  public:

    typedef typename Field::element element;
    typedef Row row_type;
    typedef typename row_type::iterator row_iter;
    typedef typename row_type::const_iterator const_row_iter;
    sparsemat_base(const Field& F, size_t m, size_t n);
    ~sparsemat_base() {}
    size_t get_rowdim(void) const { return _m; }
    size_t get_coldim(void) const { return _n; }
    element operator[] (const pair<size_t, size_t>& ind) const;
    void put_value (const pair<size_t, size_t>& ind, const element& a);
    ostream& write(ostream& os) const;
    istream& read(istream& is);
    void swaprow(size_t i, size_t j);
    void addrow(size_t i, size_t j, const element& a);

  protected:

    std::vector< row_type > _A;
    Field _F;
    size_t _m;
    size_t _n;

  private:
/*
    // Friend declarations
    friend ostream& operator<< <> (ostream&, const sparsemat_base<Field, Row, vector_categories::sparse_associative_vector_tag>&);
    friend istream& operator>> <> (istream&, sparsemat_base<Field, Row, vector_categories::sparse_associative_vector_tag>&);
*/
  };// class sparsemat_base<sparse_associative_vector_tag>

  // Implementation of matrix methods for sparse sequence rows

  template <class Field, class Row>
  inline sparsemat_base<Field, 
                        Row, 
			vector_categories::sparse_sequence_vector_tag>
  ::sparsemat_base(const Field& F, size_t m, size_t n) 
  : // constructs a matrix of the given dimensions with all 0s
    _A(m, row_type() ),
    _F(F),
    _m(m), _n(n) // set the dimensions
  { }

  template <class Field, class Row>
  void sparsemat_base<Field, 
                      Row, 
		      vector_categories::sparse_sequence_vector_tag>
  ::put_value (const pair<size_t, size_t>& ind, const element& a) 
  {
    size_t i = ind.first;
    size_t j = ind.second;
    row_iter iter;
    bool found(true);

    // Check row and column indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) 
    {
      cerr << endl
	   << "ERROR:  Element indices exceed matrix dimensions." << endl
	   << "	 Element not inserted in matrix." << endl << endl;
      return;
    } // if-statement

    // Find appropriate location of element in sparse vector.
    if( (_A[i]).begin() == (_A[i]).end() )
      iter = (_A[i]).end();
    else
      iter = lower_bound( (_A[i]).begin(), 
			  (_A[i]).end(), 
			  j,
			  comp_w_col_index() );

    // Check to see if element already exists.
    if ( (_A[i]).end() == iter )
      found = false;
    else
      if ( iter->first != j )
  	found = false;


    // If element is already in row, replace old value with new.
    // Otherwise, insert the element in the row.
    if (found) 
    {
      if (_F.isZero(a))
  	_A[i].erase(iter);
      else
  	iter->second = a;
    } // if (found)
    else
      if (!_F.isZero(a))
  	_A[i].insert(iter, make_pair(j,a));

  } // void sparsemat_base<sparse_sequence_vector_tag>::put_value(...)

  template <class Field, class Row>
  typename Field::element  
  sparsemat_base<Field, 
                 Row, 
		 vector_categories::sparse_sequence_vector_tag>
  ::operator[] (const pair<size_t, size_t>& ind) const
  {
    size_t i = ind.first;
    size_t j = ind.second;
    element zero;
    _F.init(zero, 0);

    // Check row and column indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) 
    {
      cerr << endl
	   << "ERROR:  Element indices exceed matrix dimensions." << endl
	   << endl;
      return zero;
    }

    row_type row(_A[i]);
    row_iter iter;
    bool found(true);

    // Find appropriate location of element in row.
    if ( row.begin() == row.end() )
      iter = row.end();
    else
      iter = lower_bound( row.begin(), row.end(), j, comp_w_col_index() );
    // Check to see if element exists.
    if ( row.end() == iter )
      found = false;
    else
      if ( iter->first != j )
	found = false;

    // If element is found, return non-zero value.
    // Otherwise value is zero.
    if (found)
      return iter->second;
    else
      return zero;

  } // element sparsemat_base<sparse_sequence_vector_tag>::operator[] (...) const

  template <class Field, class Row>
  inline ostream& sparsemat_base<Field, 
                                 Row, 
				 vector_categories::sparse_sequence_vector_tag>
  ::write(ostream& os) const
  {
    for (size_t i = 0; i <= _m - 1; i++)
    {
      os << "Row " << i << ": ";

      for (const_row_iter iter_i = _A[i].begin(); 
	   iter_i != _A[i].end(); 
	   iter_i++)
      {
//  	os << "(" << iter_i->first << ", ";
//  	_F.write(os, iter_i->second);
  	os << "(" << (*iter_i).first << ", ";
  	_F.write(os, (*iter_i).second);
  	os << ")";
      }
      os << endl;
    }
 
    return os;

  } // ostream& sparsemat_base<sparse_sequence_vector_tag>::write(...) const

  template <class Field, class Row>
  inline istream& sparsemat_base<Field, 
                                 Row, 
				 vector_categories::sparse_sequence_vector_tag>
  ::read(istream& is)
  {
    size_t i, j;
    element el;
    _F.init(el, 0);

    while (is >> i)
      // operator>> returns a reference to an istream
      // then istream::operator void*() is called which
      // returns !basic_ios::fail() [Stroutstrup, p.617, ???]
    {
      if(i == size_t(-1)) break; // return also if row index is -1
      is >> j;
      _F.read(is, el);
      put_value(make_pair(i,j), el);
    } // while-loop

    return is;

  } // istream& sparsemat_base<sparse_sequence_vector_tag>::read(...)

  template <class Field, class Row>
  void sparsemat_base<Field, Row, vector_categories::sparse_sequence_vector_tag>
  ::swaprow(size_t i, size_t j) 
  {
    // exchanges row i and j in A
    // note:	  uses row::swap

    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No rows exchanged." << endl << endl;
      return;
    }

    // Swap rows i and j using row::swap
    _A[i].swap( _A[j] );

    return;
  } // void sparsemat_base<sparse_sequence_vector_tag>::swaprow(...)

  /* This implementation is works for lists and deques, but it causes 
   * segmentation faults for vectors.  Inserting elements into a vector
   * invalidates iterators *before* the inserted element, which is contrary 
   * to the standard.
   * 
   * This may not be a good implementation, anyway, because according
   * to the C++ standard, insertion and erasure can invalidate all iterators 
   * and element references to the sequence.
   *
  template <class Field, class Row>
  void sparsemat_base<Field, Row, vector_categories::sparse_sequence_vector_tag>
  ::addrow(size_t i, size_t j,const element& a) 
  {
    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No row addition preformed." << endl << endl;
      return;
    }

    // Check to see if a is the zero Field element.
    // If so, no addition is performed.
    if (_F.isZero(a)) return;

    // Check to see if row i is empty.  If so, no addition is preformed.
    if ( (_A[i]).begin() == (_A[i]).end() ) return;

    size_t k;
    element value;
    _F.init(value, 0);
    LinBox::faxpy<Field> Faxpy(_F, a);
    // iterators to point to place in rows i and j respectively,
    // and extra iterator for erasing from row j;
    row_iter iter_i, iter_j, iter;

    bool found(true);

    iter_j = (_A[j]).begin(); // start at beginning of second row

    // iterate over elements in row i
    for ( iter_i = (_A[i]).begin(); iter_i != (_A[i]).end(); iter_i++) 
    {
      found = true;
      k = iter_i->first;  // marks current column.

      // Find where column k occurs in row j.
      while ( ( (_A[j]).end() != iter_j ) && ( iter_j->first < k ) )
      {
  	iter_j++;
      }

      // Check if row j has element for column k.
      if ( ( (_A[j]).end() == iter_j ) || ( iter_j->first != k ) )
  	found = false;

      // If row j contains element for column k, perform sum.
      // Otherwise, sum = a * _A[i,k]
      if (found) 
      {
  	if (_F.isZero(Faxpy.applyin(iter_j->second, iter_i->second)))
  	{
  	  iter = iter_j++;
  	  _A[j].erase(iter);
  	}
  	else
  	  iter_j++;

      } // if (found)
      else //if (!found)
  	_A[j].insert(iter_j,
  		     make_pair(k, _F.mul(value,a,iter_i->second)));
	  
    } // end for loop (iteration over elements in row i)

    return;
  } // void sparsemat_base<sparse_sequence_vector_tag>::addrow(...)
*/
  
  /* This implementation works for vectors because it avoids the insert 
   * method.  It creates a new row, using push_back insert new elements, and 
   * then copies it into the _A[j] at the end.  This is less efficient than  
   * doing an inplace row add like above, but no iterators are invalidated 
   * through the insert methods.
   */
  template <class Field, class Row>
  void sparsemat_base<Field, Row, vector_categories::sparse_sequence_vector_tag>
  ::addrow(size_t i, size_t j,const element& a) 
  {
    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No row addition preformed." << endl << endl;
      return;
    }

    // Check to see if a is the zero Field element.
    // If so, no addition is performed.
    if (_F.isZero(a)) return;

    // Check to see if row i is empty.  If so, no addition is preformed.
    if ( (_A[i]).begin() == (_A[i]).end() ) return;

    // variables used in computation
    element value;
    _F.init(value, 0);
    LinBox::faxpy<Field> Faxpy(_F, a);
    row_type row;
    row_iter iter_i, iter_j(_A[j].begin());

    for (iter_i = _A[i].begin(); iter_i != _A[i].end(); iter_i++)
    {
      while( (iter_j != _A[j].end()) && (iter_j->first < iter_i->first) )
      {
	row.push_back(*iter_j);
	iter_j++;
      } // while( (iter_j != _A[j].end()) && (iter_j->first < iter_i->first) )

      if( (iter_j != _A[j].end()) && (iter_j->first == iter_i->first) )
      {
	if (!_F.isZero(Faxpy.apply(value, iter_i->second, iter_j->second)))
	  row.push_back(make_pair(iter_i->first, value));
	
	iter_j++;
      } // if( (iter_j != _A[j].end()) && (iter_j->first == iter_i->first) )
      else
	row.push_back(make_pair(iter_i->first,
				_F.mul(value, a, iter_i->second)));
      
    } // for (iter_i = _A[i].begin(); iter_i != _A[i].end(); iter_i++)

    while(iter_j != _A[j].end())
    {
      row.push_back(*iter_j);
      iter_j++;
    } // while(iter_j != _A[j].end())
    
    _A[j] = row;

    return;
    
  } //  void sparsemat_base<sparse_sequence_vector_tag>::addrow(...)

  // Implementation of matrix methods for sparse associative rows

  template <class Field, class Row>
  inline sparsemat_base<Field, 
                        Row, 
			vector_categories::sparse_associative_vector_tag>
  ::sparsemat_base(const Field& F, size_t m, size_t n) 
  : // constructs a matrix of the given dimensions with all 0s
    _A(m, row_type() ),
    _F(F),
    _m(m), _n(n) // set the dimensions
  { }

  template <class Field, class Row>
  void sparsemat_base<Field, 
                      Row, 
		      vector_categories::sparse_associative_vector_tag>
  ::put_value (const pair<size_t, size_t>& ind, const element& a) 
  {
    size_t i = ind.first;
    size_t j = ind.second;
    row_iter iter;

    // Check row and column indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) 
    {
      cerr << endl
	   << "ERROR:  Element indices exceed matrix dimensions." << endl
	   << "	 Element not inserted in matrix." << endl << endl;
      return;
    } // if-statement

    // Find element in map.  
    // If exists, replace value if not zero, or remove if value is zero.
    // If not found, insert non-zero element
    if ( (iter = _A[i].find(j)) != _A[i].end() )
    {
      if (_F.isZero(a))
	_A[i].erase(iter);
      else
	iter->second = a;
    } // ( (iter = _A[i].find(j)) != _A[i].end() )
    else
    {
      if (!_F.isZero(a))
	_A[i].insert(make_pair(j, a));
    }

  } // void sparsemat_base<sparse_associative_vector_tag>::put_value(...)

  template <class Field, class Row>
  typename Field::element  
  sparsemat_base<Field, 
                 Row, 
		 vector_categories::sparse_associative_vector_tag>
  ::operator[] (const pair<size_t, size_t>& ind) const
  {
    size_t i = ind.first;
    size_t j = ind.second;
    element zero;
    _F.init(zero, 0);

    // Check row and column indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) 
    {
      cerr << endl
	   << "ERROR:  Element indices exceed matrix dimensions." << endl
	   << endl;
      return zero;
    }

    const_row_iter iter;

    if ( (iter = _A[i].find(j)) != _A[i].end() )
      return iter->second;
    else
      return zero;

  } // element sparsemat_base<sparse_associative_vector_tag>::operator[] (...)

  template <class Field, class Row>
  inline ostream& 
  sparsemat_base<Field, 
                 Row, 
		 vector_categories::sparse_associative_vector_tag>
  ::write(ostream& os) const
  {
    for (size_t i = 0; i <= _m - 1; i++)
    {
      os << "Row " << i << ": ";

      for (const_row_iter iter_i = _A[i].begin(); 
	   iter_i != _A[i].end(); 
	   iter_i++)
      {
  	os << "(" << iter_i->first << ", ";
  	_F.write(os, iter_i->second);
  	os << ")";
      }
      os << endl;
    }
 
    return os;

  } // ostream& sparsemat_base<sparse_associative_vector_tag>::write(...) const

  template <class Field, class Row>
  inline istream& 
  sparsemat_base<Field, 
                 Row, 
		 vector_categories::sparse_associative_vector_tag>
  ::read(istream& is)
  {
    size_t i, j;
    element el;
    _F.init(el, 0);

    while (is >> i)
      // operator>> returns a reference to an istream
      // then istream::operator void*() is called which
      // returns !basic_ios::fail() [Stroutstrup, p.617, ???]
    {
      if(i == size_t(-1)) break; // return also if row index is -1
      is >> j;
      _F.read(is, el);
      put_value(make_pair(i,j), el);
    } // while-loop

    return is;

  } // istream& sparsemat_base<sparse_associative_vector_tag>::read(...)

  template <class Field, class Row>
  void sparsemat_base<Field, 
                      Row, 
		      vector_categories::sparse_associative_vector_tag>
  ::swaprow(size_t i, size_t j) 
  {
    // exchanges row i and j in A
    // note:	  uses row::swap

    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No rows exchanged." << endl << endl;
      return;
    }

    // Swap rows i and j using row::swap
    _A[i].swap( _A[j] );

    return;
  } // void sparsemat_base<sparse_associative_vector_tag>::swaprow(...)

  /* This implementation is an inplace row addition, but it isn't clear 
   * from the standard if insertion and deletion invalidates any iterators
   * and element references other than the obvious ones refering to a
   * deleted entries.
   *
  template <class Field, class Row>
  void sparsemat_base<Field, 
                      Row, 
		      vector_categories::sparse_associative_vector_tag>
  ::addrow(size_t i, size_t j,const element& a) 
  {
    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No row addition preformed." << endl << endl;
      return;
    }

    // Check to see if a is the zero Field element.
    // If so, no addition is performed.
    if (_F.isZero(a)) return;

    // Check to see if row i is empty.  If so, no addition is preformed.
    if ( (_A[i]).begin() == (_A[i]).end() ) return;

    size_t k;
    element value;
    _F.init(value, 0);
    LinBox::faxpy<Field> Faxpy(_F, a);
    // iterators to point to place in rows i and j respectively,
    // and extra iterator for erasing from row j;
    row_iter iter_i, iter_j;

    bool found(true);

    iter_j = (_A[j]).begin(); // start at beginning of second row

    // iterate over elements in row i
    for ( iter_i = (_A[i]).begin(); iter_i != (_A[i]).end(); iter_i++) 
    {
      found = true;
      k = iter_i->first;  // marks current column.

      // Find where column k occurs in row j.
      while ( ( (_A[j]).end() != iter_j ) && ( iter_j->first < k ) )
  	iter_j++;

      // Check if row j has element for column k.
      if ( ( (_A[j]).end() == iter_j ) || ( iter_j->first != k ) )
  	found = false;

      // If row j contains element for column k, perform sum.
      // Otherwise, sum = a * _A[i,k]
      if (found) 
      {
  	if (_F.isZero(Faxpy.applyin(iter_j->second, iter_i->second)))
  	  _A[j].erase(iter_j++);
  	else
  	  iter_j++;

      } // if (found)
      else //if (!found)
  	_A[j].insert(iter_j,
  		     make_pair(k, _F.mul(value,a,iter_i->second)));
	  
    } // for ( iter_i = (_A[i]).begin(); iter_i != (_A[i]).end(); iter_i++)

    return;

  } // void sparsemat_base<sparse_associative_vector_tag>::addrow(...)
*/
  /* This implementation does not have to worry about any invalidated
   * iterators and references because the addition is not done inplace.
   * However, this means it is not as efficient since a new row has
   * to be created and then assigned to _A[j].
   */
  template <class Field, class Row>
  void sparsemat_base<Field, 
                      Row, 
		      vector_categories::sparse_associative_vector_tag>
  ::addrow(size_t i, size_t j,const element& a) 
  {
    // Check row indices and print error if they are out of range.
    if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
    {
      cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
	           << "        No row addition preformed." << endl << endl;
      return;
    }

    // Check to see if a is the zero Field element.
    // If so, no addition is performed.
    if (_F.isZero(a)) return;

    // Check to see if row i is empty.  If so, no addition is preformed.
    if ( (_A[i]).begin() == (_A[i]).end() ) return;

    // variables used in computation
    element value;
    _F.init(value, 0);
    LinBox::faxpy<Field> Faxpy(_F, a);
    row_type row;
    row_iter iter_i, iter_j(_A[j].begin());

    for (iter_i = _A[i].begin(); iter_i != _A[i].end(); iter_i++)
    {
      while( (iter_j != _A[j].end()) && (iter_j->first < iter_i->first) )
      {
	row.insert(*iter_j);
	iter_j++;
      } // while( (iter_j != _A[j].end()) && (iter_j->first < iter_i->first) )

      if( (iter_j != _A[j].end()) && (iter_j->first == iter_i->first) )
      {
	if (!_F.isZero(Faxpy.apply(value, iter_i->second, iter_j->second)))
	  row.insert(make_pair(iter_i->first, value));
	
	iter_j++;
      } // if( (iter_j != _A[j].end()) && (iter_j->first == iter_i->first) )
      else
	row.insert(make_pair(iter_i->first,
			     _F.mul(value, a, iter_i->second)));
      
    } // for (iter_i = _A[i].begin(); iter_i != _A[i].end(); iter_i++)

    while(iter_j != _A[j].end())
    {
      row.insert(*iter_j);
      iter_j++;
    } // while(iter_j != _A[j].end())
    
    _A[j] = row;

    return;
    
  } //  void sparsemat_base<sparse_associative_vector_tag>::addrow(...)

  // Input/Output Operators.

  template <class element>
  ostream& operator<<(ostream& os, pair< size_t, element > entry)
  {
  // Requires operator<<(ostream& element) which may not be provided.
    os << "(" << entry.first << ", " << entry.second << ")";
    return os;
  } // operator<< pair

  template <class Field, class Row>
  ostream& operator<<(ostream& os, const sparsemat_base<Field, Row>& A)
  { return A.write(os); }


  template <class Field, class Row>
  istream& operator>>(istream& is, sparsemat_base<Field, Row>& A)
  { return A.read(is); }

  /* Creates new sparse matrix.
   * Reads matrix from input stream.
   * If TRACE is defined, output includes all inputs.
   * @return pointer to new matrix in dynamic memory
   * @param  F field in which all arithmetic is done
   * @param  m row dimensions of matrix (defualt = 0)
   * @param  n column dimensions of matrix (default = 0)
   * @param  prompt  boolean for whether to prompt user for input
   *		      (default = true)
   * @param  is  istream from which to read input (default = cin)
   * @param  os  output stream to which to print output (default = cout)
   * @see sparsemat_base
   */
  template <class Field, class Row>
  sparsemat_base<Field, Row>* 
  newSparsemat(const Field& F, 
	       size_t m=0, 
	       size_t n=0,
	       bool prompt = true,
	       istream& is = cin, 
	       ostream& os = cout)
  {
    while ( (m <= 0) || (n <= 0) )
    {
      if (prompt)
  	cout << "What are the matrix's row and column dimenstions? ";
 
      is >> m >> n;

#ifdef TRACE
      os << endl << "The matrix has " << m << " rows and " << n << endl;
#endif
    } // if
 
    sparsemat_base<Field, Row>* A_ptr = new sparsemat_base<Field, Row>(F,m,n);

    if (prompt)
      cout << endl << "Input sparse matrix by entering row index, column"
  	   << endl << "index, and value.  Remember the matrix is indexed"
	   << endl << "starting at 0.  End with a row index of -1."
	   << endl;

    is >> (*A_ptr);

#ifdef TRACE
    os << endl << "The matrix contains the following elements: "
       << endl << (*A_ptr) << endl;
#endif

    return A_ptr;
  } // newSparsemat()

} // namespace LinBox

#endif // _SPARSEMAT_BASE_
