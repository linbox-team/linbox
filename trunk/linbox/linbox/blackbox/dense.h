/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/dense.h
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 */

#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H

#include <iostream>
#include <vector>
#include <fstream>

#include "linbox/blackbox/archetype.h"
#include "linbox/field/dense-vector-domain.h"
#include "linbox/vector/subiterator.h"
#include "linbox/vector/subvector.h"

namespace LinBox
{
  /** Blackbox dense matrix template. This is a class of dense matrices
   * templatized by the {@link Fields field} in which the elements
   * reside. The matrix is stored as a one dimensional STL vector of
   * the elements, by rows. The interface provides for iteration
   * over rows and over columns.
   *
   * The class also conforms to the {@link Archetypes archetype} for
   * \Ref{Blackbox Matrices}.
   *
   * @param Field \Ref{LinBox} field
   */
  
  template <class Field>
  class DenseMatrix
    : public BlackboxArchetype<std::vector<typename Field::Element> >
  {
  public:
    typedef typename Field::Element        Element;
    typedef std::vector<Element>           Vector;
    typedef std::vector<Element>::iterator pointer;
    
    /** Constructor.
     * @param  F the field of entries; passed so that arithmetic may be done on elements. 
     * @param  m  row dimension
     * @param  n  column dimension
     */
    DenseMatrix (Field &F, size_t m, size_t n)
      : _F (F), _rep (m*n), _rows(m), _cols(n)
    {}
    
    
    /** Constructor.
     * @param  F the field of entries; passed so that arithmetic may be done on elements. 
     * @param  m  row dimension
     * @param  n  column dimension
     * @para iter, random iterator
     */
    template<class RandIter>
    DenseMatrix (Field &F, size_t m, size_t n,RandIter& iter )
      : _F (F), _rep (m*n), _rows(m), _cols(n)
    {
      for(Vector::iterator p=_rep.begin();p!=_rep.end();++p)
	iter.random(*p);
    }
    /** Copy constructor
     */
    DenseMatrix (const DenseMatrix &M)
      : _F (M._F), _rep (M._rep),_rows(M._rows), _cols(M._cols)
    {}
    
    /// Blackbox interface
    
    BlackboxArchetype<Vector> *clone () const 
    { return new DenseMatrix<Field> (*this);}
    
    template<class Vect1, class Vect2>
    Vect1& apply (Vect1& y, const Vect2& x) const;
    
    Vector& apply (Vector &y, const Vector &x) const
    {
      return apply<Vector,Vector>(y,x);
    }
    
    template<class Iterator1, class Iterator2 >
    Iterator1& apply( Iterator1 in, 
		      const Iterator2& outbegin, 
		      const Iterator2& outend) const;
    template<class Vect1, class Vect2>
    Vect1& applyTranspose (Vect1& y, const Vect2& x) const;
    
    Vector& applyTranspose (Vector& y, const Vector& x) const
    {
      return applyTranspose<Vector,Vector>(y,x);
    }
    
    template<class Iterator1, class Iterator2>
    Iterator1& applyTranspose (Iterator1 in, 
			       const Iterator2& outbegin, 
			       const Iterator2& outend) const;
    size_t rowdim (void) const;
    
    size_t coldim (void) const;
    
    
    /// entry access raw view.  Size m*x vector in C (row major) order.
    typedef Vector::iterator RawIterator;
    typedef Vector::const_iterator ConstRawIterator;
    
    RawIterator rawBegin();		  
    RawIterator rawEnd();
    
    ConstRawIterator rawBegin() const;
    
    ConstRawIterator rawEnd() const;
    
    // col sequence of rows view
    typedef Vector::iterator RowIterator;
    typedef Vector::const_iterator ConstRowIterator;
    
    class ColOfRowsIterator;
    
    typedef Subvector<DenseMatrix<Field>::RowIterator> Row;  
    typedef Subvector<DenseMatrix<Field>::ConstRowIterator> ConstRow;
     
    class ColOfRowsIterator;    
    class ConstColOfRowsIterator;

    ColOfRowsIterator colOfRowsBegin();  
    ColOfRowsIterator colOfRowsEnd();
 
    ConstColOfRowsIterator colOfRowsBegin() const;        
    ConstColOfRowsIterator colOfRowsEnd() const;
    

    typedef Subiterator<Vector::iterator> ColIterator;
    typedef Subiterator<Vector::const_iterator> ConstColIterator;

    typedef Subvector<ColIterator> Col;
    typedef Subvector<ConstColIterator> ConstCol;
    
    class RowOfColsIterator;
    class ConstRowOfColsIterator;
    
    RowOfColsIterator rowOfColsBegin();
    RowOfColsIterator rowOfColsEnd();
    
    ConstRowOfColsIterator rowOfColsBegin() const;    
    ConstRowOfColsIterator rowOfColsEnd() const;
    
    /** Set the entry at (i, j)
     * @param i Row number, 0...rowdim () - 1
     * @param j Column number 0...coldim () - 1
     * @param a_ij Element to set
     */
    void setEntry (size_t i, size_t j, const Element& a_ij) ;
    
    Element& getEntry (size_t i, size_t j, Element& a_ij) ;
    
    /** Read the matrix from an input stream
     * @param file Input stream from which to read
     */
    void read (std::istream &file);
    
    /** Write the matrix to an output stream
     * @param os Output stream to which to write
     */
    std::ostream &write(std::ostream &os = std::cout) const;
 
    const Field& field() const
    { return _F;}
    
  protected:
    
    std::vector<Element>                 _rep;
    Field                                _F;
    size_t _rows, _cols;
  };
}
#include "dense-matrix1.C"

#endif
