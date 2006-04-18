/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/zero-one.h
 *
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __ZERO_ONE_H
#define __ZERO_ONE_H

#include "linbox/integer.h"
//#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/field/modular.h"
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/quad-matrix.h>

// For STL pair in RawIndexIterator
#include <utility>
#include <iterator>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>

namespace LinBox
{
  bool revLexLess(const std::pair<size_t,size_t>& a, const std::pair<size_t,size_t> b)
  { return a.second < b.second || (b.second == a.second && a.first < b.first); }
  
  /** \brief Time and space efficient representation of sparse {0,1}-matrices.
   *  
   * A 0-1 matrix is a matrix with all 0's and 1's as entries.  
   * We're using a NAG-sparse format. 
   * Applies can be performed fast, using only additions.
   * When initalizing this class, you only need to build 2 arrays of equal length:
   * an array of the row indices for the non-zero (1's) entries, and an array of the column
   * indices for the non-zero (1's) entries.
   
   A {0, 1,-1} matrix can be effecively represented as the \ref Dif of two ZeroOne's.
   \ingroup blackbox
   */
  
  template<class _Field>
  class ZeroOne : public BlackboxInterface
  {
  public:
    friend class ZOQuad;
    
    typedef size_t Index;
    typedef ZeroOne<_Field> Self_t;
    typedef _Field Field;
    typedef typename _Field::Element Element;
    typedef std::vector<Index> IndexVector;
    typedef std::vector<IndexVector::iterator> PointerVector;
    //to denote by which way we sort our matrix
    enum howToSort { sortedByRow, sortedByCol };

    // DEFAULT CONSTRUCTOR, do nothing. Matrix will be uninitialized.
    ZeroOne(){};

    // basic constructor, can be used with subsequent read.
    ZeroOne(const Field& F): _F(F), sorted(sortedByRow){}

    // constructor for use by ZOQuad.  Needs work.
    ZeroOne
    (Field& F, IndexVector& index, PointerVector& indexP, Index coldim, howToSort sortedBy)
      : _F(F), _index(index), _indexP(indexP), _rowdim(indexP.size()-1), _coldim(coldim), sorted(sortedBy)
      //(Field& F, IndexVector& col, PointerVector& rowP, IndexVector& row, PointerVector& colP)
      //  : _F(F), _col(col), _rowP(rowP), _row(row), _colP(colP), _rowdim(rowP.size()-1), _coldim(cols), sorted(sortedByRow)
    {}

    /** The real constructor /todo give docs here
	assuming entries are sorted in lexicographic order by (row,col) pair.
    */
    ZeroOne
    (Field& F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz) 
      : _F(F), _rowdim(rows), _coldim(cols), sorted(sortedByRow)
    { 
      std::vector<std::pair<Index, Index> > indexPairs;
      for (Index i = 0; i < NNz; ++i, ++rowP, ++colP) 
	indexPairs.push_back(std::pair<Index,Index>(static_cast<Index>(*rowP), static_cast<Index>(*colP)));
      init(indexPairs);
    }

   //switching the way in which the matrix is sorted
    void switch_sort() const
    {
      Index dim;

      if( sorted == sortedByRow ) dim = _coldim;
      else dim = _rowdim;

      IndexVector temp[dim];
      for( PointerVector::iterator i = _indexP.begin(); i != _indexP.end() - 1; ++i )
	{
	  IndexVector::iterator j = *i;
	  for( ; j != *(i+1); ++j ) temp[*j].push_back( (Index)(i - _indexP.begin()) );
	}
      _index.clear(); _indexP.clear();
      std::back_insert_iterator < std::vector<Index> > colend( _index) ;
      for( size_t k = 0; k < dim; ++k )
	{
	  _indexP.push_back( _index.end() );
	  copy( temp[k].begin(), temp[k].end(), colend );
	}
      _indexP.push_back( _index.end() );
      if( sorted == sortedByRow ) sorted = sortedByCol;
      else sorted = sortedByRow;
      
      return;
    }

  protected: 
    
    void init(std::vector<std::pair<Index, Index> >& ip) 
    { 
      Index NNz = ip.size();
      sort(ip.begin(), ip.end());
      
      // set up _index
      for (Index i = 0; i < NNz; ++i) 
	_index.push_back(ip[i].second);
      
      // set up _indexP
      IndexVector::iterator p = _index.begin(); //how about if the first row of the matrix is empty
      _indexP.push_back(p);                     //we haven't met this case up to now
      
      std::vector<std::pair<Index, Index> >::iterator q = ip.begin();
      Index i = q->first;
      p++;q++; //start from the second place
      
      for (; q != ip.end(); ++q, ++p)
	if (i != q->first)
	  {
	    for (Index j = i; j < (q+1)->first; j++)//difference may be more than 1
	      _indexP.push_back(p);
	    i = q->first; //we should change i after the for loop, otherwise the 
	                  //for loop will not run at all
	  }
      _indexP.push_back(_index.end());
      
      /* keep another copy is not needed if we can switch sort between row and col
      // sort by cols first, then do the same as above
      sort(ip.begin(), ip.end(), revLexLess);
      
      // set up _row
      for (Index i = 0; i < NNz; ++i)
	_row.push_back(ip[i].first);
      
      // set up _colP
      p = _row.begin();
      _colP.push_back(p);
      
      q =ip.begin();
      i = q->second;
      p++;q++;
      
      for ( ; q != ip.end(); ++q, ++p)
	if (i != q->second)
	  {
	    for (Index j = i; j < (q+1)->second; j++)
	      _colP.push_back(p);
	    i = q->second; 
	  }
      _colP.push_back(_row.end());
      */

    }
  public:
    
    // Destructor, once again do nothing
    //~ZeroOne(){};
    
    /** \brief
     *
     * Uses one of the three
     * private utility functions. It calls the generalized utility function
     * _apply if there is no special ordering, _fyapply if there is C_ordering
     * or _fxapply if there is fortran_ordering
     */
    template<class OutVector, class InVector>
    OutVector& apply(OutVector& y, const InVector& x) const; // y = Ax;
    //OutVector& apply(OutVector& y, const InVector& x); // y = Ax;
    
    /** \brief
     *
     * Uses one of the three
     * private utility functions, in the manner described above.  Worthy of
     * note is the fact that applyTranspose works by passing the column
     * positions to the _apply functions as if they were rows, and row positions
     * as if they were columns, as if the matrix had been transposed.
     */
    
    template<class OutVector, class InVector>
    OutVector& applyTranspose(OutVector& y, const InVector& x) const; // y = ATx
    //OutVector& applyTranspose(OutVector& y, const InVector& x); // y = ATx
    
    size_t rowdim() const { return _rowdim; }
    
    size_t coldim() const { return _coldim; }      
    
    
    template<typename _Tp1>
    struct rebind 
    { 
      typedef ZeroOne<_Tp1> other;
      void operator() (other *& Ap,
		       const Self_t& A, 
		       const _Tp1& F) {
	Ap = new other(F, A._indexP, A._index, A._rowdim, A._coldim);
      }
    };
    
    /** Read the matrix from a stream in the JGD's SMS format
     *  @param is Input stream from which to read the matrix
     *  @return Reference to input stream 
     */
    std::istream &read (std::istream &is)
    {
      std::vector<std::pair<Index, Index> > indexPairs;
      Index r, c; 
      Element v;
      MatrixStream<Field> S(_F, is);
      
      while (S.nextTriple(r, c, v) ) 
	//indexPairs.push_back(std::pair<Index,Index>(r, c));
	indexPairs.push_back(std::pair<Index,Index>(static_cast<Index>(r), static_cast<Index>(c)));
      
      S.getRows(_rowdim); S.getColumns(_coldim);
      //_nnz = indexPairs.size();
      init(indexPairs);
	//init(indexPairs, _rowdim, _coldim);
      return is;
    }
    std::ostream& write(std::ostream& out =std::cout) const
    {
      out << "ZeroOne Matrix: _index.size() " << _index.size();
      out << ", _indexP.size() " << _indexP.size();
      out << ", _rowdim " << _rowdim;
      out << ", _coldim " << _coldim;
      return out; 
    }
    
    const Field& field() const { return _F; }
    
    /* Non blackbox function.  Tells the number of nonzero entries
     */
    size_t nnz() const { return _index.size(); };
    
    
  protected:
    
    Field _F; // The field used by this class
    
    /* _indexP is a pointer to an array of row indexes.  _colP is a pointer
     * to an array of column indexes. These two are the other arrays of a
     * NAGSparse format Matrix.  _rowdim and _coldim are the number of rows and
     * columns of the Matrix if it were in dense format.  _nnz is the Number of
     * Non-Zero elements in the Matrix.  It also happens to be the length of
     * the three NAGSparse arrays.
     *
     * Note: I changed their names since now how the matrix is sorted is 
     * decided by the member "sorted" of type howToSort and we no longer
     * keep two copies of zo matrix by having (_col, _rowP) and (_row, _colP)
     *
     * _index stores the real row or column index, and _indexP always points to 
     * some position in _index which is the beginning or a row or a column, row 
     * or column depending on how the matrix is sorted
     *
     *
     */

    mutable IndexVector _index; // The nnz indices sorted by row or by col
    mutable PointerVector _indexP; // the pointers to beginning of each row if sorted by row 
                           // and to beginning of each col if sorted by col

    //IndexVector _col; // The nnz column indices
    //PointerVector _rowP; // the _rowdim+1 pointers to beginning of row in _col

    /* Keep another copy of the Matrix for applyTranspose */
    //IndexVector _row; // The nnz row indices
    //PointerVector _colP; // the _coldim+1 pointers to beginning of col in _row

    Index _rowdim, _coldim;
    mutable howToSort sorted;

  }; //ZeroOne

       
}//End of LinBox

#include "zo.inl"

#endif // __ZERO_ONE_H
