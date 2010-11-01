/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/zero-one.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong, -bds
 * Time-stamp: <22 Jun 10 17:34:01 Jean-Guillaume.Dumas@imag.fr> 
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_zero_one_H
#define __LINBOX_zero_one_H

#include "linbox/integer.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/modular.h"
#include <linbox/blackbox/blackbox-interface.h>

// For STL pair in RawIndexIterator
#include <utility>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>
#include <iostream>

namespace LinBox
{
  
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
  protected:  
    typedef size_t Index;
  public:
    typedef ZeroOne<_Field> Self_t;
    typedef _Field Field;
    typedef typename _Field::Element Element;
    
    // Default constructor, do nothing.
    ZeroOne(const Field& F);
    // The real constructor /todo give docs here
    ZeroOne(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort = false, bool colSort = false);
    // Destructor, once again do nothing
    ~ZeroOne();
   
    /** \brief
     *
     * Uses one of the three
     * private utility functions. It calls the generalized utility function
     * _apply if there is no special ordering, _fyapply if there is C_ordering
     * or _fxapply if there is fortran_ordering
     */
    template<class OutVector, class InVector>
    OutVector& apply(OutVector& y, const InVector& x) const // y = Ax;
    { return applySpecialization(y,x,getType(_F)); }
    /** \brief
     *
     * Uses one of the three
     * private utility functions, in the manner described above.  Worthy of
     * note is the fact that applyTranspose works by passing the column
     * positions to the _apply functions as if they were rows, and row positions
     * as if they were columns, as if the matrix had been transposed.
     */
    
    template<class OutVector, class InVector>
    OutVector& applyTranspose(OutVector& y, const InVector& x) const // y = ATx
    { return applyTransposeSpecialization(y,x,getType(_F));}
    
    size_t rowdim() const { return _rows; }
    
    size_t coldim() const { return _cols; }      


    template<typename _Tp1>
    struct rebind 
    { 
      typedef ZeroOne<_Tp1> other;
      void operator() (other & Ap,
		       const Self_t& A, 
		       const _Tp1& F) {
              // ZeroOne does not store any field element
      }
    };

    template<typename _Tp1>
    ZeroOne(const ZeroOne<_Tp1>& Z, const Field& F) :
            _F(F), 
            _rows(Z.rowdim()), _cols(Z.coldim()), _nnz(Z.nnz()), 
            _rowP(new Index[Z.nnz()]), _colP(new Index[Z.nnz()]), 
            _rowSort(Z.isRowSorted()), _colSort(Z.isColSorted()) {

            Index * rowit = _rowP;
            Index * colit = _colP;

            for(typename ZeroOne<_Tp1>::RawIndexIterator it = Z.indexBegin();
                it != Z.indexEnd(); ++it,++rowit,++colit) {
                *rowit = (*it).first;
                *colit = (*it).second;
            }
    }
    
    /** RawIterator class.  Iterates straight through the values of the matrix
     */
    class RawIterator;

    RawIterator rawBegin();
    RawIterator rawEnd();
    const RawIterator rawBegin() const;
    const RawIterator rawEnd() const;

     /** RawIndexIterator - Iterates through the i and j of the current element
   * and when accessed returns an STL pair containing the coordinates
   */
    class RawIndexIterator;
    RawIndexIterator indexBegin();
    const RawIndexIterator indexBegin() const;
    RawIndexIterator indexEnd();
    const RawIndexIterator indexEnd() const;

    /** Read the matrix from a stream in the JGD's SMS format
     *  @param is Input stream from which to read the matrix
     *  @return Reference to input stream 
     */
    std::istream &read (std::istream &is){
      size_t i, j, k, m, n;
      
      char buf[80];
      buf[0]=0;
      is.getline (buf, 80);
      std::istringstream str (buf);	
      str >> m >> n >> k;
      _rows = m;
      _cols = n;
	  std::vector<size_t> rowP, colP;
      size_t x;
      while (is >> i >> j >> x) {
		if (i == 0 || i == (size_t) -1) break;
		if (x == 1UL) {
	  		rowP.push_back(i-1);
	  		colP.push_back(j-1);
		}
	  }
	  _nnz = rowP.size();
	  _rowP = new size_t[_nnz];
	  _colP = new size_t[_nnz];
	  copy(rowP.begin(), rowP.end(), _rowP);
	  copy(colP.begin(), colP.end(), _colP);
      return is;
    }
 
    std::ostream& write(std::ostream& out =std::cout)
    {
      size_t* i=_rowP;
      size_t* j=_colP;
      std::cout<<"Row dim: "<<rowdim()
	       <<" Col dim: "<<coldim()
	       <<" Total nnz: "<<nnz()<<"\n";
      for(;i<_rowP+nnz();++i,++j)
	std::cout<<*i<<" "<<*j<<"\n";     
	return out;
    }

    const Field& field() const { return _F; }

    bool isRowSorted() const { return _rowSort; }
    bool isColSorted() const { return _colSort; }
              
    size_t nnz() const { return _nnz; };
    

  protected:
   

    Field _F; // The field used by this class
    
    /* A temporary element used for initalization for the rawBegin() and
     * rawEnd() methods of the ZeroOne class.  Is used to initalize a 1
     * so that the RawIterator returned stores a 1 
     */

    Element _tmp; 
    
    
    /* _rowP is a pointer to an array of row indexes.  _colP is a pointer
     * to an array of column indexes. These two are the other arrays of a
     * NAGSparse format Matrix.  _rows and _cols are the number of rows and
     * columns of the Matrix if it were in dense format.  _nnz is the Number of
     * Non-Zero elements in the Matrix.  It also happens to be the length of
     * the three NAGSparse arrays.
     */
    
    Index _rows, _cols, _nnz;
    mutable Index* _rowP ;
	mutable Index* _colP;
    mutable bool _rowSort ;
	mutable bool _colSort; // status flags for sorting state          
    bool dynamic;

     /* Non blackbox function.  Tells the number of nonzero entries
     */
    void rowSort() const;
    void colSort() const;

    void _qsort(size_t start, size_t endp1, int &mode) const; // QuickSort function for when there is no sorting
    size_t _part( size_t start, size_t endp1, int &mode) const; // Partition for quicksort

  private:

    class FieldType {};
    class NormField : public FieldType {};
    class Mod32Field : public FieldType {}; 

    template<class F>
    NormField getType(const F &  f) const
    {
      return NormField();
    }
        
    Mod32Field getType(const Modular<uint32> &) const
    {
      return Mod32Field();
    }

    template<class OutVector, class InVector>
    OutVector& applySpecialization(OutVector &, const InVector &,const NormField& ) const;
    template<class OutVector, class InVector>
    OutVector& applySpecialization(OutVector &, const InVector &, const Mod32Field& )const;
    template<class OutVector, class InVector>
    OutVector& applyTransposeSpecialization(OutVector &, const InVector &,const NormField& ) const;
    template<class OutVector, class InVector>
    OutVector& applyTransposeSpecialization(OutVector &, const InVector &, const Mod32Field& )const;

  }; //ZeroOne

} //LinBox
       
#include "linbox/blackbox/zero-one.inl"

#endif // __LINBOX_zero_one_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
