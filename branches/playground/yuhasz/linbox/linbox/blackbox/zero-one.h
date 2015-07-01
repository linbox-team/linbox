/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong 
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __ZERO_ONE_H
#define __ZERO_ONE_H

#include "linbox/integer.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/modular.h"

// For STL pair in RawIndexIterator
#include <utility>
#include <vector> // For vectors in _col2row and _row2col
#include <cstdlib> // For randomness in randomized quicksort
#include <ctime>

using std::vector;

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>

using std::istream;
using std::ostream;
using std::string;

#endif


// home of linbox functionality
namespace LinBox
{
  
  /*- BlackBox representation of the 0-1's matrix
   *  
   * A 0-1's matrix is a matrix with all 0's and 1's as entries.  In
   * a nag-spasre format, applies could be performed with lightening speed
   * When initalizing this class, you only need to build 2 arrays of equal length:
   * an array of the row indices for the non-zero (1's) entries, and an array of the column
   * indices for the non-zero (1's) entries.
   */
  
  template<class Field>
  class ZeroOneBase
  {
  protected:  
    typedef typename Field::Element Element;
    typedef size_t Index;
    typedef LinBox::uint32 uint32;
    typedef LinBox::uint64 uint64;
  public:
    
    // Default constructor, do nothing.
    ZeroOneBase();
    // The real constructor
    ZeroOneBase(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort = false, bool colSort = false);
    // Destructor, once again do nothing
    ~ZeroOneBase();
   
#ifdef __LINBOX_XMLENABLED
	ZeroOneBase(Reader &);
	ZeroOneBase(const ZeroOneBase<Field>&);
#endif 
     
    /** Apply function.
     * Take constant vector x and
     * vector y, and perform the calculation y = Ax.  Uses one of the three
     * private utility functions. It calls the generalized utility function
     * _apply if there is no special ordering, _fyapply if there is C_ordering
     * or _fxapply if there is fortran_ordering
     */
    template<class OutVector, class InVector>
    OutVector& apply(OutVector& y, const InVector& x) const // y = Ax;
    { return applySpecialization(y,x,getType(_F)); }
    /** ApplyTranspose function. Take constant vector x and
     * vector y, and perform the calculation y = ATx.  Uses one of the three
     * private utility functions, in the manner described above.  Worthy of
     * note is the fact that applyTranspose works by passing the column
     * positions to the _apply functions as if they were rows, and row positions
     * as if they were columns, as if the matrix had been transposed.
     */
    
    template<class OutVector, class InVector>
    OutVector& applyTranspose(OutVector& y, const InVector& x) const // y = ATx
    { return applyTransposeSpecialization(y,x,getType(_F));}
    /** BlackBoxArchetype rowdim function.  Passes back the number of rows of
     * the matrix represented.  Note that's the number of rows of the matrix
     * as if it were in dense format, not in it's actual representation here.
     */
    
    size_t rowdim() const;
    
    /** BlackBoxArchetype coldim function.  Passes back the number of columns
     * of the matrix represented.  Not much more to say about this.
     */
    
    size_t coldim() const;      


    
    
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
    mutable Index* _rowP, *_colP;
    mutable bool _rowSort, _colSort; // status flags for sorting state          
    bool dynamic;

     /* Non blackbox function.  Tells the number of nonzero entries
     */
    size_t nnz() const;
    
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

  };


  /// Time and space efficient representation of sparse \{0,1}-matrices.
  template<class _Field>
  class ZeroOne : public ZeroOneBase<_Field>  {
	 
	  typedef typename ZeroOneBase<_Field>::Index Index;
  public:
	  typedef _Field Field;
	  typedef typename Field::Element Element;
	  
    ZeroOne(){}
    // The real constructor
    ZeroOne(const Field& F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, 
	    bool rowSort = false, bool colSort = false) 
    : ZeroOneBase<Field>(F,rowP,colP,rows,cols,NNz,rowSort,colSort)
    {}

    // Destructor, once again do nothing
    virtual ~ZeroOne() {};

#ifdef __LINBOX_XMLENABLED

	  ZeroOne(Reader &R);
	  ZeroOne(const ZeroOne<Field, Vector>&);
#endif


    template<class OutVector, class InVector>
    OutVector& apply(OutVector& y, const InVector& x) const
    {
      return ZeroOneBase<Field>::apply(y,x);
    }

    
    template<class OutVector, class InVector>
    OutVector& applyTranspose(OutVector&y, const InVector& x) const
    {
      return ZeroOneBase<Field>::apply(y,x);
    }

    virtual size_t coldim() const
    {
      return ZeroOneBase<Field>::coldim();
    }

    virtual size_t rowdim() const
    {
      return ZeroOneBase<Field>::rowdim();
    }

#ifndef __LINBOX_XMLENABLED

    std::ostream& write(std::ostream& out =std::cout)
    {
      size_t* i=_rowP;
      size_t* j=_colP;
      std::cout<<"Row dim: "<<rowdim()
	       <<" Col dim: "<<coldim()
	       <<" Total nnz: "<<nnz()<<"\n";
      for(;i<_rowP+nnz();++i,++j)
	std::cout<<*i<<" "<<*j<<"\n";     
    }
#else 
	ostream &write(ostream &) const;
	bool toTag(Writer &) const;

#endif

    const Field& field() const
    {
      return _F;
    }

  };
}

//End of namespace LinBox
 
       
#include "linbox/blackbox/zero-one.inl"


#endif // __ZERO_ONE_H
