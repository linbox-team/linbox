/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/dense.h
 *
 * evolved from dense-matrix.h by -bds, Zhendong Wan
 */

#ifndef __DENSE_MATRIX_C
#define __DENSE_MATRIX_C

#include <iostream>
#include <vector>
#include <fstream>
#include "linbox/blackbox/dense-matrix1.h"
#include "linbox/util/debug.h"

namespace LinBox
{
  template<class Field>
    template<class Iterator1, class Iterator2>
    Iterator1& DenseMatrix<Field>::apply (Iterator1& in, Iterator2& outbegin, Iterator2& outend) const
    {
      ConstColOfRowsIterator colp;
      Iterator2 p_out;
      p_y=y.begin();
      RowIterator prow;
      for (colp=colOfRowsBegin();colp!=colOfRowsEnd();++colp,++in)
	{
	  _field.init(*in,0);
	  for(prow=colp->begin(),p_out=outbegin;prow!=colp->end();++prow,++p_out)
	    _field.axpyin(*in,*prow,*p_out);
	}
	      
      return in;
    }

  template<class Field>
    template<class Vect1, class Vect2>
    Vect1& DenseMatrix<Field>::apply (Vect1& y, const Vect2& x) const
    {
      ConstColOfRowsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();

      for (p=colOfRowsBegin();p!=colOfRowsEnd();++p,++p_y)
	_VD.dotproduct (*p_y,*p, x);
      
      return y;
    }
  template<class Field>
    template<class Vect1, class Vect2>
    Vect1& DenseMatrix<Field>::applyTranspose (Vect1& y, const Vect2& x) const
    {
      ConstRowOfColsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();
      
      for (p=rowOfColsBegin();p!=rowOfColsEnd();++p,++p_y)
	_VD.dotproduct (*p_y,*p, x);

      return y;
    }

   template<class Field>
    template<class Iterator1, class Iterator2>
    Iterator1& DenseMatrix<Field>::applyTranspose (Iterator1& in, Iterator2& outbegin, Iterator2& outend) const
     {
       ConstRowOfColsIterator rowp;
      Iterator2 p_out;
      p_y=y.begin();
      ColIterator pcol;
      for (rowp=rowOfColsBegin();rowp!=rowOfColsEnd();++rowp,++in)
	{
	  _field.init(*in,0);
	  for(pcol=rowp->begin(),p_out=outbegin;pcol!=rowp->end();++prow,++p_out)
	    _field.axpyin(*in,*prow,*p_out);
	}
	      
      return in;
     }

  template<class Field>
    size_t DenseMatrix<Field>::rowdim (void) const 
    { return _rows;}
  
  template<class Field>
  size_t DenseMatrix<Field>::coldim (void) const
    { return _cols; }
  
  // End, Blackbox interface
  
  /// entry access raw view.  Size m*x vector in C (row major) order.
  template<class Field>
    DenseMatrix<Field>::RawIterator DenseMatrix<Field>::rawBegin()
    {return _rep.begin();}

  template<class Field>
    DenseMatrix<Field>::RawIterator DenseMatrix<Field>::rawEnd()
    {return _rep.end();}

  template<class Field>
    DenseMatrix<Field>::ConstRawIterator DenseMatrix<Field>::rawBegin() const
    {return _rep.begin();}

  template<class Field>
    DenseMatrix<Field>::ConstRawIterator DenseMatrix<Field>::rawEnd() const
    {return _rep.end();}
 	  
  template<class Field>
    class DenseMatrix<Field>::ColOfRowsIterator
    {
    public:
      ColOfRowsIterator(const Vector::iterator& p,size_t len =0)
	:_row(p,p+len){}

      ColOfRowsIterator() {}
      
      ColOfRowsIterator(const ColOfRowsIterator& colp)
	:_row(colp._row){}
      
      ColOfRowsIterator& operator=(const ColOfRowsIterator& colp)
      {
	_row=colp._row;
	return *this;
      }

      const ColOfRowsIterator& operator=(const ColOfRowsIterator& colp) const
      {
	const_cast<ColOfRowsIterator*>(this)->_row=colp._row;
	return *this;
      }
      
      ColOfRowsIterator& operator++()
	{
	  ++_row;
	  return *this;
	}
      
      const ColOfRowsIterator& operator++() const
	{
	  ++_row;
	  return *this;
	}
      
      ColOfRowsIterator  operator++(int)
	{
	  ColOfRowsIterator tmp(*this);
	  ++_row;
	  return tmp;
	}
      
      ColOfRowsIterator operator++(int) const
	{
	  ColOfRowsIterator tmp(*this);
	  ++_row;
	  return tmp;
	}

      Row* operator->()
	{return &_row;}

      const Row* operator->() const
	{return &_row;}
      
      Row& operator*()
	{return _row;}

      const Row& operator*() const
	{return _row;}
      
      bool operator!=(const ColOfRowsIterator& c)
      {return _row!=c._row;}
       
      bool operator!=(const ColOfRowsIterator& c) const
      {return _row!=c._row;}
      
    private:
      Row _row;
    };
   
  template<class Field>
    DenseMatrix<Field>::ColOfRowsIterator DenseMatrix<Field>::colOfRowsBegin()
    {
      return ColOfRowsIterator(_rep.begin(),_cols);
    }
  
  template<class Field>
    DenseMatrix<Field>::ConstColOfRowsIterator DenseMatrix<Field>::colOfRowsBegin() const
    {
      return ColOfRowsIterator(const_cast<DenseMatrix<Field>*>(this)->_rep.begin(),_cols);
    }
  
  template<class Field>
    DenseMatrix<Field>::ColOfRowsIterator DenseMatrix<Field>::colOfRowsEnd()
    {return ColOfRowsIterator(_rep.end(),_cols);}
  
  template<class Field>
    DenseMatrix<Field>::ConstColOfRowsIterator DenseMatrix<Field>::colOfRowsEnd() const
    {return ColOfRowsIterator(const_cast<DenseMatrix<Field>*>(this)->_rep.end(),_cols);}

  template<class Field>
    class DenseMatrix<Field>::RowOfColsIterator
    {
    public:
      RowOfColsIterator(Vector::iterator p,size_t stride =0, size_t len =0)
	:_col(ColIterator(p,stride), ColIterator(p+len*stride, stride)){}

      RowOfColsIterator() {}
      
      RowOfColsIterator(const RowOfColsIterator& rowp)
	:_col(rowp._col){}
      
      RowOfColsIterator& operator=(const RowOfColsIterator& rowp)
      {
	_col=rowp._col;
	return *this;
      }

      const RowOfColsIterator& operator=(const RowOfColsIterator& rowp) const
      {
	const_cast<RowOfColsIterator*>(this)->_col=rowp._col;
	return *this;
      }

      RowOfColsIterator& operator++()
	{
	  ++_col;
	  return *this;
	}
      
      const RowOfColsIterator& operator++() const
	{
	  ++_col;
	  return *this;
	}
      
      RowOfColsIterator  operator++(int)
	{
	   Col tmp(_col);
	  ++_col;
	  return tmp;
	}
      
      const RowOfColsIterator operator++(int)  const
	{
	  Col tmp(_col);
	  const_cast<RowOfColsIterator*>(this)->_col+1;
	  return tmp;
	}
      
      Col* operator->()
	{return &_col;}
      ConstCol* operator->() const
	{return &_col;}
      
      Col& operator*()
	{return _col;}
      
      const Col& operator*() const
	{return _col;}
      
      bool operator!=(const RowOfColsIterator& c)
      {return _col!=c._col;}
      
      bool operator!=(const RowOfColsIterator& c) const
      {return _col!=c._col;}
      
    private:
      Col _col;
    };
 
  template<class Field>
    DenseMatrix<Field>::RowOfColsIterator DenseMatrix<Field>::rowOfColsBegin()
    {
      return  DenseMatrix<Field>::RowOfColsIterator(_rep.begin(),_cols,_rows);
    }
  
  template<class Field>
    DenseMatrix<Field>::RowOfColsIterator DenseMatrix<Field>::rowOfColsEnd()
    {
      return  DenseMatrix<Field>::RowOfColsIterator(_rep.begin()+_cols,_cols,_rows);
    }
  
  template<class Field>
    DenseMatrix<Field>::ConstRowOfColsIterator DenseMatrix<Field>::rowOfColsBegin() const
    {
      return  DenseMatrix<Field>::RowOfColsIterator(const_cast<DenseMatrix<Field>*>(this)->_rep.begin(),_cols,_rows);
    }

  template<class Field>
    DenseMatrix<Field>::ConstRowOfColsIterator DenseMatrix<Field>::rowOfColsEnd() const
    {
      return  DenseMatrix<Field>::RowOfColsIterator(const_cast<DenseMatrix<Field>*>(this)->_rep.begin()+_cols,_cols,_rows);
    }



  /** Set the entry at (i, j)
   * @param i Row number, 0...rowdim () - 1
   * @param j Column number 0...coldim () - 1
   * @param a_ij Element to set
   */
  template<class Field>
    void DenseMatrix<Field>::setEntry (size_t i, size_t j, const Element& a_ij) 
    { _rep[i*_cols+j] = a_ij; }
  
  template<class Field>
    DenseMatrix<Field>::Element& DenseMatrix<Field>::getEntry (size_t i, size_t j, Element& a_ij) 
    { return a_ij = _rep[i*_cols+j]; }
		/** Read the matrix from an input stream
		 * @param file Input stream from which to read
		 */
  template<class Field>
    void DenseMatrix<Field>::read (std::istream &file)
    {
      RawIterator p;
      for (p=begin();p!=end();++p) 
	{
	  file.ignore (1);
	  _VD.read (file, *p);
	}
    }
  
  /** Write the matrix to an output stream
   * @param os Output stream to which to write
   */
  template<class Field>
  std::ostream& DenseMatrix<Field>::write(std::ostream &os)
    {
      RawIterator p;
      for (p=begin();p!=end();++p) 
	{
	  _VD.write (os, *p);
	  os << " ";
	}		  
      os << endl;
      return os;
    }
}

#endif
