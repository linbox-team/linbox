/* -*- mode: c; style: linux -*- */

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

namespace LinBox
{
  template<class Field,class Vect>
    template<class Iterator1, class Iterator2>
    Iterator1& DenseMatrix<Field,Vect>::apply (Iterator1& in, Iterator2& outbegin, Iterator2& outend) const
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

  template<class Field, class Vect>
    template<class Vect1, class Vect2>
    Vect1& DenseMatrix<Field,Vect>::apply (Vect1& y, const Vect2& x) const
    {
      ConstColOfRowsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();

      for (p=colOfRowsBegin();p!=colOfRowsEnd();++p,++p_y)
	_VD.dotproduct (*p_y,*p, x);
      
      return y;
    }
  template<class Field, class Vect>
    template<class Vect1, class Vect2>
    Vect1& DenseMatrix<Field,Vect>::applyTranspose (Vect1& y, const Vect2& x) const
    {
      ConstRowOfColsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();
      
      for (p=rowOfColsBegin();p!=rowOfColsEnd();++p,++p_y)
	_VD.dotproduct (*p_y,*p, x);
      
      return y;
    }

   template<class Field, class Vect>
    template<class Iterator1, class Iterator2>
    Iterator1& DenseMatrix<Field, Vect>::applyTranspose (Iterator1& in, Iterator2& outbegin, Iterator2& outend) const
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

  template<class Field, class Vect>
    size_t DenseMatrix<Field,Vect>::rowdim (void) const 
    { return _rows;}
  
  template<class Field, class Vect>
  size_t DenseMatrix<Field,Vect>::coldim (void) const
    { return _cols; }
  
  // End, Blackbox interface
  
  /// entry access raw view.  Size m*x vector in C (row major) order.
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::RawIterator DenseMatrix<Field,Vect>::rawBegin()
    {return _rep.begin();}

  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::RawIterator DenseMatrix<Field,Vect>::rawEnd()
    {return _rep.end();}

  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstRawIterator DenseMatrix<Field,Vect>::rawBegin() const
    {return _rep.begin();}

  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstRawIterator DenseMatrix<Field,Vect>::rawEnd() const
    {return _rep.end();}

  template<class Field, class Vect>
    class DenseMatrix<Field,Vect>::Row
    {
    public:
      typedef RowIterator iterator;
      typedef ConstRowIterator const_iterator;
      Row(const Vector::iterator& p, size_t len):_rep(p),_len(len)
      {}
      
      Row(const Row& r):_rep(r._rep),_len(r._len)
      {}
      
      iterator begin()
	{return _rep;}
      
      const_iterator begin() const
	{return _rep;}
      
      iterator end()
	{return _rep+_len;}			
      
      const_iterator end() const
	{return _rep+_len;}
      
      Row& operator++()
	{
	  _rep+=_len;
	  return *this;
	}
      
      const Row& operator++() const
	{
	  const_cast<Row*>(this)->_rep=_rep+_len;
	  return *this;
	}

            
      Row operator++(int)
	{
	  Row tmp(*this);
	  _rep+=_len;
	  return tmp;
	}
      
      Row operator++(int) const
	{
	  Row tmp(*this);
	  const_cast<Row*>(this)->_rep=_rep+_len;
	  return tmp;
	}

      size_t size()
	{return _len;}
	  
      size_t size() const
	{ return _len;}

      bool operator!=(const Row& r)
      {return (_rep!=r._rep)||(_len!=r._len);}
      
      bool operator!=(const Row& r) const
      {return (_rep!=r._rep)||(_len!=r._len);}
      
    private:
      Vector::iterator _rep;
      size_t _len;  
    };

  template<class Field, class Vect>
    class DenseMatrix<Field,Vect>::ColOfRowsIterator
    {
    public:
      ColOfRowsIterator(const Vector::iterator& p =0,size_t len =0)
	:_row(p,len){}
      
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

      ConstRow* operator->() const
	{return &_row;}
      
      Row& operator*()
	{return _row;}

      ConstRow& operator*() const
	{return _row;}
      
      bool operator!=(const ColOfRowsIterator& c)
      {return _row!=c._row;}
       
      bool operator!=(const ColOfRowsIterator& c) const
      {return _row!=c._row;}
      
    private:
      Row _row;
    };

  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ColOfRowsIterator DenseMatrix<Field,Vect>::colOfRowsBegin()
    {
      return ColOfRowsIterator(_rep.begin(),_cols);
    }
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstColOfRowsIterator DenseMatrix<Field,Vect>::colOfRowsBegin() const
    {
      return ColOfRowsIterator(const_cast<DenseMatrix<Field,Vect>*>(this)->_rep.begin(),_cols);
    }
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ColOfRowsIterator DenseMatrix<Field,Vect>::colOfRowsEnd()
    {return ColOfRowsIterator(_rep.end(),_cols);}
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstColOfRowsIterator DenseMatrix<Field,Vect>::colOfRowsEnd() const
    {return ColOfRowsIterator(const_cast<DenseMatrix<Field,Vect>*>(this)->_rep.end(),_cols);}
  
  
  // row sequence of cols view
  template<class Field, class Vect>
    class DenseMatrix<Field,Vect>::ColIterator
    {
    public:
      ColIterator(Vector::iterator p =0, size_t step =0)
	{
	  _rep=p;
	  _step=step;
	}
      ColIterator(const ColIterator& ptr)
	{
	  _rep=ptr._rep;
	  _step=ptr._step;
	}

      ColIterator& operator=(const ColIterator& colp)
      {
	_rep=colp._rep;
	_step=colp._step;
	return *this;
      }
      
      const ColIterator& operator=(const ColIterator& colp) const
      {
	const_cast<ColIterator*>(this)->_rep=colp._rep;
	const_cast<ColIterator*>(this)->_step=colp._step;
	return *this;
      }
      ColIterator& operator++()
	{
	  _rep+=_step;
	  return *this;
	}
      const ColIterator& operator++ () const
	{
	  const_cast<ColIterator*>(this)->_rep=_rep+_step;
	  return *this;
	}
      
      ColIterator operator++(int)
	{
	  ColIterator tmp(*this);
	  ++this;
	  return tmp;
	}
      const ColIterator operator++(int) const
	{
	  ColIterator tmp(*this);
	  ++this;
	  return tmp;
	}
      
      Element& operator*()
	{
	  return *_rep;
	}
      
      const Element& operator*() const
	{
	  return *_rep;
	}
      
      bool operator!=(const ColIterator& p)
      { return (_rep!=p._rep)||(_step!=p._step);}

      bool operator!=(const ColIterator& p) const
      { return (_rep!=p._rep)||(_step!=p._step);}
	
    private:
      Vector::iterator _rep;
      size_t _step;
    };
  
  template<class Field, class Vect>
    class DenseMatrix<Field,Vect>::Col
    {
    public:
      typedef ColIterator iterator;
      typedef ConstColIterator const_iterator;
      Col(Vector::iterator p =0, size_t len =0, int stride =0):_rep(p),_len(len),_stride(stride)
	{}
      
      Col(const Col& r):_rep(r._rep),_len(r._len),_stride(r._stride)
	{}
      
      Col& operator=(const Col& colp)
      {
	_rep=colp._rep;
	_len=colp._len;
	_stride=colp._stride;
	return *this;
      }
      const Col& operator=(const Col& colp) const
      {
	const_cast<Col*>(this)->_rep=colp._rep;
	const_cast<Col*>(this)->_len=colp._len;
	const_cast<Col*>(this)->_stride=colp._stride;
	return *this;
      }
      iterator begin()
	{return ColIterator(_rep,_stride);}
      
      const_iterator begin() const
	{return ColIterator(_rep,_stride);}
      
      iterator end()
	{return ColIterator(_rep+_len*_stride,_stride);}
      
      const_iterator end() const
	{return ColIterator(_rep+_len*_stride,_stride);}
   
      Col& operator++()
	{
	  ++_rep;
	  return *this;
	}
      
      const Col& operator++() const
	{
	  ++const_cast<Col*>(this)->_rep;
	  return *this;
	}
      size_t size() 
	{ return _len;}

      size_t size() const
	{ return _len;}
      
      bool operator!=(const Col& r)
      { return (_rep!=r._rep)||(_len!=r._len)||(_stride!=r._stride);}
      
      bool operator!=(const Col& r) const
      { return (_rep!=r._rep)||(_len!=r._len)||(_stride!=r._stride);}
      
    private:
      Vector::iterator _rep;
      size_t _len;  
      int _stride;
    };
  
  template<class Field, class Vect>
    class DenseMatrix<Field,Vect>::RowOfColsIterator
    {
    public:
      RowOfColsIterator(Vector::iterator p =0,size_t len =0, int stride =0)
	:_col(p,len,stride){}
      
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
      
      Col operator*()
	{return _col;}
      
      ConstCol operator*() const
	{return _col;}
      
      bool operator!=(const RowOfColsIterator& c)
      {return _col!=c._col;}
      
      bool operator!=(const RowOfColsIterator& c) const
      {return _col!=c._col;}
      
    private:
      Col _col;
    };
 
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::RowOfColsIterator DenseMatrix<Field,Vect>::rowOfColsBegin()
    {
      return  DenseMatrix<Field,Vect>::RowOfColsIterator(_rep.begin(),_cols,_rows);
    }
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::RowOfColsIterator DenseMatrix<Field,Vect>::rowOfColsEnd()
    {
      return  DenseMatrix<Field,Vect>::RowOfColsIterator(_rep.begin()+_rows,_cols,_rows);
    }
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstRowOfColsIterator DenseMatrix<Field,Vect>::rowOfColsBegin() const
    {
      return  DenseMatrix<Field,Vect>::RowOfColsIterator(const_cast<DenseMatrix<Field,Vect>*>(this)->_rep.begin(),_cols,_rows);
    }

  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::ConstRowOfColsIterator DenseMatrix<Field,Vect>::rowOfColsEnd() const
    {
      return  DenseMatrix<Field,Vect>::RowOfColsIterator(const_cast<DenseMatrix<Field,Vect>*>(this)->_rep.begin()+_rows,_cols,_rows);
    }



  /** Set the entry at (i, j)
   * @param i Row number, 0...rowdim () - 1
   * @param j Column number 0...coldim () - 1
   * @param a_ij Element to set
   */
  template<class Field, class Vect>
    void DenseMatrix<Field,Vect>::setEntry (size_t i, size_t j, Element& a_ij) 
    { _rep[i*_cols+j] = a_ij; }
  
  template<class Field, class Vect>
    DenseMatrix<Field,Vect>::Element& DenseMatrix<Field,Vect>::getEntry (size_t i, size_t j, Element& a_ij) 
    { return a_ij = _rep[i*_cols+j]; }
		/** Read the matrix from an input stream
		 * @param file Input stream from which to read
		 */
  template<class Field, class Vect>
    void DenseMatrix<Field,Vect>::read (std::istream &file)
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
  template<class Field, class Vect>
  std::ostream& DenseMatrix<Field,Vect>::write(std::ostream &os = std::cout)
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
