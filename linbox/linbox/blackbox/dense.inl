/* -*- mode: c; style: linux -*- */

/* linbox/blackbox/dense-matrix1.h
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
  template<class Field>
    /* try later
    template<class Vect1, class Vect2>
    Vect1& DenseMatrix1<Field>::apply (Vect1& y, const Vect2& x) const
    */
    std::vector<typename Field::Element>& 
    DenseMatrix1<Field>::apply 
	(std::vector<typename Field::Element>& y, 
	 const std::vector<typename Field::Element>& x) const
    {
      /* */ typedef std::vector<typename Field::Element> Vect1; 
	    typedef Vect1 Vect2;
      ColOfRowsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();

      for (p=colOfRowsBegin();p!=colOfRowsEnd();++p,++p_y)
	_VD.dotprod (*p_y,*p, x);
      
      return y;
    }
  template<class Field>
    //template<class Vect1, class Vect2>
    //Vect1& DenseMatrix1<Field>::applyTranspose (Vect1& y, const Vect2& x) const
    std::vector<typename Field::Element>& 
    DenseMatrix1<Field>::applyTranspose 
	(std::vector<typename Field::Element>& y, 
	 const std::vector<typename Field::Element>& x) const
    {
      /* */ typedef std::vector<typename Field::Element> Vect1; 
	    typedef Vect1 Vect2;
      RowOfColsIterator p;
      typename Vect1::iterator  p_y;
      p_y=y.begin();
      
      for (p=rowOfColsBegin();p!=rowOfColsEnd();++p,++p_y)
	_VD.dotprod (*p_y,*p, x);
      
      return y;
    }
		
  template<class Field>
    size_t DenseMatrix1<Field>::rowdim (void) const 
    { return _rows;}
  
  template<class Field>
  size_t DenseMatrix1<Field>::coldim (void) const
    { return _cols; }
  
		// End, Blackbox interface

		/// entry access raw view.  Size m*x vector in C (row major) order.
  template<class Field>
    DenseMatrix1<Field>::RawIterator DenseMatrix1<Field>::rawBegin()
    {return _rep.begin();}

  template<class Field>
    DenseMatrix1<Field>::RawIterator DenseMatrix1<Field>::rawEnd()
    {return _rep.end();}

  template<class Field>
    DenseMatrix1<Field>::ConstRawIterator DenseMatrix1<Field>::rawBegin() const
    {return _rep.begin();}

  template<class Field>
    DenseMatrix1<Field>::ConstRawIterator DenseMatrix1<Field>::rawEnd() const
    {return _rep.end();}

  template<class Field>
    class DenseMatrix1<Field>::Row
    {
    public:
      typedef RowIterator iterator;
      typedef ConstRowIterator const_iterator;
      Row(RowIterator p, size_t len):_len(len)
      {
	_rep=p;	
      }
      
      Row(const Row& r):_len(r._len)
	{
	  _rep=r._rep;
	}
      
      iterator begin()
	{return _rep;}
      
      const_iterator begin() const
	{return _rep;}
      
      iterator end()
	{return _rep+_len;}			
      
      const_iterator end() const
	{return _rep+_len;}
      
	// does this belong here?
      Row& operator+(size_t step)
	{
	  _rep+=step;
	  return *this;
	}
      
      size_t size()
	{ return _len;}

      bool operator!=(const Row& r)
      {return (_rep!=r._rep)||(_len!=r._len);}
      
    private:
      RowIterator _rep;
      size_t _len;  
    };

  template<class Field>
    class DenseMatrix1<Field>::ColOfRowsIterator
    {
    public:
      ColOfRowsIterator(RowIterator p,size_t len)
	:_row(p,len),_len(len){}
      
      ColOfRowsIterator(const ColOfRowsIterator& colp)
	:_row(colp._row),_len(colp._len){}
      
      ColOfRowsIterator& operator++()
	{
	  _row += _len;
	  return *this;
	}
      
      const ColOfRowsIterator& operator++() const
	{
	  const_cast<ColOfRowsIterator*>(this)->_row+_len;
	  return *this;
	}
      
      ColOfRowsIterator  operator++(int)
	{
	  Row tmp(_row);
	  _row += _len;
	  return tmp;
	}
      
      const ColOfRowsIterator operator++(int)  const
	{
	  Row tmp(_row);
	  const_cast<ColOfRowsIterator*>(this)->_row+_len;
	  return tmp;
	}
      Row* operator->()
	{return &_row;}
      ConstRow* operator->() const
	{return &_row;}
      
      Row& operator*()
	{return _row;}
      ConstRow operator*()
	{return _row;}
      
      
      bool operator!=(const ColOfRowsIterator& c) const
      {return (_row!=c._row)||(_len!=c._len);}
      
    private:
      Row _row;
      size_t  _len;
    };

  template<class Field>
    DenseMatrix1<Field>::ColOfRowsIterator DenseMatrix1<Field>::colOfRowsBegin()
    {
      return ColOfRowsIterator(_rep.begin(),_cols);
    }
  
  template<class Field>
    DenseMatrix1<Field>::ConstColOfRowsIterator DenseMatrix1<Field>::colOfRowsBegin() const
    {
      return ColOfRowsIterator(_rep.begin(),_cols);
    }
  
  template<class Field>
    DenseMatrix1<Field>::ColOfRowsIterator DenseMatrix1<Field>::colOfRowsEnd()
    {return ColOfRowsIterator(_rep.end(),_cols);}
  
  template<class Field>
    DenseMatrix1<Field>::ConstColOfRowsIterator DenseMatrix1<Field>::colOfRowsEnd() const
    {return ColOfRowsIterator(_rep.end(),_cols);}
  
  
  // row sequence of cols view
  template<class Field>
    class DenseMatrix1<Field>::ColIterator
    {
    public:
      ColIterator(Vector::iterator p, size_t step)
	{
	  _rep=p;
	  _step=step;
	}
      ColIterator(const ColIterator& ptr)
	{
	  _rep=ptr._rep;
	  _step=ptr._step;
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
      
    private:
      Vector::iterator _rep;
      size_t _step;
    };
  
  template<class Field>
    class DenseMatrix1<Field>::Col
    {
    public:
      typedef ColIterator iterator;
      typedef ConstColIterator const_iterator;
      Col(ColIterator p, size_t len, int stride):_len(len),_stride(stride)
      {
	_rep=p;	
      }
      
      Col(const Col& r):_len(r._len),_stride(r.stride)
	{
	  _rep=r._rep;
	}
      
      iterator begin()
	{return ColIterator(_rep,_step);}
      
      const_iterator begin() const
	{return ColIterator(_rep,_step);}
      
      iterator end()
	{return ColIterator(_rep+_len*_step,_step);}			
      
      const_iterator end() const
	{return Col(_rep+_len*_step,_step);}
      
      Col& operator+(size_t stride)
	{
	  _rep+=stride;
	  return *this;
	}
      
      size_t size() 
	{ return _len;}

      bool operator!=(const Col& r)
      { return (_rep!=r._rep)||(_len!=r._len)||(_stride!=r._stride);}
      
    private:
      ColIterator _rep;
      size_t _len;  
      int _stride;
    };
  
  template<class Field>
    class DenseMatrix1<Field>::RowOfColsIterator
    {
    public:
      RowOfColsIterator(ColIterator p,size_t len, int stride)
	:_col(p,len,stride){}
      
      RowOfColsIterator(const RowOfColsIterator& rowp)
	:_col(rowp._col){}
      
      RowOfColsIterator& operator++()
	{
	  _col+1;
	  return *this;
	}
      
      const RowOfColsIterator& operator++() const
	{
	  const_cast<RowOfColsIterator*>(this)->_col+1;
	  return *this;
	}
      
      RowOfColsIterator  operator++(int)
	{
	   Col tmp(_col);
	  _col+1;
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
      
      bool operator!=(const RowOfColsIterator& c) const
      {return _col!=c._col;}
      
    private:
      Col _col;
    };
 
  template<class Field>
    DenseMatrix1<Field>::RowOfColsIterator DenseMatrix1<Field>::rowOfColsBegin()
    {
      return ColIterator(_rep.begin(),_cols,_rows);
    }
  
  template<class Field>
    DenseMatrix1<Field>::RowOfColsIterator DenseMatrix1<Field>::rowOfColsEnd()
    {
      return ColIterator(_rep.begin+_rows,_cols,_rows);
    }
  
  template<class Field>
    DenseMatrix1<Field>::ConstRowOfColsIterator DenseMatrix1<Field>::rowOfColsBegin() const
    {
      return ColIterator(_rep.begin(),_col,_rows);
    }

  template<class Field>
    DenseMatrix1<Field>::ConstRowOfColsIterator DenseMatrix1<Field>::rowOfColsEnd() const
    {
      return ColIterator(_rep.begin()+_rows,_col,_rows);
    }



  /** Set the entry at (i, j)
   * @param i Row number, 0...rowdim () - 1
   * @param j Column number 0...coldim () - 1
   * @param a_ij Element to set
   */
  template<class Field>
    void DenseMatrix1<Field>::setEntry (size_t i, size_t j, Element& a_ij) 
    { _rep[i*_cols+j] = a_ij; }
  
  template<class Field>
    DenseMatrix1<Field>::Element& DenseMatrix1<Field>::getEntry (size_t i, size_t j, Element& a_ij) 
    { return a_ij = _rep[i*_cols+j]; }
		/** Read the matrix from an input stream
		 * @param file Input stream from which to read
		 */
  template<class Field>
    void DenseMatrix1<Field>::read (std::istream &file)
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
  std::ostream& DenseMatrix1<Field>::write(std::ostream &os = std::cout)
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
