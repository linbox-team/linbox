/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
  class DenseMatrix<Field>::ConstColOfRowsIterator
  {
  public:
    ConstColOfRowsIterator(const Vector::const_iterator& p,size_t len, size_t d)
      :_row(p,p+len),_dis(d){}
    
    ConstColOfRowsIterator() {}
    
    ConstColOfRowsIterator(const ConstColOfRowsIterator& colp)
      :_row(colp._row),_dis(colp._dis){}
    
    ConstColOfRowsIterator& operator=(const ConstColOfRowsIterator& colp)
    {
      _row=colp._row;
      _dis=colp._dis;
      return *this;
    }
   
    ConstColOfRowsIterator& operator++()
    {
      _row=ConstRow(_row.begin()+_dis,_row.end()+_dis);
      return *this;
    }
    
    ConstColOfRowsIterator  operator++(int)
    {
      ColOfRowsIterator tmp(*this);
      ++_row;
      return tmp;
    }

    ConstColOfRowsIterator& operator+(int i)
    {
      _row=ConstRow(_row.begin()+_dis*i,_row.end()+_dis*i);
      return *this;
    }

    ConstRow operator[] (int i) const
    { return ConstRow(_row.begin()+_dis*i,_row.end()+_dis*i);}
    
    ConstRow* operator->()
    {return &_row;}
    
    ConstRow& operator*()
    {return _row;}
    
    bool operator!=(const ConstColOfRowsIterator& c) const
    {return (_row.begin()!=c._row.begin())||(_row.end()!=c._row.end())||(_dis!=c._dis);}
    
  private:
    ConstRow _row;
    size_t _dis;
  };

  template<class Field>
  class DenseMatrix<Field>::ColOfRowsIterator
  {
  public:
    ColOfRowsIterator(const Vector::iterator& p,size_t len, size_t d)
      :_row(p,p+len),_dis(d){}
    
    ColOfRowsIterator() {}
    
    ColOfRowsIterator(const ColOfRowsIterator& colp)
      :_row(colp._row),_dis(colp._dis){}
    
    ColOfRowsIterator& operator=(const ColOfRowsIterator& colp)
    {
      _row=colp._row;
      _dis=colp._dis;
      return *this;
    }
    
    
    ColOfRowsIterator& operator++()
    {
      _row=Row(_row.begin()+_dis,_row.end()+_dis);
      return *this;
    }
    
    ColOfRowsIterator  operator++(int)
    {
      ColOfRowsIterator tmp(*this);
      ++_row;
      return tmp;
    }
    
    ColOfRowsIterator& operator+(int i)
    {
      _row=tRow(_row.begin()+_dis*i,_row.end()+_dis*i);
      return *this;
    }

    Row operator[] (int i) const
    { return Row(const_cast<Row&>(_row).begin()+_dis*i,const_cast<Row&>(_row).end()+_dis*i);}
      
    Row* operator->()
    {return &_row;}
    
    Row& operator*()
    {return _row;}
 
    bool operator!=(const ColOfRowsIterator& c) const
    {return (_row.begin()!=c._row.begin())||(_row.end()!=c._row.end())||(_dis!=c._dis);}
    
    operator ConstColOfRowsIterator()
    { return ConstColOfRowsIterator(_row.begin(),_row.size(),_dis); }
    
  private:
    Row _row;
    size_t _dis;
  };

  template<class Field>
  class DenseMatrix<Field>::ConstRowOfColsIterator
  {
  public:
    ConstRowOfColsIterator(Vector::const_iterator p,size_t stride, size_t len)
      :_col(ConstColIterator(p,stride), ConstColIterator(p+len*stride, stride)),_stride(stride){}
    
    ConstRowOfColsIterator(const ConstCol& col, size_t stride)
      :_col(col),_stride(stride){}

    ConstRowOfColsIterator() {}
    
    ConstRowOfColsIterator(const ConstRowOfColsIterator& rowp)
      :_col(rowp._col){}
    
    ConstRowOfColsIterator& operator=(const ConstRowOfColsIterator& rowp)
    {
      _col=rowp._col;
      _stride=rowp._stride;
      return *this;
    }
      
    ConstRowOfColsIterator& operator++()
    {
      _col=ConstCol(ConstColIterator(_col.begin().operator->()+1,_stride),ConstColIterator(_col.end().operator->()+1,_stride));
      return *this;
    }
    
    ConstRowOfColsIterator  operator++(int)
    {
      Col tmp(_col);
      this->operator++();
      return tmp;
    }
    
    ConstRowOfColsIterator& operator+(int i)
    { 
      _col=ConstCol(ConstColIterator(_col.begin().operator->()+i,_stride),
		    ConstColIterator(_col.end().operator->()+i,_stride));
      return *this;
    }


    ConstCol operator[](int i) const
    { return ConstCol(ConstColIterator(_col.begin().operator->()+i,_stride),
		      ConstColIterator(_col.end().operator->()+i,_stride)); }

    ConstCol* operator->()
    {return &_col;}
 
    ConstCol& operator*()
    {return _col;}
    
    bool operator!=(const ConstRowOfColsIterator& c) const
    {return (_col.begin()!=c._col.begin())||(_col.end()!=c._col.end());}
    
  private:
    ConstCol _col;
    size_t _stride;
  };

  template<class Field>
  class DenseMatrix<Field>::RowOfColsIterator
  {
  public:
    RowOfColsIterator(Vector::iterator p,size_t stride, size_t len)
      :_col(ColIterator(p,stride), ColIterator(p+len*stride, stride)),_stride(stride){}
    
    RowOfColsIterator() {}
    
    RowOfColsIterator(const RowOfColsIterator& rowp)
      :_col(rowp._col){}
    
    RowOfColsIterator& operator=(const RowOfColsIterator& rowp)
    {
      _col=rowp._col;
      _stride=rowp._stride;
      return *this;
    }
    
    const RowOfColsIterator& operator=(const RowOfColsIterator& rowp) const
    {
      const_cast<RowOfColsIterator*>(this)->_col=rowp._col;
      return *this;
    }
    
    RowOfColsIterator& operator++()
    {
      _col=Col(ColIterator(_col.begin().operator->()+1,_stride),ColIterator(_col.end().operator->()+1,_stride));
      return *this;
    }
    
    RowOfColsIterator  operator++(int)
    {
      Col tmp(_col);
      this->operator++();
      return tmp;
    }
        
    RowOfColsIterator& operator+(int i)
    { 
      _col=Col(ConstColIterator(_col.begin().operator->()+i,_stride),
	       ConstColIterator(_col.end().operator->()+i,_stride));
      return *this;
    }
    

    Col operator[](int i) const
    { return Col(ColIterator(const_cast<Col&>(_col).begin().operator->()+i,_stride),
		 ColIterator(const_cast<Col&>(_col).end().operator->()+i,_stride)); }
    
    Col* operator->()
    {return &_col;}
    
    Col& operator*()
    {return _col;}
    
    bool operator!=(const RowOfColsIterator& c) const
    { return (_col.begin()!=c._col.begin())||(_col.end()!=c._col.end());}
    
    operator ConstRowOfColsIterator()
    {
      ConstCol tmp;
      tmp=_col;
      return ConstRowOfColsIterator(tmp, _stride);
    }
    
  private:
    Col _col;
    size_t _stride;
  };

  template<class Field>
  template<class Vect1, class Vect2>
  Vect1& DenseMatrix<Field>::apply (Vect1& y, const Vect2& x) const
  {
    ConstColOfRowsIterator p;
    ConstRow::const_iterator pe;
    typename Vect1::iterator  p_y=y.begin();  
    typename Vect2::const_iterator p_x;

    for (p=colOfRowsBegin();p!=colOfRowsEnd();++p,++p_y)
      {
	_F.init(*p_y,0);

	for(pe=p->begin(),p_x=x.begin();pe!=p->end();++pe,++p_x)
	  _F.axpyin(*p_y,*pe,*p_x);
      }
    
    return y;

  }
 
 
  template<class Field>
  template<class Iterator1, class Iterator2>
  Iterator1& DenseMatrix<Field>::apply (Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const
  {
    linbox_check(coldim()==(outend-outbegin));
    ConstColOfRowsIterator rowp;
    Iterator2 p_out;
    ConstRow::const_iterator pe;
    for (rowp=colOfRowsBegin();rowp!=colOfRowsEnd();++rowp,++in)
      {
        _F.init(*in,0);
        for(pe=rowp->begin(),p_out=outbegin;pe!=rowp->end();++pe,++p_out)
          _F.axpyin(*in,*pe,*p_out);
      }
    
    return in;

  }
  
 
  template<class Field>
  template<class Vect1, class Vect2>
  Vect1& DenseMatrix<Field>::applyTranspose (Vect1& y, const Vect2& x) const
  {
    ConstRowOfColsIterator colp;
    ConstCol::const_iterator pe;
    typename Vect1::iterator  p_y=y.begin();  
    typename Vect2::const_iterator p_x;

    for (colp=rowOfColsBegin();colp!=rowOfColsEnd();++colp,++p_y)
      {
	_F.init(*p_y,0);
	for(pe=colp->begin(),p_x=x.begin();pe!=colp->end();++pe,++p_x)
	  _F.axpyin(*p_y,*pe,*p_x);
      }
    
    return y;
  }
  
  template<class Field>
  template<class Iterator1, class Iterator2>
  Iterator1& DenseMatrix<Field>::applyTranspose (Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const
  {
    linbox_check(rowdim()==(outend-outbegin));
    ConstRowOfColsIterator colp;
    Iterator2 p_out;
    ConstCol::const_iterator pe;
    for (colp=rowOfColsBegin();colp!=rowOfColsEnd();++colp,++in)
      {
        _F.init(*in,0);
        for(pe=colp->begin(),p_out=outbegin;pe!=colp->end();++pe,++p_out)
          _F.axpyin(*in,*pe,*p_out);
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
  DenseMatrix<Field>::ColOfRowsIterator DenseMatrix<Field>::colOfRowsBegin()
  {return ColOfRowsIterator(_rep.begin(),_cols,_cols);}
  template<class Field>
  DenseMatrix<Field>::ColOfRowsIterator DenseMatrix<Field>::colOfRowsEnd()
  {return ColOfRowsIterator(_rep.end(),_cols,_cols);}
  
  template<class Field>
  DenseMatrix<Field>::ConstColOfRowsIterator DenseMatrix<Field>::colOfRowsBegin() const
  { return ConstColOfRowsIterator(_rep.begin(),_cols,_cols);}  
  template<class Field>
  DenseMatrix<Field>::ConstColOfRowsIterator DenseMatrix<Field>::colOfRowsEnd() const
  {return ConstColOfRowsIterator(_rep.end(),_cols,_cols);}
  
 
  
  template<class Field>
  DenseMatrix<Field>::RowOfColsIterator DenseMatrix<Field>::rowOfColsBegin()
  { return  DenseMatrix<Field>::RowOfColsIterator(_rep.begin(),_cols,_rows);}
  template<class Field>
  DenseMatrix<Field>::RowOfColsIterator DenseMatrix<Field>::rowOfColsEnd()
  { return  DenseMatrix<Field>::RowOfColsIterator(_rep.begin()+_cols,_cols,_rows);}
  
  template<class Field>
  DenseMatrix<Field>::ConstRowOfColsIterator DenseMatrix<Field>::rowOfColsBegin() const
  { return  DenseMatrix<Field>::ConstRowOfColsIterator(_rep.begin(),_cols,_rows); }
  
  template<class Field>
  DenseMatrix<Field>::ConstRowOfColsIterator DenseMatrix<Field>::rowOfColsEnd() const
  { return  DenseMatrix<Field>::ConstRowOfColsIterator(_rep.begin()+_cols,_cols,_rows);}
  
  /** Set the entry at (i, j)
   * @param i Row number, 0...rowdim () - 1
   * @param j Column number 0...coldim () - 1
   * @param a_ij Element to set
   */
  template<class Field>
  void DenseMatrix<Field>::setEntry (size_t i, size_t j, const Element& a_ij) 
  { _rep[i*_cols+j] = a_ij; }
  
  template<class Field>
  DenseMatrix<Field>::Element& DenseMatrix<Field>::getEntry (Element& a_ij, size_t i, size_t j ) const
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
	_F.read (file, *p);
      }
  }
  
  /** Write the matrix to an output stream
   * @param os Output stream to which to write
   */
  template<class Field>
  std::ostream& DenseMatrix<Field>::write(std::ostream &os) const
  {
    ConstColOfRowsIterator p;
    for (p=colOfRowsBegin();p!=colOfRowsEnd();++p) 
      {
	ConstRow::const_iterator pe;
	for(pe=p->begin();pe!=p->end();++pe)
	  {
	    _F.write (os, *pe);
	    os << " ";
	  }
	os<<"\n";
      }		  
    os << endl;
    return os;
  }

  template<class Field>
  DenseMatrix<Field>::Row DenseMatrix<Field>::operator[] (size_t i)
  { return Row(_rep.begin()+i*_cols,_rep.begin()+i*_cols+_cols);}

  template<class Field>
  DenseMatrix<Field>::ConstRow DenseMatrix<Field>::operator[] (size_t i) const
  { return Row(_rep.begin()+i*_cols,_rep.begin()+i*_cols+_cols);}
}

#endif
