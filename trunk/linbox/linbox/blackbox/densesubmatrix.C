#ifndef DENSE_SUB_MATRIX_C
#define DENSE_SUB_MATRIX_C

#include <linbox/util/debug.h>
#include <linbox/blackbox/dense-matrix1.h>
#include <linbox/blackbox/densesubmatrix.h>

namespace LinBox
{
 
  template<class Field>
  DenseSubMatrix<Field>::DenseSubMatrix<Field>(DenseMatrix<Field>* _M, size_t _beg_row, size_t _end_row, size_t _beg_col, size_t _end_col):M(_M),beg_row(_beg_row), end_row(_end_row+1), beg_col(_beg_col), end_col(_end_col+1)
  {
    linbox_check((beg_row<end_row)&&(end_row<=M->rowdim()));
    linbox_check((beg_col<end_col)&&(end_col<=M->coldim()));
  }

  template<class Field>
  DenseSubMatrix<Field>::DenseSubMatrix<Field>(const DenseSubMatrix<Field>& SM,size_t _beg_row, size_t _end_row,size_t _beg_col, size_t _end_col) : M(SM.M),beg_row(SM.beg_row+_beg_row),end_row(SM.beg_row+_end_row+1),beg_col(SM.beg_col+_beg_col),end_col(SM.beg_col+_end_col+1)
  {
    linbox_check((beg_row<end_row)&&(beg_col<end_col));
    linbox_check((_end_row<SM.rowdim())&&(_end_col<SM.coldim()));
  }
  
  
  template<class Field>
  DenseSubMatrix<Field>::DenseSubMatrix<Field>(const DenseSubMatrix<Field>& SM) : M(SM.M),beg_row(SM.beg_row),end_row(SM.end_row),beg_col(SM.beg_col),end_col(SM.end_col){}

  template<class Field>
  DenseSubMatrix<Field>& DenseSubMatrix<Field>::operator=(const DenseSubMatrix<Field>& SM)
  {
    M=SM.M;
    beg_row=SM.beg_row;
    end_row=SM.end_row;
    beg_col=SM.beg_col;
    end_col=SM.end_col;
    return *this;
  }
  
  template<class Field>
  class DenseSubMatrix<Field>::RawIterator
  {
  public:
    RawIterator(){}
    RawIterator(const DenseMatrix<Field>::RawIterator& beg,
		const DenseMatrix<Field>::RawIterator& cur,
		size_t cont_len, size_t gap_len) : _beg(beg),_cur(cur),_cont_len(cont_len),_gap_len(gap_len){}
    
    RawIterator& operator=(const RawIterator& r)
    {
      _cur=r._cur;
      _beg=r._beg;
      _cont_len=r._cont_len;
      _gap_len=r._gap_len;
      return *this;
    }
    
    RawIterator& operator++()
    {
      if(((_cur-_beg+1)%_cont_len)!=0)
	++_cur;
      else
	{
	  _cur=_cur+_gap_len+1;
	  _beg=_beg+_gap_len+_cont_len;
	}
      return *this;  
    }
    
    RawIterator& operator++(int)
    {
      RawIterator tmp=*this;
      this->operator++();
      return tmp;
    }
    
    bool operator!=(const RawIterator& r) const
    {
      return (_cur!=r._cur)||(_beg!=r._beg)||(_cont_len!=r._cont_len)||(_gap_len!=r._gap_len);
    } 
    
    Element& operator*()
    {
      return *_cur;
    }
    
    const Element& operator*() const
    {
      return *_cur;
    }
      
  protected:
    DenseMatrix<Field>::RawIterator _cur;
    DenseMatrix<Field>::RawIterator _beg;
    size_t _cont_len;
    size_t _gap_len;
  };
  
  template<class Field>
  class DenseSubMatrix<Field>::ConstRawIterator
  {   
  public:
    ConstRawIterator(){}
    
    ConstRawIterator(const DenseMatrix<Field>::ConstRawIterator& beg,
		     const DenseMatrix<Field>::ConstRawIterator& cur,
		     size_t cont_len, size_t gap_len) : _beg(beg),_cur(cur),_cont_len(cont_len),_gap_len(gap_len){}
    
    ConstRawIterator& operator=(const RawIterator& r)
      {
	_cur=r._cur;
	_beg=r._beg;
	_cont_len=r._cont_len;
	_gap_len=r._gap_len;
	return *this;
      }
    
    ConstRawIterator& operator=(const ConstRawIterator& r)
    {
      _cur=r._cur;
	_beg=r._beg;
	_cont_len=r._cont_len;
	_gap_len=r._gap_len;
	return *this;
    }
    
    ConstRawIterator& operator++()
    {
      if(((_cur-_beg+1)%_cont_len)!=0)
	++_cur;
      else
	{
	  _cur=_cur+_gap_len+1;
	  _beg=_beg+_gap_len+_cont_len;
	}
      return *this;
    }
    
    ConstRawIterator& operator++(int)
    {
      RawIterator tmp=*this;
      this->operator++();
      return tmp;
    }
    
    bool operator!=(const ConstRawIterator& r) const
    {
      return (_cur!=r._cur)||(_beg!=r._beg)||(_cont_len!=r._cont_len)||(_gap_len!=r._gap_len);
    }
    
    Element& operator*()
    {
      return *_cur;
    }
    
    const Element& operator*() const
    {
      return *_cur;
    }
  protected:
    DenseMatrix<Field>::ConstRawIterator _cur;
    DenseMatrix<Field>::ConstRawIterator _beg;
    size_t _cont_len;
    size_t _gap_len;
  };
  
  template<class Field>
  DenseSubMatrix<Field>::RawIterator DenseSubMatrix<Field>::rawBegin()
  {
    return RawIterator(M->rawBegin()+beg_row*M->coldim()+beg_col, 
		       M->rawBegin()+beg_row*M->coldim()+beg_col,
		       coldim(), M->coldim()-coldim());
  }
  
  template<class Field>
  DenseSubMatrix<Field>::RawIterator DenseSubMatrix<Field>::rawEnd()
  {
    return RawIterator(M->rawBegin()+end_row*M->coldim()+beg_col,
		       M->rawBegin()+end_row*M->coldim()+beg_col, 
		       coldim(), M->coldim()-coldim());
  }
    
  template<class Field>
  DenseSubMatrix<Field>::ConstRawIterator DenseSubMatrix<Field>::rawBegin() const
  {
    return ConstRawIterator(M->rawBegin()+beg_row*M->coldim()+beg_col, 
			    M->rawBegin()+beg_row*M->coldim()+beg_col,
			    coldim(), M->coldim()-coldim());
  }
  
  template<class Field>
  DenseSubMatrix<Field>::ConstRawIterator DenseSubMatrix<Field>::rawEnd() const
  {
    return ConstRawIterator(M->rawBegin()+end_row*M->coldim()+beg_col,
			    M->rawBegin()+end_row*M->coldim()+beg_col, 
			    coldim(), M->coldim()-coldim());
  }

  template<class Field>
  size_t DenseSubMatrix<Field>::rowdim(void) const
  {
    return end_row-beg_row;
  }
  
  template<class Field>
  size_t DenseSubMatrix<Field>::coldim(void) const
  {
    return end_col-beg_col;
  }
  
  template<class Field>
  DenseSubMatrix<Field>::ColOfRowsIterator DenseSubMatrix<Field>::colOfRowsBegin()
  {
    return ColOfRowsIterator(M->rawBegin()+beg_row*M->coldim()+beg_col,end_col-beg_col,M->coldim());
  }
 
  template<class Field>
  DenseSubMatrix<Field>::ColOfRowsIterator DenseSubMatrix<Field>::colOfRowsEnd()
  {
    return ColOfRowsIterator(M->rawBegin()+end_row*M->coldim()+beg_col,end_col-beg_col,M->coldim());
  }

  template<class Field>
  DenseSubMatrix<Field>::ConstColOfRowsIterator DenseSubMatrix<Field>::colOfRowsBegin() const
  {
    return ConstColOfRowsIterator(M->rawBegin()+beg_row*M->coldim()+beg_col,end_col-beg_col,M->coldim());
  }
  
  template<class Field>
  DenseSubMatrix<Field>::ConstColOfRowsIterator DenseSubMatrix<Field>::colOfRowsEnd() const
  {
    return ConstColOfRowsIterator(M->rawBegin()+end_row*M->coldim()+beg_col,end_col-beg_col,M->coldim());
  }

  template<class Field>
  DenseSubMatrix<Field>::RowOfColsIterator DenseSubMatrix<Field>::rowOfColsBegin()
  {
    return RowOfColsIterator(M->rawBegin()+beg_col+beg_row*M->coldim(),M->coldim(),rowdim());
  }
  
  template<class Field>
  DenseSubMatrix<Field>::RowOfColsIterator DenseSubMatrix<Field>::rowOfColsEnd()
  {
    return RowOfColsIterator(M->rawBegin()+end_col+beg_row*M->coldim(),M->coldim(),rowdim());
  }
  
  template<class Field>
  DenseSubMatrix<Field>::ConstRowOfColsIterator DenseSubMatrix<Field>::rowOfColsBegin() const
  {
    return ConstRowOfColsIterator(M->rawBegin()+beg_col+beg_row*M->coldim(),M->coldim(),rowdim());
  }    
 
  template<class Field>
  DenseSubMatrix<Field>::ConstRowOfColsIterator DenseSubMatrix<Field>::rowOfColsEnd() const
  {
    return ConstRowOfColsIterator(M->rawBegin()+end_col+beg_row*M->coldim(),M->coldim(),rowdim());
  }      
  
  template<class Field>
  void DenseSubMatrix<Field>::setEntry(size_t i, size_t j, const Element& a_ij)
  {
    M->setEntry(beg_row+i, beg_col+j, a_ij);
  }
  
  template<class Field>
  DenseSubMatrix<Field>::Element& DenseSubMatrix<Field>::getEntry(size_t i, size_t j, Element& a_ij)
  {
    return M->getEntry(i+beg_row, j+beg_col,a_ij);
  } 
  
  template<class Field>
  std::ostream& DenseSubMatrix<Field>::write(std::ostream & os) const
  {
      ConstColOfRowsIterator p;
    for (p=colOfRowsBegin();p!=colOfRowsEnd();++p) 
      {
        ConstRow::const_iterator pe;
        for(pe=p->begin();pe!=p->end();++pe)
          {
            M->field().write (os, *pe);
            os << " ";
          }
        os<<"\n";
      }           
    os << endl;
    return os;
  }

  template<class Field>
  template<class Vect1, class Vect2>
  Vect1& DenseSubMatrix<Field>::apply (Vect1& y, const Vect2& x) const
  {
    ConstColOfRowsIterator p;
    typename Vect1::iterator  p_y=y.begin();  
    typename Vect2::const_iterator p_x;

    for (p=colOfRowsBegin();p!=colOfRowsEnd();++p,++p_y)
      {
        M->field().init(*p_y,0);
        for(ConstRow::iterator pe=p->begin(),p_x=x.begin();pe!=p->end();++pe,++p_x)
          M->field().axpyin(*p_y,*pe,*p_x);
      }
    
    return y;
  }
  
  template<class Field>
  template<class Iterator1, class Iterator2>
  Iterator1& DenseSubMatrix<Field>::apply (Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const
  {
    linbox_check(coldim()==(outend-outbegin));
    ConstColOfRowsIterator rowp;
    Iterator2 p_out;
    ConstRow::const_iterator pe;
    for (rowp=colOfRowsBegin();rowp!=colOfRowsEnd();++rowp,++in)
      {
        M->field().init(*in,0);
        for(pe=rowp->begin(),p_out=outbegin;pe!=rowp->end();++pe,++p_out)
          M->field().axpyin(*in,*pe,*p_out);
      }
    
    return in;
  }

  template<class Field>
  template<class Vect1, class Vect2>
  Vect1& DenseSubMatrix<Field>::applyTranspose (Vect1& y, const Vect2& x) const
  {
    ConstRowOfColsIterator colp;
    typename Vect1::iterator  p_y=y.begin();  
    typename Vect2::const_iterator p_x;

    for (colp=rowOfColsBegin();colp!=rowOfColsEnd();++colp,++p_y)
      {
        M->field().init(*p_y,0);
        for(ConstCol::const_iterator pe=colp->begin(),p_x=x.begin();pe!=colp->end();++pe,++p_x)
          M->field().axpyin(*p_y,*pe,*p_x);
      }
    
    return y;
  }

  template<class Field>
  template<class Iterator1, class Iterator2>
  Iterator1& DenseSubMatrix<Field>::applyTranspose(Iterator1 in, const Iterator2& outbegin, const Iterator2& outend) const
  {
    linbox_check(rowdim()==(outend-outbegin));
    ConstRowOfColsIterator colp;
    Iterator2 p_out;
    ConstCol::const_iterator pe;
    for (colp=rowOfColsBegin();colp!=rowOfColsEnd();++colp,++in)
      {
        M->field().init(*in,0);
        for(pe=colp->begin(),p_out=outbegin;pe!=colp->end();++pe,++p_out)
          M->field().axpyin(*in,*pe,*p_out);
      }
    
    return in;
  }

}

#endif
