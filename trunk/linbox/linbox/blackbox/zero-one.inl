/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/nag-sparse.h
 * Copyright (C) 2002 Rich Seagraves
 *
 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * Modified by Zhendong Wan
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __ZERO_ONE_INL
#define __ZERO_ONE_INL

#include "linbox/blackbox/zero-one.h"

namespace LinBox
{
  template<class Field>
  class ZeroOneBase<Field>::RawIterator
  {	
  public:
    typedef Element value_type;
    
    RawIterator(size_t pos, Element elem) :
      _pos(pos), _elem(elem) {}
    
    RawIterator(const RawIterator &In) :
      _pos(In._pos), _elem(In._elem) {}
    
    const RawIterator& operator=(const RawIterator& rhs) 
    {
      _pos = rhs._pos;
      _elem = rhs._elem;
      return *this;
    }
    
        
    bool operator==(const RawIterator &rhs) 
    {
      return ( _pos == rhs._pos && _elem == rhs._elem);
    }
    
    bool operator!=(const RawIterator &rhs) 
    {
      return ( _pos != rhs._pos || _elem != rhs._elem );
    }
    
    RawIterator & operator++() 
    {
      ++_pos;
      return *this;
    }
    
    RawIterator operator++(int) 
    {
      
      RawIterator tmp = *this;
      _pos++;
      return tmp;
    }
    
    value_type operator*() 
    {
      return _elem;
    }
    
    const value_type operator*() const 
    {
      return _elem;
    }
    
  private:
    value_type _elem;
    size_t _pos;
  };  
  
  /* STL standard Begin and End functions.  Used to get
   * the beginning and end of the data.  So that RawIterator
   * can be used in algorithms like a normal STL iterator.
   */
  template<class Field> typename
  ZeroOneBase<Field>::RawIterator ZeroOneBase<Field>::rawBegin()
  { return RawIterator( 0, _F.init(_tmp, 1) ); }
  
  template<class Field> typename
  ZeroOneBase<Field>::RawIterator ZeroOneBase<Field>::rawEnd() 
  { return RawIterator( _nnz, _F.init(_tmp, 1) ); }
  
  template<class Field> 
  const typename ZeroOneBase<Field>::RawIterator ZeroOneBase<Field>::rawBegin() const
  { return RawIterator(0, _F.init(_tmp, 1) ); }

  template<class Field> 
  const typename ZeroOneBase<Field>::RawIterator ZeroOneBase<Field>::rawEnd() const 
  { return RawIterator(_nnz, _F.init(_tmp, 1) ); } 
  
  /* RawIndexIterator - Iterates through the i and j of the current element
   * and when accessed returns an STL pair containing the coordinates
   */
  template<class Field>
  class ZeroOneBase<Field>::RawIndexIterator 
  {
  public:
    typedef std::pair<size_t, size_t> value_type;
    
    RawIndexIterator() {}
    
    RawIndexIterator(size_t* row, size_t* col):
      _row(row), _col(col) {}
    
    RawIndexIterator(const RawIndexIterator &In):
      _row(In._row), _col(In._col) {}
    
    const RawIndexIterator &operator=(const RawIndexIterator &rhs) 
    {
      _row = rhs._row;
      _col = rhs._col;
      return *this;				
    }
    
    bool operator==(const RawIndexIterator &rhs) 
    {
      _row == rhs._row && _col == rhs._col;     
    }
    
    bool operator!=(const RawIndexIterator &rhs) 
    {
      return _row != rhs._row || _col != rhs._col;
    }
    
    const RawIndexIterator& operator++() 
    {
      ++_row; ++_col;	
      return *this;
    }
    
    const RawIndexIterator operator++(int) 
    {
      RawIndexIterator tmp = *this;
      ++_row; ++_col;
      return tmp;
    }
    
    value_type operator*() 
    {
      return std::pair<size_t,size_t>(*_row, *_col);
    }
    
    const value_type operator*() const 
    {
      return std::pair<size_t,size_t>(*_row, *_col);
    }
  private:
    size_t* _row, *_col;
  };
  
  template<class Field> typename
  ZeroOneBase<Field>::RawIndexIterator ZeroOneBase<Field>::indexBegin() 
  {
    return RawIndexIterator(_rowP, _colP);
  }
  
  template<class Field> 
  const typename ZeroOneBase<Field>::RawIndexIterator ZeroOneBase<Field>::indexBegin() const
  {
    return RawIndexIterator(_rowP, _colP);
  }

  template<class Field> typename
  ZeroOneBase<Field>::RawIndexIterator ZeroOneBase<Field>::indexEnd() 
  {
    return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
  }

  template<class Field> 
  const typename ZeroOneBase<Field>::RawIndexIterator ZeroOneBase<Field>::indexEnd() const 
  {
    return RawIndexIterator(_rowP + _nnz, _colP + _nnz);
  }
 
  template<class Field>
  ZeroOneBase<Field>::ZeroOneBase() { srand( time(NULL) ); dynamic = false;}
    
  
  template<class Field>
  ZeroOneBase<Field>::ZeroOneBase(Field F, Index* rowP, Index* colP, Index rows, Index cols, Index NNz, bool rowSort, bool colSort):
    _F(F), _rows(rows), _cols(cols), _nnz(NNz), _rowP(rowP), _colP(colP), _rowSort(rowSort), _colSort(colSort) , dynamic(false) { srand(time(NULL)); }
  
  template<class Field>
  ZeroOneBase<Field>::~ZeroOneBase() 
  {
	  if(dynamic) {
		  delete [] _rowP;
		  delete [] _colP;
	  }
  }

   
  template<class Field>
  size_t ZeroOneBase<Field>::rowdim() const
  {
    return _rows;
  }
  
  
  template<class Field>
  size_t ZeroOneBase<Field>::coldim() const
  {
    return _cols;
  }
    
  template<class Field>
  void ZeroOneBase<Field>::rowSort() const
  {
    int mode = 0;
    if( _rowSort) return;  // Already sorted, we're done   
    else _qsort( (size_t) 0, _nnz, mode);
    _rowSort = true;
    return;
  }

  template<class Field>
  void ZeroOneBase<Field>::colSort() const
  {
    int mode = 1;
    if( _colSort) return; // Already sorted, good to go  
    else _qsort( (size_t) 0, _nnz, mode);
    _colSort = true; _rowSort = false;
    return;
  }   
  
  template<class Field>
  void ZeroOneBase<Field>::_qsort(size_t p, size_t e, int &mode) const
  {
    int i;
    if( (e - p) <= 1) ;
    else 
      {
	i = 1 + _part(p, e, mode);
	_qsort(p, i, mode);
	_qsort(i, e, mode);
      }
  }
  
  template<class Field>
  size_t ZeroOneBase<Field>::_part(size_t p, size_t e, int &mode) const
  {
    size_t rtemp, ctemp, rowval, colval;
    int i = p + rand() % (e - p), j = e;
    rtemp = _rowP[p];
    ctemp = _colP[p];
    _rowP[p] = _rowP[i];
    _colP[p] = _colP[i];
    _rowP[i] = rtemp;
    _colP[i] = ctemp;
    rowval = _rowP[p];
    colval = _colP[p];
    i = p - 1;
    
    if(mode == 0) 
      { // Row mode, go by row order, then column
	while(true) 
	  {
	    do j--; while( _rowP[j] > rowval || ( _rowP[j] == rowval && _colP[j] > colval ));
	    do i++; while( _rowP[i] < rowval || ( _rowP[i] == rowval && _colP[i] < colval ));
	    if( i < j) 
	      {
		rtemp = _rowP[j];
		ctemp = _colP[j];
		_rowP[j] = _rowP[i];
		_colP[j] = _colP[i];
		_rowP[i] = rtemp;
		_colP[i] = ctemp;
	      }
	    else return j;
	  }
      }
    else 
      { // Col mode, go by col order, then row
	while(true) 
	  {
	    do j--; while( _colP[j] > colval || ( _colP[j] == colval && _rowP[j] > rowval ));
	    do i++; while( _colP[i] < colval || ( _colP[i] == colval && _rowP[i] < rowval ));
	    if( i < j) {
	      rtemp = _rowP[j];
	      ctemp = _colP[j];
	      _rowP[j] = _rowP[i];
	      _colP[j] = _colP[i];
	      _rowP[i] = rtemp;
	      _colP[i] = ctemp;
	    }
	    else return j;
	  }
      }
  }
      
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOneBase<Field>::applySpecialization(OutVector & y, const InVector & x, const NormField& n) const
  {
    //std::cout<<"Call general case\n";
    linbox_check((y.size()==rowdim())&&(x.size()==coldim()));         
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    Index* ip, *jp;
    
    // 0 out y.  Note, this implementation assumes a dense vector.
    for(yp = y.begin(); yp != y.end(); ++yp)
      _F.init(*yp , 0);
    
    rowSort();
    
    yp=y.begin();
    xp=x.begin();
    ip=_rowP;
    jp=_colP;
    size_t rowI =0;
    
    for(; ip <_rowP+nnz(); ++ip,++jp) 
      {       
	if( *ip == rowI)
	  _F.addin(*yp,*(xp + *jp));
	else
	  {
	    if((*ip-rowI)==1)
	      ++yp;
	    else	      
	      yp=y.begin()+*ip;
	    
	    rowI=*ip;
	    _F.addin(*yp,*(xp + *jp));
	  }
      }  
    return y;
  }
  
   
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOneBase<Field>::applySpecialization(OutVector & y, const InVector & x, const Mod32Field& m) const
  {
    //std::cout<<"Called specialization\n";
    linbox_check((y.size()==rowdim())&&(x.size()==coldim()));
    
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    Index* ip, *jp;
        
    for(yp = y.begin(); yp != y.end(); ++yp)
      _F.init(*yp , 0);
    
    rowSort();
    
    yp=y.begin();
    xp=x.begin();
    ip=_rowP;
    jp=_colP;
    size_t rowI =0;
    integer _prime;

    _F.characteristic(_prime);
    
    uint32 prime = static_cast<uint32>(_prime);
    
    uint64 accum =0;
    
    for(; ip <_rowP+nnz(); ++ip,++jp) 
      {       
	if( *ip == rowI)
	  accum=accum+*(xp + *jp);
	else
	  {
	    *yp= accum % prime;
	    if((*ip-rowI)==1)
	      ++yp;
	    else	      
	      yp=y.begin()+*ip;
	    
	    rowI=*ip;	    
	    accum=*(xp+*jp);
	  }
      }
    if(rowI)
      *yp=accum % prime;
    
    return y;
  }
 
  
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOneBase<Field>::applyTransposeSpecialization(OutVector & y, const InVector & x, const NormField& n) const
  {
    //std::cout<<"Call general case\n";
    linbox_check((y.size()==coldim())&&(x.size()==rowdim()));   
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    Index* ip, *jp;
    
    // 0 out y.  Note, this implementation assumes a dense vector.
    for(yp = y.begin(); yp != y.end(); ++yp)
      _F.init(*yp , 0);
    
    rowSort();
    
    yp=y.begin();
    xp=x.begin();
    ip=_rowP;
    jp=_colP;
    size_t rowI =0;
    
    for(; ip <_rowP+nnz(); ++ip,++jp) 
      {       
	if( *ip == rowI)
	  _F.addin(*(yp+*jp),*xp);
	else
	  {
	    if((*ip-rowI)==1)
	      ++xp;
	    else	      
	      xp=x.begin()+*ip;
	    
	    rowI=*ip;
	    _F.addin(*(yp+*jp),*xp);
	  }
      }	
  
    return y;
  }
  
  
  template<class Field>
  template<class OutVector, class InVector>
  OutVector & ZeroOneBase<Field>::applyTransposeSpecialization(OutVector & y, const InVector & x, const Mod32Field& m) const
  {
    //std::cout<<"Called specialization\n";
    linbox_check((y.size()==coldim())&&(x.size()==rowdim()));
    
    std::vector<uint64> y_c (y.size(),0);
    
    typename OutVector::iterator yp;
    typename InVector::const_iterator xp;
    Index* ip, *jp;      
    
    rowSort();
     
    xp=x.begin();
    ip=_rowP;
    jp=_colP;
    size_t rowI =0;
    std::vector<uint64>::iterator y_cp;
    y_cp=y_c.begin();
    
    for(; ip <_rowP+nnz(); ++ip,++jp) 
      {       
	if( *ip == rowI)
	  *(y_cp+*jp) += *xp;
	else
	  {
	    if((*ip-rowI)==1)
	      ++xp;
	    else	      
	      xp=x.begin()+*ip;
	    
	    rowI=*ip;	    
	    *(y_cp+*jp) += *xp;
	  }
      }
    
    integer _prime;
    _F.characteristic(_prime);    
    uint32 prime = static_cast<uint32>(_prime);
    
    yp=y.begin();
    y_cp=y_c.begin();
    for(;yp!=y.end();++yp,++y_cp)
      *yp = (*y_cp) % prime;
    
    return y;
  }
       
  template<class Field>
  size_t ZeroOneBase<Field>::nnz() const
  {
    return _nnz;
  }

#ifdef __LINBOX_XMLENABLED

 template<class Field, class Vector>
 ostream &ZeroOne<Field, Vector>::write(ostream &out) const
 {
	 Writer W;
	 if( toTag(W)) 
		 W.write(out);

	 return out;
 }

 template<class Field, class Vector>
 bool ZeroOne<Field, Vector>::toTag(Writer &W) const
 {
	 size_t i;
	 vector<size_t> rows, cols;
	 string s;
	 W.setTagName("MatrixOver");
	 W.setAttribute("rows", Writer::numToString(s, _rows));
	 W.setAttribute("cols", Writer::numToString(s, _cols));
	 W.setAttribute("implDetail", "zero-one");

	 W.addTagChild();
	 _F.toTag(W);
	 W.upToParent();

	 W.addTagChild();
	 W.setTagName("zero-one");
	 
	 for(i = 0; i < _nnz; ++i) {
		 rows.push_back(_rowP[i]);
		 cols.push_back(_colP[i]);
	 }

	 W.addTagChild();
	 W.setTagName("index");
	 W.addNumericalList(rows);
	 W.upToParent();

	 W.addTagChild();
	 W.setTagName("index");
	 W.addNumericalList(cols);
	 W.upToParent();

	 W.upToParent();

	 return true;
 }



 template<class Field>
 ZeroOneBase<Field>::ZeroOneBase(Reader &R) : _F(R.Down(1))
 {
	 vector<size_t> rows, cols;
	 size_t i;

	 R.Up(1);
	 if(!R.expectTagName("MatrixOver") ) return;
	 if(!R.expectAttributeNum("rows", _rows) || !R.expectAttributeNum("cols", _cols)) return;


	 if(!R.expectChildTag()) return;

	 R.traverseChild();
	 if(!R.expectTagName("field")) return;
	 R.upToParent();

	 if(!R.getNextChild()) {
		 R.setErrorString("Got a matrix with a field and no data.");
		 R.setErrorCode(Reader::OTHER);
		 return;
	 }

	 if(!R.expectChildTag()) return;
	 
	 R.traverseChild();
	 if(!R.expectTagName("zero-one") || !R.expectChildTag()) return;

	 R.traverseChild();
	 if(!R.expectTagName("index") || !R.expectTagNumVector(rows)) return;
	 R.upToParent();

	 if(!R.getNextChild()) {
		 R.setErrorString("Didn't get columnar indices for zero-one matrix");
		 R.setErrorCode(Reader::OTHER);
		 return;
	 }

	 if(!R.expectChildTag()) return;
	 R.traverseChild();
	 if(!R.expectTagName("index") || !R.expectTagNumVector(cols)) return;
	 R.upToParent();
	 R.upToParent();
	 R.getPrevChild();

	 dynamic = true;
	 _rowP = new size_t[rows.size()];
	 _colP = new size_t[rows.size()];
	 _nnz = rows.size();
	 _rowSort = _colSort = false;

	 for(i = 0; i < _nnz; ++i) {
		 _rowP[i] = rows[i];
		 _colP[i] = cols[i];
	 }
	 
	 return;

 }

 template<class Field>
 ZeroOneBase<Field>::ZeroOneBase(const ZeroOneBase<Field> &M) : _F(M._F)
 {
	 size_t i;

	 dynamic = true;
	 _nnz = _M._nnz;
	 _rowP = new size_t[_nnz];
	 _colP = new size_t[_nnz];
	 _rows = M._rows;
	 _cols = M._cols;

	 for(i = 0; i < _nnz; ++i) {
		 _rowP[i] = M._rowP[i];
		 _colP[i] = M._colP[i];
	 }
 }

 template<class Field, class Vector>
 ZeroOne<Field, Vector>::ZeroOne(Reader &R) : ZeroOneBase<Field>(R) {}

 template<class Field, class Vector>
 ZeroOne<Field, Vector>::ZeroOne(const ZeroOne<Field, Vector> &M) : ZeroOneBase<Field>(M) {}

#endif	 
	 
	 
   
}

//End of LinBox
#include "linbox/blackbox/zero-one.h"

#endif
