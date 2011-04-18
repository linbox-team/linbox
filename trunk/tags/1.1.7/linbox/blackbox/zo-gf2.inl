/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_zo_gf2_INL
#define __LINBOX_zo_gf2_INL

#include <givaro/givintfactor.h>

namespace LinBox
{
// Dot product structure enabling std::transform call
template<class Blackbox, class InVector>
struct dotp {
    const typename Blackbox::Field& _F;
    const InVector& _x;
    dotp(const typename Blackbox::Field& F, const InVector& x) : _F(F), _x(x) {}
        
    bool operator()(const typename Blackbox::Row_t& row) const {
            bool tmp(false);
            for(typename Blackbox::Row_t::const_iterator loc = row.begin(); loc != row.end(); ++loc) {
                _F.addin(tmp,_x[*loc]);
            }
            return tmp;
    }
};

#include <algorithm>
template<class OutVector, class InVector>
inline OutVector & ZeroOne<GF2>::apply(OutVector & y, const InVector & x) const {
    dotp<Self_t,InVector> mydp(this->_F, x);
    std::transform(this->begin(), this->end(), y.begin(), mydp );
    return y;
}
    
/*
template<class OutVector, class InVector>
inline OutVector & ZeroOne<GF2>::apply(OutVector & y, const InVector & x) const {
    typename OutVector::iterator yit = y.begin();
    Self_t::const_iterator row = this->begin();
    for( ; row != this->end(); ++yit, ++row) {
        bool tmp(false);
        for(Row_t::const_iterator loc = row->begin();loc != row->end(); ++loc)
            _F.addin(tmp,x[*loc]);
        *yit = tmp;
    }
    return y;
}
*/
    
template<class OutVector, class InVector>
inline OutVector & ZeroOne<GF2>::applyTranspose(OutVector & y, const InVector & x) const {
    std::fill(y.begin(),y.end(),false);
    typename InVector::const_iterator xit = x.begin();
    Self_t::const_iterator row = this->begin();
    for( ; row != this->end(); ++row, ++xit) {
        for(typename Self_t::Row_t::const_iterator loc = row->begin(); loc != row->end(); ++loc) {
            _F.addin(y[*loc],*xit);
        }
    }
    return y;
}
    

inline void ZeroOne<GF2>::setEntry(size_t i, size_t j, const Element& v) {
    Row_t& rowi = this->operator[](i);
    Row_t::iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
    if (! _F.isZero(v) ) {
        if ( (there == rowi.end() ) || (*there != j) ) {
            rowi.insert(there, j);
            ++_nnz;
        }
    } else {
        if ( (there != rowi.end() ) && (*there == j) ) {
            rowi.erase(there);
            --_nnz;           
        }
    }
}

inline ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(Element& r, size_t i, size_t j) const {
	const Row_t& rowi = this->operator[](i);
	Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
	if (there != rowi.end() )
		return r=*there;
	else
		return r=_F.zero;
}

inline const ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(size_t i, size_t j) const {
	const Row_t& rowi = this->operator[](i);
	Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
	if (there != rowi.end() )
		return reinterpret_cast<const ZeroOne<GF2>::Element&>(*there);
	else
		return _F.zero;
}

inline std::istream &ZeroOne<GF2>::read (std::istream &is) {
	// Reads a long int and take it mod 2 afterwards (v&1)
	UnparametricField<long> Ints;
	MatrixStream<UnparametricField<long> > S(Ints, is);
	S.getDimensions( _rowdim, _coldim );
	this->resize(_rowdim);
	Index r, c; 
	long v;
        _nnz = 0;
	while( S.nextTriple(r, c, v) ) {
		if (v&1) {
                    this->operator[](r).push_back(c);
                    ++_nnz;
                }
	}
	for(Father_t::iterator i=this->begin(); i!=this->end(); ++i)
		std::sort(i->begin(),i->end());
	return is;
}

inline std::ostream& ZeroOne<GF2>::write (std::ostream& out, FileFormatTag format) const {
    if (format == FORMAT_GUILLAUME) {
	out << _rowdim << ' ' << _coldim << " M\n";
	for(size_t i=0; i<_rowdim; ++i) {
            const Row_t& rowi = this->operator[](i);
            for(Row_t::const_iterator it=rowi.begin(); it != rowi.end(); ++it)
                out << (i+1) << ' ' << (*it+1) << " 1\n";
	}
	return out << "0 0 0" << std::endl;
    } else if (format == FORMAT_MAPLE) {
        out << '[';
        bool firstrow=true;
        for (const_iterator i = begin (); i != end (); ++i) {
            if (firstrow) {
                out << '[';
                firstrow =false;
            } else 
                out << ", [";
            
            Row_t::const_iterator j = i->begin ();
            for (long j_idx = 0; j_idx < static_cast<long>(_coldim); j_idx++) {
                if (j == i->end () || j_idx != static_cast<long>(*j) )
                    out << '0';
                else {
                    out << '1';
                    ++j;
                }
                if (j_idx < (static_cast<long>(_coldim)-1) )
                    out << ',';
            }

            out << ']';
        }
        return out << ']';
    } else
        return out << "ZeroOne over GF(2), format other than SMS or Maple not implemented" << std::endl;
}

  class ZeroOne<GF2>::RawIterator
  {	
  public:
    typedef Element value_type;
    
    RawIterator(size_t pos, Element elem) :
      _elem(elem),_pos(pos)  {}
    
    RawIterator(const RawIterator &In) :
      _elem(In._elem),_pos(In._pos) {}
    
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
    
    value_type operator*() { return _elem; }
    
    const value_type operator*() const { return _elem; }
    
  private:
    value_type _elem;
    size_t _pos;
  };  
  
  /* STL standard Begin and End functions.  Used to get
   * the beginning and end of the data.  So that RawIterator
   * can be used in algorithms like a normal STL iterator.
   */
  inline ZeroOne<GF2>::RawIterator ZeroOne<GF2>::rawBegin()
  { return RawIterator( 0, _F.one ); }
  
  inline ZeroOne<GF2>::RawIterator ZeroOne<GF2>::rawEnd() 
  { return RawIterator( _nnz, _F.one ); }
  
  inline const ZeroOne<GF2>::RawIterator ZeroOne<GF2>::rawBegin() const
  { return RawIterator(0, _F.one ); }

  inline const ZeroOne<GF2>::RawIterator ZeroOne<GF2>::rawEnd() const 
  { return RawIterator(_nnz, _F.one ); } 
  
  /* RawIndexIterator - Iterates through the i and j of the current element
   * and when accessed returns an STL pair containing the coordinates
   */
  class ZeroOne<GF2>::RawIndexIterator 
  {
  public:
    typedef std::pair<size_t, size_t> value_type;
    
    RawIndexIterator() {}
    
      RawIndexIterator(size_t rowidx, 
                       LightContainer<LightContainer<size_t> >::const_iterator rowbeg, 
                       LightContainer<LightContainer<size_t> >::const_iterator rowend, 
                       size_t colidx, 
                       LightContainer<size_t>::const_iterator colbeg)
              : _rowbeg( LightContainer<LightContainer<size_t> >::iterator(rowbeg) ), 
                _rowend( LightContainer<LightContainer<size_t> >::iterator(rowend) ), 
                _colbeg( LightContainer<size_t>::iterator(colbeg) ), 
                _row(rowidx), 
                _col(colidx) {
          
          if( _rowbeg == _rowend ) return;
          
          while ( _colbeg == _rowbeg->end() ) {
              
              if (++_rowbeg == _rowend) return;

              _colbeg = _rowbeg->begin();

          }

    }
    
    RawIndexIterator(const RawIndexIterator &In):
      _rowbeg(In._rowbeg), _rowend(In._rowend), _colbeg(In._colbeg), _row(In._row), _col(In._col) {}
    
    const RawIndexIterator &operator=(const RawIndexIterator &rhs) 
    {
      _rowbeg = rhs._rowbeg;
      _rowend = rhs._rowend;
      _colbeg = rhs._colbeg;
      _row = rhs._row;
      _col = rhs._col;
      return *this;				
    }
    
    bool operator==(const RawIndexIterator &rhs) 
    {
      return _rowbeg == rhs._rowbeg && _colbeg == rhs._colbeg;     
    }
    
    bool operator!=(const RawIndexIterator &rhs) 
    {
      return _rowbeg != rhs._rowbeg || _colbeg != rhs._colbeg;
    }
    
    const RawIndexIterator& operator++() {

        

        ++_colbeg;
        while(_colbeg == _rowbeg->end()) {
            if (++_rowbeg == _rowend) return *this;
            ++_row;	
            _colbeg = _rowbeg->begin();
        }
        _col = *_colbeg;
	

        return *this;
    }
      
    const RawIndexIterator operator++(int) 
    {
      RawIndexIterator tmp = *this;
      this->operator++();
      return tmp;
    }
    
    value_type operator*() 
    {
      return std::pair<size_t,size_t>(_row, _col);
    }
    
    const value_type operator*() const 
    {
      return std::pair<size_t,size_t>(_row, _col);
    }
  private:
      LightContainer<LightContainer<size_t> >::iterator _rowbeg, _rowend;
      LightContainer<size_t>::iterator _colbeg;
      size_t _row, _col;
  };
  
  inline ZeroOne<GF2>::RawIndexIterator ZeroOne<GF2>::indexBegin() 
  {
    return RawIndexIterator(0, this->begin(), this->end(), 0, this->front().begin() );
  }
  
  inline const ZeroOne<GF2>::RawIndexIterator ZeroOne<GF2>::indexBegin() const
  {
    return RawIndexIterator(0, this->begin(), this->end(), 0, this->front().begin() );
  }

  inline ZeroOne<GF2>::RawIndexIterator ZeroOne<GF2>::indexEnd() 
  {
    return RawIndexIterator(_rowdim, this->end(), this->end(), this->back().size(),this->back().end() );
  }

  inline const ZeroOne<GF2>::RawIndexIterator ZeroOne<GF2>::indexEnd() const 
  {
    return RawIndexIterator(_rowdim, this->end(), this->end(), this->back().size(),this->back().end() );
  }
 




}; // end of namespace LinBox


// Specialization of getentry
#include "linbox/solutions/getentry.h"
namespace LinBox
{
template<> struct GetEntryCategory<ZeroOne<GF2> > { typedef GetEntryTags::Local Tag; };
} // end of namespace LinBox

#endif //__LINBOX_zo_gf2_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
