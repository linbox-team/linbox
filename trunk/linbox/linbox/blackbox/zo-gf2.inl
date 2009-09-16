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
	if (! _F.isZero(v) ) {
		Row_t& rowi = this->operator[](i);
		rowi.insert(
			std::lower_bound(rowi.begin(), rowi.end(), j),
			j );
	}	
}

inline ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(Element& r, size_t i, size_t j) const {
	const Row_t& rowi = this->operator[](i);
	Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
	if (there != rowi.end() )
		return r=*there;
	else
		return _F.init(r,0);
}

inline const ZeroOne<GF2>::Element& ZeroOne<GF2>::getEntry(size_t i, size_t j) const {
	static Element zero;
	const Row_t& rowi = this->operator[](i);
	Row_t::const_iterator there = std::lower_bound(rowi.begin(), rowi.end(), j);
	if (there != rowi.end() )
		return reinterpret_cast<const ZeroOne<GF2>::Element&>(*there);
	else
		return zero;
}

inline std::istream &ZeroOne<GF2>::read (std::istream &is) {
	// Reads a long int and take it mod 2 afterwards (v&1)
	UnparametricField<long> Ints;
	MatrixStream<UnparametricField<long> > S(Ints, is);
	S.getDimensions( _rowdim, _coldim );
	this->resize(_rowdim);
	Index r, c; 
	long v;
	while( S.nextTriple(r, c, v) ) {
		if (v&1) this->operator[](r).push_back(c);
	}
	for(Father_t::iterator i=this->begin(); i!=this->end(); ++i)
		std::sort(i->begin(),i->end());
	return is;
}

inline std::ostream& ZeroOne<GF2>::write (std::ostream& out, FileFormatTag format) const {
	if (format != FORMAT_GUILLAUME) 
	out << "Format other than SMS not implemented" << std::endl;
	out << _rowdim << ' ' << _coldim << " M\n";
	for(size_t i=0; i<_rowdim; ++i) {
		const Row_t& rowi = this->operator[](i);
		for(Row_t::const_iterator it=rowi.begin(); it != rowi.end(); ++it)
			out << (i+1) << ' ' << (*it+1) << " 1\n";
	}
	return out << "0 0 0" << std::endl;
}

}; // end of namespace LinBox


// Specialization of getentry
#include "linbox/solutions/getentry.h"
namespace LinBox
{
template<> struct GetEntryCategory<ZeroOne<GF2> > { typedef GetEntryTags::Local Tag; };
}; // end of namespace LinBox
