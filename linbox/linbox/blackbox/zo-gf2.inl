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
    
}; // end of namespace LinBox
