// ======================================================================= //
// Invariant factors : - contains size (ni x nj) and rank
//                     - list of values and multiplicities
//
// Time-stamp: <01 Jul 99 15:52:42 Jean-Guillaume.Dumas@imag.fr> 
// (C) 1999 The Linbox group
// ======================================================================= //

#ifndef _LIN_INVARIANTS_H_
#define _LIN_INVARIANTS_H_
#include <vector.h>

template<class T>
class Invariants {
private:
    typedef long Internal;
    vector<T> _values;
    vector<Internal> _exponents;
    Internal _ni, _nj, _rank;
    
    friend ostream& operator<< <>(ostream& o, const Invariants<T>& v);

public:
    Invariants() : _rank(0),_ni(0), _nj(0) {}
    Invariants(Internal ni, Internal nj) : _rank(0), _ni(ni), _nj(nj) {}
    void PushInvariant(const T& val, Internal exp) {
        _values.push_back(val);
        _exponents.push_back(exp);
        _rank += exp;
    }
    Internal rank() const { return _rank; }
    Internal n_row() const {return _ni; }
    Internal n_col() const {return _nj; }
    Internal size() const {return _values.size(); }
    const T& getvalue( const Internal i) const { return _values[i]; }
    Internal getexponent( const Internal i) const { return _exponents[i]; }
            
};
    
    
template<class T> 
ostream& operator<< (ostream& o, const Invariants<T>& v){
        if (v._values.size())
            for(Indice i=0;i<v._values.size();i++)
                o << v._values[i] << ":" << v._exponents[i] << " ";
        return o;
};

#endif
