// ======================================================================= // (C) Linbox 2000
// Sparse Vector      : vector< Pair<T> > and an additional actual size
// Time-stamp: <19 Sep 03 14:23:40 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= 
#ifndef _SPARSE_VECTOR_H_
#define _SPARSE_VECTOR_H_
#include <iostream>

#include "linbox/vector/pair.h"
// ---------------------------------------------------
//
/// Default container for sparse vectors is STL vector
#ifndef _IBB_VECTOR_
#include <vector>
#define _IBB_VECTOR_ std::vector
#endif // _IBB_VECTOR_

// ---------------------------------------------------
//
/// Sparse Vector : vector< Pair<T> > and actualsize
template<class T, class I = unsigned long>
class Sparse_Vector : public _IBB_VECTOR_< Pair<T, I> > {
public:
    typedef Pair<T, I>             Element;
    typedef T                      Type_t;
    typedef Sparse_Vector<T, I>    Self_t;

    Sparse_Vector() {};
    Sparse_Vector(size_t n) : _IBB_VECTOR_< Pair<T, I> >(n), _rsize(0) {};
    Sparse_Vector(size_t n, size_t rn) : _IBB_VECTOR_< Pair<T, I> >(n), _rsize(rn) {};
    ~Sparse_Vector() {};
    
            

        /// Actual dimension of the vector, 0 for infinite or unknown
    inline size_t actualsize() const { return _rsize; };            
    inline size_t reactualsize( const size_t s ) { return _rsize = s; };
    template<class XX> inline size_t reactualsize( const XX s ) { return _rsize = (size_t)s; };

    friend inline std::ostream& operator<< (std::ostream& o, const Sparse_Vector<T, I> v) {
        if (v.size())
            for(long i=0;i<v.size();i++)
                o << v[i] << endl;
        return o;
    }

private:
    size_t _rsize;
};    


#endif // _SPARSE_VECTOR_H_
