// ======================================================================= // (C) Linbox 2000
// Sparse Vector      : vector< Pair<T> > and an additional actual size
// Time-stamp: <11 Sep 08 13:54:23 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= 
#ifndef _SPARSE_VECTOR_H_
#define _SPARSE_VECTOR_H_
#include <iostream>

#include <linbox/vector/vector-traits.h>

// ---------------------------------------------------
//
/// Default container for sparse vectors is LightContainer
#ifndef _IBB_VECTOR_
// #include <vector>
// #define _IBB_VECTOR_ std::vector
#include "linbox/vector/light_container.h"
#define _IBB_VECTOR_ LightContainer
#endif // _IBB_VECTOR_
#ifndef _IBB_PAIR_
#include <utility>
#define _IBB_PAIR_ std::pair
// #include "linbox/vector/pair.h"
// #define _IBB_PAIR_ Pair
#endif // _IBB_PAIR_



namespace LinBox{
// ---------------------------------------------------
//
/** \brief vector< Pair<T,I> > and actualsize
\ingroup vector
*/
template<class T, class I = unsigned int>
class Sparse_Vector : public _IBB_VECTOR_< _IBB_PAIR_<I, T> > {
public:
    typedef _IBB_PAIR_<I, T>             Element;
    typedef T                      Type_t;
    typedef Sparse_Vector<T, I>    Self_t;

    // Dan Roche 6-30-04
    typedef VectorCategories::SparseSequenceVectorTag VectorCategory;


    Sparse_Vector() {};
    Sparse_Vector(size_t n) : _IBB_VECTOR_< _IBB_PAIR_<I, T> >(n), _rsize(0) {};
    Sparse_Vector(size_t n, size_t rn) : _IBB_VECTOR_< _IBB_PAIR_<I, T> >(n), _rsize(rn) {};
    ~Sparse_Vector() {};
    
            

        /// Actual dimension of the vector, 0 for infinite or unknown
    inline size_t actualsize() const { return _rsize; };            
    inline size_t reactualsize( const size_t s ) { return _rsize = s; };
    template<class XX> inline size_t reactualsize( const XX s ) { return _rsize = (size_t)s; };

    friend inline std::ostream& operator<< (std::ostream& o, const Sparse_Vector<T, I> v) {
        if (v.size())
            for(long i=0;i<v.size();i++)
                o << v[i] << std::endl;
        return o;
    }

private:
    size_t _rsize;
};    

} //end of namespace LinBox
#endif // _SPARSE_VECTOR_H_
