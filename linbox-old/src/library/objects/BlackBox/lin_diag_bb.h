// =========================================================
// (C) The LinBox Group 1999
// Linbox implementation of diagonal BlackBoxes
// file : lin_dom_diag_bb.h
// Time-stamp: <06 Mar 00 17:42:35 Jean-Guillaume.Dumas@imag.fr> 
// =========================================================


template<class Domain>
class DiagDom {
public:
    typedef typename Domain::element                    Type_t;
    typedef typename vector<Type_t>                     PreferredInMatrix_t;
    typedef typename vector<Type_t>                     PreferredOutMatrix_t;
    typedef          (Type_t *)                         element;
    typedef typename Domain                             Domain_t;
    typedef          DiagDom< BlackBoxDomain >          Self_t;
private:
    Type_t * _values;
    Indice _size;
    Domain_t _domain;
public:
    
    BBDiag() : _values(0), _size(0) {};

    BBDiag(Indice n, Type_t * v, const Domain_t& D) : _values(v), _size(n), _domain(D)  {};        

    long n_row() const { return _size; }
    long n_col() const { return _size; }
    long n_elem() const { return _size; }

    template<class OutMatrix, class InMatrix>
    OutMatrix& Apply(OutMatrix& res,  const InMatrix& vect ) {
        res.resize(n_row());
        Type_t * it(_values);
        Indice i(0);
        for(;i!=_size;++i,++it)
            _domain.mul(res[i], *it, vect[i]);
        return res;
    }

    template<class OutMatrix, class InMatrix>
    OutMatrix& ApplyTrans(OutMatrix& res,  const InMatrix& vect ) {
        return Apply(res,vect);
    }

};

