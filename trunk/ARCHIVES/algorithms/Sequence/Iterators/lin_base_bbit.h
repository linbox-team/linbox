// ================================================================
// LinBox Project 1999
// Base ForwardIterator wrapper for BlackBoxes
// Have to be provided :
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// Time-stamp: <15 May 00 15:06:13 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================
#ifndef __Base_BB_ITERATOR_H__
#define __Base_BB_ITERATOR_H__

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif


template<class BlackBoxDomain, class Vecteur = typename BlackBoxDomain::PreferredInMatrix_t>
class Base_BB_Container {
public:
    typedef typename BlackBoxDomain::Domain_t             Domain_t;
    typedef typename BlackBoxDomain::Type_t               Type_t;
    typedef          BlackBoxDomain                       BlackBoxDomain_t;
    typedef          Base_BB_Container< BlackBoxDomain >  Self_t;

        //-- Constructors
    Base_BB_Container() {} 

    Base_BB_Container(BlackBoxDomain_t * BD) 
            : _domain(BD->getdomain()), _BB_domain(BD), _size(GIVMIN(BD->n_row(),BD->n_col()) << 1) {}
    
    Base_BB_Container(BlackBoxDomain_t * BD, const Domain_t& D) 
            : _domain(D), _BB_domain(BD), _size(GIVMIN(BD->n_row(),BD->n_col()) << 1) {}
    
    class const_iterator {
        Self_t& _c;
    public:
        const_iterator() {}
        const_iterator( Self_t& C) : _c(C) {}
        const_iterator& operator++() { _c._launch(); return *this; }
        const Type_t& operator*() { _c._wait(); return _c.getvalue(); }
    };

    const_iterator begin() { return const_iterator(*this); }
    const_iterator end() { return const_iterator(); }

    long size() { return _size; }
    const Domain_t& getdomain() const { return _domain; }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

protected:

    friend class const_iterator;
    
    virtual void _wait() = 0;
    virtual void _launch() = 0;

        /// Generic dot product using the container domain
    template <class A1, class A2>
    Type_t& DOTPROD(Type_t& coeff, const A1& u, const A2& v) {
        _domain.mul(coeff, u[0], v[0]);
        for (long k=u.size()-1 ; k>0; --k)
            _domain.axpyin(coeff, u[k], v[k]);
        return coeff;
    }

        /// User Left and Right vectors 
    Type_t& init(const Vecteur& uu, const Vecteur& vv) {
        even = 1;
        u = uu;
        v = vv;
        return DOTPROD(_value, u, u);
    }

        /// Random Left vectors, Zero Right vector
    template<class RandIter>
    Type_t& init(RandIter& g) {
        even = 1;
        u.resize(_BB_domain->n_col());
        for(long i=u.size();i--;)
            _domain.random(g, u[i]);
        v.resize(_BB_domain->n_row());
        return DOTPROD(_value, u, u);
    }
    

    const Type_t& getvalue() { return _value; }

    Domain_t _domain;
    BlackBoxDomain_t * _BB_domain;
    
    long _size;

    bool even;
    Vecteur u, v;
    Type_t _value;
};
            

#endif // __Base_BB_ITERATOR_H__
