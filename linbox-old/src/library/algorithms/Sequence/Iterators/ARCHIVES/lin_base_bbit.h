// ================================================================
// LinBox Project 1999
// Base ForwardIterator wrapper for BlackBoxes
// Have to be provided :
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// Time-stamp: <08 Mar 00 14:56:20 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================
#ifndef __Base_BB_ITERATOR_H__
#define __Base_BB_ITERATOR_H__

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif

template<class BlackBoxDomain>
class Base_BB_Iterator {
public:
    typedef typename BlackBoxDomain::Domain_t             Domain_t;
    typedef typename BlackBoxDomain::Type_t               Type_t;
    typedef          BlackBoxDomain                       BlackBoxDomain_t;
    typedef          Base_BB_Iterator< BlackBoxDomain >         Self_t;
        //-- Constructors
    Base_BB_Iterator() {} 

    Base_BB_Iterator(BlackBoxDomain_t * D) 
            : _domain(D->getdomain()), _BB_domain(D) { }
    
    Base_BB_Iterator(BlackBoxDomain_t * BD, const Domain_t& D) 
            : _domain(D), _BB_domain(BD) {  }

    void operator++() { _launch(); }
            
    const Type_t& operator*() { _wait(); return _value; }
    
    Domain_t getdomain() const { return _domain; }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

protected:
    Domain_t _domain;
    BlackBoxDomain_t * _BB_domain;
    Type_t _value;

    virtual void _wait() = 0;
    virtual void _launch() = 0;

    template <class A1, class A2>
    Type_t& DOTPROD(Type_t& coeff, const A1& u, const A2& v) {
        _domain.mul(coeff, u[0], v[0]);
        for (long k=u.size()-1 ; k>0; --k)
            _domain.axpyin(coeff, u[k], v[k]);
        return coeff;
    }

};


template<class BlackBoxDomain>
class Base_BB_Container {
public:
    typedef typename BlackBoxDomain::Domain_t             Domain_t;
    typedef typename BlackBoxDomain::Type_t               Type_t;
    typedef          BlackBoxDomain                       BlackBoxDomain_t;
    typedef          Base_BB_Container< BlackBoxDomain >  Self_t;
    typedef          Base_BB_Iterator< BlackBoxDomain >   const_iterator;

        //-- Constructors
    Base_BB_Container() {} 

    Base_BB_Container(BlackBoxDomain_t * BD) 
            : _domain(BD->getdomain()), _BB_domain(BD), _size(GIVMIN(BD->n_row(),BD->n_col()) << 1) {}
    
    Base_BB_Container(BlackBoxDomain_t * BD, const Domain_t& D) 
            : _domain(D), _BB_domain(BD), _size(GIVMIN(BD->n_row(),BD->n_col()) << 1) {}
    
    virtual const_iterator& begin() const = 0;
    virtual const_iterator& end() const = 0;

    long size() { return _size; }
    Domain_t getdomain() const { return _domain; }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

    
protected:
    Domain_t _domain;
    BlackBoxDomain_t * _BB_domain;
    long _size;
};
            

#endif // __Base_BB_ITERATOR_H__
