// ================================================================
// LinBox Project 1999
// Base ForwardIterator wrapper for BlackBoxes
// Have to be provided :
// - init   : sets the first value
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// Time-stamp: <06 Mar 00 18:28:47 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================
#ifndef __Base_BB_ITERATOR_H__
#define __Base_BB_ITERATOR_H__

template<class BlackBoxDomain>
class Base_BB_Iterator {
    Domain_t _domain;
    BlackBoxDomain_t * _BB_domain;
public:
        //-- Constructors
    Base_Iterator() {} 

    Base_Iterator(BlackBoxDomain_t * D) 
            : _domain(D->getdomain()), _BB_domain(D) { }
    
    BlackBoxDomain_t(BlackBoxDomain_t * BD, const Domain_t& D) 
            : _domain(D), _BB_domain(BD) {  }


    virtual Type_t& init() = 0;

        // Principal method : next
    void operator++() { _launch(); }
            
    const Type_t& operator*() { _wait(); return _value; }
    
    Domain_t getdomain() const { return _domain; }
    BlackBoxDomain_t * getBBdomain() const { return _BB_domain; }

protected:
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


#endif // __Base_BB_ITERATOR_H__
