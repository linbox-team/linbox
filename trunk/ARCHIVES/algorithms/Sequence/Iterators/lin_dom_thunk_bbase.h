// ================================================================
// LinBox Project 1999
// Base Black Box Thunk
// Next is ready for the diffenrent blackbox implementations
// Have to be provided :
// - init   : sets the first value
// - launch : launches the following computation
// - wait   : waits for the end of the current computation
// Time-stamp: <20 Dec 99 11:04:51 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================


/* A Thunk<T> is a nullary function object returning a T.
   Normally a thunk has state and may return different values on 
   successive calls.  Thus the following code outputs "2 4 16 256".
     Thunk<int> t(square, 2);
     cout << t() << t() << t() << t();

   An example of a thunk is random(), 
   A purely functional thunk, such as one(), represents a constant.  
   and returns the same value each time,  The term thunk is seldom 
   applied to such objects.  A more proper example is random(), whose
   semantics and internal state you can imagine.  Note that thunks
   are not in general thread safe, although it would be reasonable
   to compute a copy constructor such that both copies would compute
   the same sequence of values on future calls but independently.
   
   A common way to create a thunk is with a unary function f and 
   initial value v.  The sequence of values returned by successive 
   calls to the thunk is v, f(v), f(f(v)), ...  Such thunks may
   be lazy, diligent, or eager.  A lazy thunk stores the previously
   output value v.  When called it executes the steps 
     v = f(v); return u;
   A diligent thunk stores the next value to return.  When called it executes
     return v; v = f(v); 
   (What is meant is that it can launch computation of the next value and 
   quickly return the current value, if ready).
   
   An eager thunk might speculatively compute several future values
   in parallel with the calling routine.  I know of no application for
   this feature.

   The usefulness of thunks is the encapsulation and the opportunity
   for parallelism in the diligent implementation.
*/

#ifndef __THUNK_B_B_H__
#define __THUNK_B_B_H__

template<class BlackBoxDomain>
class BBThunkDom {
public:
    typedef typename BlackBoxDomain::Domain_t             Domain_t;
    typedef typename BlackBoxDomain::Type_t               Type_t;
    typedef          BlackBoxDomain                       BlackBoxDomain_t;
    typedef          BBThunkDom< BlackBoxDomain >      Self_t;
    typedef typename BlackBoxDomain::PreferredInMatrix_t  PreferredInMatrix_t;
    typedef typename BlackBoxDomain::PreferredInMatrix_t  PreferredOutMatrix_t;
protected:
    typedef typename BlackBoxDomain::PreferredInMatrix_t  Vecteur;

    Domain_t _domain;
    BlackBoxDomain_t * _BB_domain;
public:
        //-- Constructors
    BBThunkDom() {} 

    BBThunkDom(BlackBoxDomain_t * D) 
            : _domain(D->getdomain()), _BB_domain(D) { }
    
    BBThunkDom(BlackBoxDomain_t * BD, const Domain_t& D) 
            : _domain(D), _BB_domain(BD) {  }


    virtual Type_t& init() = 0;

        // Principal method : next
    Type_t& next(Type_t& r) { 
        _wait();
        r = _value; 
        _launch(); 
        return r;
    }

    Type_t& operator()(Type_t& u) { return next(u); }
 
    long n_row() const { return _BB_domain->n_row(); }
    long n_col() const { return _BB_domain->n_col(); }
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


#endif
