// =================================================================== //
// Givaro : Prime numbers
//              Factor sets :
//              Pollard's rho method for factorization
//              Elliptic curves factorization by Lenstra
// Needs Container structures : stl ones for instance
// Time-stamp: <31 Aug 00 20:47:00 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef _GIVARO_FACTORISATION_H_
#define _GIVARO_FACTORISATION_H_

#include "givtimer.h"
#include "givinteger.h"
#include "givintprime.h"
#include "lin_random.h"

// #define BOUNDARY_factor TABMAX2

#define factor_first_primes(tmp,n) (tmp = iszero(mod(tmp,n,23))?23:( iszero(mod(tmp,n,19))?19:( iszero(mod(tmp,n,17))?17:  (iszero(mod(tmp,n,2))?2:( iszero(mod(tmp,n,3))?3:( iszero(mod(tmp,n,5))?5:( iszero(mod(tmp,n,7))?7: ( iszero(mod(tmp,n,11))?11:13 ))))))))

#define factor_second_primes(tmp,n) (tmp = iszero(mod(tmp,n,31))?31:( iszero(mod(tmp,n,29))?29: ( iszero(mod(tmp,n,37))?37: ( iszero(mod(tmp,n,41))?41:( iszero(mod(tmp,n,43))?43:  ( iszero(mod(tmp,n,71))?71:( iszero(mod(tmp,n,67))?67:( iszero(mod(tmp,n,61))?61:( iszero(mod(tmp,n,59))?59: ( iszero(mod(tmp,n,53))?53:( iszero(mod(tmp,n,47))?47: ( iszero(mod(tmp,n,97))?97: ( iszero(mod(tmp,n,89))?89:( iszero(mod(tmp,n,83))?83:( iszero(mod(tmp,n,79))?79:73)))))))))))))))


// =================================================================== //
// Set or Container of divisors, factors.
// =================================================================== //

template<class RandIter = Random>
class IntFactorDom : public IntPrimeDom {
private:
    // 2*3*5*7*11*13*17*19*23
    const int PROD_first_primes;
    // 29*31*37*41*43*47*53*59*61*67*71*73*79*83*89*97
    const Rep PROD_second_primes;
protected:
    RandIter _g;

public:
    typedef RandIter random_generator;

    IntFactorDom(RandIter& g = *(new RandIter())) :  IntPrimeDom(),PROD_first_primes(223092870), PROD_second_primes("10334565887047481278774629361"), _g(g) {}

    Rep& factor(Rep& r, const Rep& n) const {
        if (isone(gcd(r,n,PROD_first_primes)))
            if (isone(gcd(r,n,PROD_second_primes))) {
#ifdef GIVARO_LENSTRA
                return Lenstra((RandIter&)_g, r, n);
#else
                return Pollard((RandIter&)_g, r, n);
#endif
            } else
                return factor_second_primes(r,n);
        else 
            return factor_first_primes(r,n);
    }
       
        ///
    template< template<class> class Container> unsigned long set
        ( Container<Rep>& setint, Container<unsigned long>& setpwd,  const Rep& a) const ;
        ///
    template< template<class> class Container> void set
        ( Container<Rep>& setint, Container<long>& setpwd,  const Rep& a) const ;
        ///
    template< template<class> class Container> void set( Container<Rep>&,  const Rep&) const ;
        ///
    template< template<class> class Container> void Erathostene(Container<Rep>&, const Rep&) const ;
        /// returns a small factor
    Rep& Erathostene(Rep&,  const Rep& p ) const ;

        // Pollard with a bound on the number of loops
        // Bound 0 is without bound
    Rep& Pollard(RandIter&, Rep&, const Rep& n, unsigned long threshold = 0) const ;
        // returns a factor by Lenstra's elliptic curves method
    Rep& Lenstra(RandIter&, Rep&, const Rep& n, const Rep& B1 = 1000000, const unsigned long curves = 30) const ;
        ///
    template< template<class> class Container> Container<Rep>& divisors(Container<Rep>& L, const Container<Rep>& Lf, const Container<unsigned long>& Le)  const;
    template< template<class> class Container> Container<Rep>& divisors(Container<Rep>&, const Rep& ) const ;
    
    ostream& write(ostream& o, const Rep& n) const;


private:
// Those are parameters for Pollard's algorithms
// Pollard_fctin : must be somewhat a "random" function in Z/nZ
// Pollard_cst can be a slight alternative for the Pfct x^2+1
#define Pollard_cst 1
    
    Rep& Pollard_fctin(Rep & x, const Rep& n) const {
        mulin(x,x);
        addin(x,Pollard_cst);
        return modin(x,n);
    }

};

#include "givintfactor.inl"

#endif _GIVARO_FACTORISATION_H_
