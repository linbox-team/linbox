// =================================================================== //
// Givaro : Prime numbers
//              Modular powering,
//              Fermat numbers,
//              Primality tests, Factorization one by one :
//                      (There are parameters to fix)
// Time-stamp: <30 Aug 00 18:32:49 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#include <math.h>
#include "givintprime.h"
#include "givcompgmp.h"
#include "LinBox/givtimer.h"

// =================================================================== //
// Fermat numbers
// =================================================================== //

FermatDom::Rep& FermatDom::fermat(Rep& f,  const long  n ) const
{
        // fermat(n) = 2^(2^n) + 1
    Rep z;
    assign(z,2);
    pow(f,z, 1 << n);
    return addin(f,one);
}
   
int FermatDom::pepin (const long n) const
{
        // Fermat number primality test
    Rep fn;
    fermat(fn,n);
    Rep y,z,t;
    sub(z,fn,one);
    assign(t,2);
    divin(z,t);
    assign(t,3);
    powmod(y,t,z,fn);
    subin(fn,y);
    return isone(fn);
}


// =================================================================== //
// Primality tests and factorization algorithms
// =================================================================== //

// =================================================================== //
// Primality tests 
// =================================================================== //

template<class RandIter>
int IntPrimeDom::Miller(RandIter& g, const Rep& n) const
{
        // Monte Carlo algorithm
        // returns 1    : n prime with probability 3/4
        // returns 0    : n composite
    if (n < 2) return 0;
    if (n <= 3) return 1;
    IntPrimeDom::Rep r=0,t=n-1,a,q;
    random(g,a,n);
    long s=0;
    while(r == 0) {
        t = t / 2;
        r = t % 2;
        s++;
    }
    powmod(q,a,t,n);
    if ( (q==1) || (q == (n-1))) return 1;
    for(;s>1;s--) {
        q = (q*q) % n;
        if (q == (n-1)) return 1;
    }
    return 0;
}

    
template<class RandIter>
IntPrimeDom::Rep& IntPrimeDom::test_Lehmann(RandIter& g, Rep& r, const Rep& n) const {
        // Monte Carlo algorithm
        // returns n-1  : n prime with probability 1/2
        // returns 1    : n composite with probability 1/2
        // else         : n composite
    IntPrimeDom::Rep A;
    random(g,A,n);
    return powmod(r,A,(n-1)/2,n);
}

template<class RandIter>
int IntPrimeDom::Lehmann(RandIter& g, const Rep& n)  const 
{
    if (n < 2) return 0;
    if (n <= 3) return 1;
    IntPrimeDom::Rep tmp;
    IntPrimeDom::test_Lehmann(g,tmp,n);
    if (tmp == (n-1))
        return 1;
    return 0;
}


IntPrimeDom::Rep& IntPrimeDom::nextprimein(Rep& n)  const {
    if (isleq( n,1)) return n=2;
    Rep::Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
    mpz_tdiv_r_ui( (mpz_ptr)& gmp_res, (mpz_ptr)&(n.gmp_rep), 2);
    if ( mpz_cmp_ui( (mpz_ptr)& gmp_res, 0UL) )
        addin(n,2);
    else
        addin(n,1);
    mpz_clear((mpz_ptr)&gmp_res);
    while (! isprime(n) )
        addin(n,2);
    return n;
}

IntPrimeDom::Rep& IntPrimeDom::nextprime(Rep& n, const Rep& p)  const {
    if (isleq( p,1)) return n=2;
    if (&n == &p) return nextprimein(n);
    if (iszero( mod(n,p,2)))
        add(n,p,1);
    else
        add(n,p,2);
    while (! isprime(n) )
        addin(n,2);
    
    return n;
}

IntPrimeDom::Rep& IntPrimeDom::prevprimein(Rep& n)  const {
    if (isleq( n,2)) return n=2;
    Rep::Rep gmp_res; mpz_init((mpz_ptr)&gmp_res);
    mpz_tdiv_r_ui( (mpz_ptr)& gmp_res, (mpz_ptr)&(n.gmp_rep), 2);
    if ( mpz_cmp_ui( (mpz_ptr)& gmp_res, 0UL) )
        subin(n,2);
    else
        subin(n,1);
    mpz_clear((mpz_ptr)&gmp_res);
    while (! isprime(n) )
        subin(n,2);
    return n;
}

IntPrimeDom::Rep& IntPrimeDom::prevprime(Rep& n, const Rep& p)  const {
    if (isleq( p,2)) return n=2;
    if (&n == &p) return prevprimein(n);
    if (iszero( mod(n,p,2)) )
        sub(n,p,1);
    else
        sub(n,p,2);
    while (! isprime(n) )
        subin(n,2);
    return n;
}


