// =================================================================== //
// Givaro : Euler's phi function
//          Primitive roots.
// Needs list structures : stl ones for instance
// Time-stamp: <31 Aug 00 22:38:20 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#include "givintnumtheo.h"
#include <list.h>

// =================================================================== //
// Euler's phi function
// =================================================================== //
template<class RandIter>
IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::phi(Rep& res, const Rep& n) const {
    if (isleq(n,1)) return res=n;
    if (isleq(n,3)) return sub(res,n,one);
    list<Rep> Lf;
    set(Lf,n);
    return phi(res,Lf,n);
}


template<class RandIter>
template< template<class> class Container> IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::phi(Rep& res, const Container<Rep>& Lf, const Rep& n) const {
    if (isleq(n,1)) return res=n;
    if (isleq(n,3)) return sub(res,n,one);
    res = n; Rep t,m;
    for(typename Container<Rep>::const_iterator f=Lf.begin(); f!=Lf.end(); ++f) 
        mul(res, divexact(t,res,*f), sub(m, *f, one));
    return res;
}

// =================================================================== //
// Möbius function
// =================================================================== //
template<class RandIter>
template< template<class> class Container> short IntNumTheoDom<RandIter>::mobius(const Container<unsigned long>& lpow) const {
    if (lpow.size()) {
        short mob = 1;
        for(typename Container<unsigned long>::const_iterator i=lpow.begin();i != lpow.end(); ++i) {
            if (*i > 1) {
                 return 0;
            } else
                mob = -mob;
        }
        return mob;
    } else
        return 1;
}
    
template<class RandIter>
short IntNumTheoDom<RandIter>::mobius(const Rep& a) const {
    list<Rep> lr;
    list<unsigned long> lp;
    set(lr, lp, a);
    return mobius(lp);
}


// =================================================================== //
// Primitive Root
// =================================================================== //

template<class RandIter>
IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root(Rep& A, int& runs, const Rep& n) const {
        // n must be in {2,4,p^m,2p^m} where p is an odd prime
        // else infinite loop
    if (isleq(n,4)) return sub(A,n,one);
    if (iszero(mod(A,n,4))) return A=zero;
    Rep p,ismod2, q, no2, root; 
    if (iszero(mod(ismod2,n,2))) divexact(no2,n,2); else no2=n;
    p=no2;
    bool isp; int k = 1; 
    while (! isprime(p) ) {
        sqrt(root, p);
        while (mul(q,root,root) == p) {
            p = root;
            sqrt(root,p);
        }
        if (! isprime(p) ) {
            q=p;
            while( p == q ) factor(p, q);
            divin(q,p);
            if (q < p) p = q;
        }
    }
    if (iszero(ismod2)) mul(q,p,2); else q=p;
    for(;q != n;++k,q*=p);
    Rep phin, tmp; 
    phi(phin,p);
    list<Rep> Lf;
    set(Lf,phin);
    list<Rep>::iterator f;
    for(f=Lf.begin();f!=Lf.end();++f)
            div(*f,phin,*f);
    int found; runs = 0;
    A=2;
    found = ++runs;
    for(f=Lf.begin();(f!=Lf.end() && found);f++)
        found = (! isone( powmod(tmp,A,*f,p)) );
    if (! found) {
        A=3;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (! found) {
        A=5;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (! found) {
        A=6;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    while (! found) {
       do {
            random(_g, A, p);
            addin( modin(A,sub(tmp,p,7)) , 7);
        } while ( ! isone(gcd(tmp,A,p)) );
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (k == 1) {
        if (iszero(ismod2) && iszero(mod(ismod2, A, 2)))
            return A+=p;
        else
            return A;
    } else {
        if (! is_prim_root(A,no2))
            A+=p;
        if (iszero(ismod2) && iszero(mod(ismod2, A, 2)))
            return A+=no2;
        else
            return A;
    }
}

template<class RandIter>
IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root(Rep& A, const Rep& n) const { int runs; return prim_root(A, runs, n); }

    
        
template<class RandIter>
IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::lowest_prim_root(Rep& A, const Rep& n) const {
        // n must be in {2,4,p^m,2p^m} where p is an odd prime
        // else returns zero
    if (isleq(n,4)) return sub(A,n,one);
    if (iszero(mod(A,n,4))) return A=zero;
    Rep phin, tmp; 
    phi(phin,n);
    list<Rep> Lf;
    set(Lf,phin);
    list<Rep>::iterator f;
//             *f = phin / (*f);
    for(f=Lf.begin();f!=Lf.end();++f)
            div(*f,phin,*f);
    int found=0;
    for(A = 2;(isleq(A,n) && (! found));addin(A,1)) {
        if (isone(gcd(tmp,A,n))) {
            found = 1;
            for(f=Lf.begin();(f!=Lf.end() && found);f++)
                found = (! isone( powmod(tmp,A,*f,n)) );
//                 found = ( powmod(A,*f,n) != 1);
        }
    }
    if (isleq(A,n))
        return subin(A,1);
    else
        return A=zero; 
}

template<class RandIter>
int IntNumTheoDom<RandIter>::is_prim_root(const Rep& p, const Rep& n) const {
        // returns 0 if failed
    int found=0;
    Rep phin, tmp; 
    phi(phin,n);
    list<Rep> Lf;
    set(Lf,phin);
    list<Rep>::iterator f=Lf.begin();
    Rep A; mod(A,p,n);
    if (isone(gcd(tmp,A,n))) {
        found = 1;
        for(;(f!=Lf.end() && found);f++) {
//             found = ( powmod(A,phin / (*f),n) != 1);
            found = (! isone( powmod(tmp,A, div(tmp,phin,*f),n)) );
        }
    }
    return found;
}

template<class RandIter>
int IntNumTheoDom<RandIter>::isorder(const Rep& g, const Rep& p, const Rep& n) const {
        // returns 1 if p is of order g in Z/nZ
    Rep tmp;
    return (isone( pow(tmp, p, Integer2long(g)) ) && isequal( g, order(tmp,p,n) ) );
}


template<class RandIter>
IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::order(Rep& g, const Rep& p, const Rep& n) const {
        // returns 0 if failed
    Rep A; mod(A,p,n);
    if (iszero(A))
	return g = zero;
    if (isone(A))
	return g = one;
    int primroot=0;
    Rep phin,gg,tmp;
    phi(phin,n);
    list<Rep> Lf;
    set(Lf,phin);
    Lf.sort();
    list<Rep>::iterator f=Lf.begin();
    if (isone(gcd(tmp,A,n))) {
        primroot = 0;
        for(;f!=Lf.end();++f)
//             if (! (primroot = ( powmod(A, g = (phin / (*f)),n) != un) ) )
            if ( primroot = isone(powmod(tmp,A, div(g,phin,*f),n)) )
                break;
        if (primroot) {
            for(;f!=Lf.end();++f)
                while (iszero(mod(tmp,g,*f)) && isone( powmod(tmp,A, div(gg,g,*f),n) ) )
                    g = gg;
            return g;
        } else
            return g=phin;
    }
    return g=zero;
}
