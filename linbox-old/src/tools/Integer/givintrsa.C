// =================================================================== //
// Givaro : Prime numbers
//              RSA public-key cipher codes
// Time-stamp: <21 Mar 00 14:32:31 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_RSA_Public_KEY_
#define _GIVARO_RSA_Public_KEY_

#include <iostream.h>
#include "givinteger.h"
#include "givintrsa.h"

// =================================================================== //
// log[10]
// =================================================================== //

template<class RandIter>
long IntRSADom<RandIter>::log(const Rep& n, const long b = 10) const {
    Rep p = n, log(b);
    long res = 0;
    do {
        divin(p,log);
        ++res;
    } while (! iszero(p) );
    return res;
}   
    

// =================================================================== //
// Text conversions
// =================================================================== //

                
template<class RandIter>
ostream& IntRSADom<RandIter>::ecriture_str(ostream& o, const Rep& n) const {
    Rep p = n,mil(1000), r;
    do {
        o << char( Integer2long( mod(r,p,mil) ) );
//         p = p / 1000;
        divin(p,mil);
    } while (p>0);
    return o;
}


template<class RandIter>
ostream& IntRSADom<RandIter>::ecriture_Int(ostream& o, const Rep& p, const long size) const {
    if (p>0) {
        for(int i=size-log(p);i--;)
            o << "0";
        o << p;
    }
    return o;
}


template<class RandIter>
ostream& IntRSADom<RandIter>::encipher(ostream& o, istream& in, const Rep& m, const Rep& k) const {
    long lm = log(m);
    int x;
    Rep pas, res, mil(1000), r;
    do {
        res = 0;
        pas = 1;
//         for(int i = 0; i<lm/3; i++) {
        for(int i = lm/3; i--;) {
            x = in.get();
            if (in.eof())
                break;
            axpyin(res,x,pas);
            mulin(pas,1000);
//             res += x * pas;
//             pas *= 1000;
        }
        ecriture_Int(o, powmod(r, res,k,m),lm);
    } while (! in.eof());
}           

template<class RandIter>
ostream& IntRSADom<RandIter>::decipher(ostream& o, istream& in, const Rep& m, const Rep& u) const {
    long lm = log(m);
    Rep pas, res, dix(10), r;
    int * TC = new int[9];
    for(int j=0;j<10;j++)
        TC[j] = j;
    int * TR = &TC[-48];
    int * envers = new int[lm];
    do {
        res = 0;
        pas = 1;
	int i=0;
        for(; i<lm; ++i)
            envers[i] = TR[in.get()];
        for(i=1;i<=lm;++i) {
            axpyin(res,envers[lm-i],pas);
            mulin(pas,dix);
//             res += envers[lm-i]*pas;
//             pas *= 10;
        }
        ecriture_str(o,powmod(r, res,u,m));
    } while (! in.eof());
}           

// =================================================================== //
// Keys generation
// public keys are m and k, the secret key is u.
// ciphering is computing	: x^k mod m
// deciphering is computing	: b^u mod m
// since for any x, x^(k.u) = x mod m
// =================================================================== //
template<class RandIter>
void IntRSADom<RandIter>::keys_gen(random_generator& g, long psize, long qsize, Rep& m, Rep& k, Rep& u) const {
    Rep p1,q1;
    random(g,p1,long((2*psize)/3));
    random(g,q1,long((2*qsize)/3));
    cerr << "p1: " << p1 << endl;
    cerr << "q1: " << q1 << endl;
    
    nextprime(p1, p1);
    nextprime(q1, q1);
    cerr << "p1: " << p1 << endl;
    cerr << "q1: " << q1 << endl;

    

    Rep d,l,deux(2), ttm;
    
    random(g,d,long(psize/3));
    random(g,l,long(qsize/3));

    if ( isone( mod(ttm, d, deux) ) )
        addin(d,one);
//         d = d + 1;
    if ( isone( mod(ttm, l, deux) ) )
        addin(l,one);
//         l = l + 1;
    
    Rep p,q;

    do {
//         p = p1 * d + 1;
        axpy(p,p1,d,one);
//         d = d + 2;
        addin(d,deux);
    } while (! (isprime(p) && isprime(p)) );
    
    do {
//         q = q1 * l + 1;
        axpy(q,q1,l,one);
        addin(l,deux);
//         l = l + 2;
    } while (! (isprime(q) && isprime(q)) );
    cerr << "p: " << p << endl;
    cerr << "q: " << q << endl;
    
//     Rep phim = (p-1)*(q-1);
    Rep phim; mul(phim, sub(d,p,one), sub(l,q,one));
//     m = p*q;
    mul(m, p, q);

    Rep v, gd;

    do {
        random(g,k,phim);
    } while (gcd(gd,u,v,k,phim) != 1);
    modin(u,phim);
//     if (u < 0) u = u + phim;
    if ( islt(u,zero) ) addin(u,phim);
}

// =================================================================== //
// Breaking codes
// =================================================================== //
template<class RandIter>
IntRSADom<RandIter>::Rep& IntRSADom<RandIter>::point_break(Rep& u, const Rep& m, const Rep& k) const {
    Rep p,v,d;
    factor(p, m);
    gcd(d,u,v,k,(p-1)*((m/p)-1));
    return u;
}

#endif
