// =================================================================== //
// Givaro : Prime numbers
//              RSA public-key cipher codes
// Time-stamp: <21 Mar 00 17:47:21 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_RSA_Public_KEY_
#define _GIVARO_RSA_Public_KEY_

#include <iostream.h>
#include <string>
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
//     Rep p = n,mil(1000), r;
//     do {
//         o << char( Integer2long( mod(r,p,mil) ) );
//         divin(p,mil);
//     } while (p>0);

    string st; Integer2string(st,n);
    int j = st.length()-(st.length()/3)*3;
    if (j) 
        o << char( atoi( (st.substr(0,j)).c_str() ) );
    
    for(; j<st.length();j+=3)
        o << char( atoi( (st.substr(j,3)).c_str() ) );
    return o;
}


template<class RandIter>
ostream& IntRSADom<RandIter>::ecriture_Int(ostream& o, const Rep& p, const long size) const {
//     if (p>0) {
//         for(int i=size-log(p);i--;)
//             o << "0";
//         o << p;
//     }
    return o << p << endl;
}


template<class RandIter>
ostream& IntRSADom<RandIter>::encipher(ostream& o, istream& in) const {
    int x;
    char * sx = new char[3]; string st;
    do {
        st = "";
        for(int i=3; i<_lm; i+=3) {
            x = in.get();
            if (in.eof())
                break;
// cerr << "x  : " << x << endl;
            if (x < 10) 
                st += "00";
            else if (x < 100)
                st += '0';
            st += lltostr(x,sx);
        }
// cerr << "sx : " << sx << endl;
// cerr << "st : " << st << endl;
        Rep res(st.c_str()), r;
        ecriture_Int(o, powmod(r, res,_k,_m),_lm);
// cerr << "enc: " << res << " ---> " << r << endl;
    } while (! in.eof());
    delete [] sx;
}           

template<class RandIter>
ostream& IntRSADom<RandIter>::decipher(ostream& o, istream& in) {
//     Rep pas, res, dix(10), r;
//     int * TC = new int[9];
//     for(int j=0;j<10;j++)
//         TC[j] = j;
//     int * TR = &TC[-48];
//     int * envers = new int[_lm];
    string st;
    do {
        in >> st;
        if (in.eof()) break;
        Rep res(st.c_str()),r;
//         res = 0;
//         pas = 1;
// 	int i=0;
//         for(; i<_lm; ++i)
//             envers[i] = TR[in.get()];
//         for(i=1;i<=_lm;++i) {
//             axpyin(res,envers[_lm-i],pas);
//             mulin(pas,dix);
//         }
//         if ( iszero(_u) ) 
//             point_break(_u);
        ecriture_str(o,powmod(r, res,_u,_m));
// cerr << "dec: " << res << " ---> " << r << endl;
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
    
    nextprime(p1, p1);
    nextprime(q1, q1);

    

    Rep d,l,deux(2), ttm;
    
    random(g,d,long(psize/3));
    random(g,l,long(qsize/3));

    if ( isone( mod(ttm, d, deux) ) )
        addin(d,one);
    if ( isone( mod(ttm, l, deux) ) )
        addin(l,one);
    
    Rep p,q;

    do {
        axpy(p,p1,d,one);
        addin(d,deux);
    } while (! (isprime(p) && isprime(p)) );
    
    do {
        axpy(q,q1,l,one);
        addin(l,deux);
    } while (! (isprime(q) && isprime(q)) );
    Rep phim; mul(phim, sub(d,p,one), sub(l,q,one));
    mul(m, p, q);

    Rep v, gd;

    do {
        random(g,k,phim);
    } while (gcd(gd,u,v,k,phim) != 1);
    modin(u,phim);
    if ( islt(u,zero) ) addin(u,phim);
}

// =================================================================== //
// Breaking codes
// =================================================================== //
template<class RandIter>
IntRSADom<RandIter>::Rep& IntRSADom<RandIter>::point_break(Rep& u) {
    if ( iszero(_u) ) {
        Rep p,v,d, pm;
        factor(p, _m);
        mul(pm, sub(v,p,one), subin( div(d,_m,p), one ) );
        gcd(d,_u,v,_k,pm);
        if (islt(_u,zero)) addin(_u, pm);
    }
    return u = _u;
}

#endif
