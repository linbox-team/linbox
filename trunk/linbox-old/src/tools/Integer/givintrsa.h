// =================================================================== //
// Givaro : RSA scheme.
// Time-stamp: <21 Mar 00 15:58:14 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_RSA_
#define _GIVARO_RSA_

#include <iostream.h>
#include "givinteger.h"
#include "givintprime.h"
#include "givintfactor.h"

// =================================================================== //
// RSA public-key cipher codes
// =================================================================== //

template<class RandIter = Random>
class IntRSADom : public IntFactorDom<RandIter> {
    Rep _m, _k;
    Rep _u;
    long _lm;
public:

    IntRSADom(RandIter& g = *(new RandIter()) ) { keys_gen(g, 100, 90, _m, _k, _u); _lm = log(m,10); }
    IntRSADom(const long p, const long q, RandIter& g = *(new RandIter()) )  { keys_gen(g, p, q, _m, _k, _u); _lm = log(m,10); }
    IntRSADom(const Rep& m, const Rep& k, const Rep& u) : _m(m), _k(k), _u(u), _lm(log(m,10)) {}
    IntRSADom(const Rep& m, const Rep& k) : _m(m), _k(k), _u(0), _lm(log(m,10)) {}

// =================================================================== //
// Accesses
// =================================================================== //
    const Rep& getm() const { return _m; }
    const Rep& getk() const { return _k; }

// =================================================================== //
// Text conversions
// =================================================================== // 
    ostream& encipher(ostream&, istream&) const ;
    ostream& encipher(ostream& o, istream& in, const Rep& m, const Rep& k) const ;
 
    ostream& decipher(ostream&, istream&) ;

protected:
// =================================================================== //
// Keys generation
// public keys are m and k, the secret key is u.
// ciphering is computing       : x^k mod m
// deciphering is computing     : b^u mod m
// since for any x, x^(k.u) = x mod m
// =================================================================== //

// =================================================================== //
// Here m = p*q
// p and q are prime numbers of respective sizes psize, qsize
// Moreover p-1 and q-1 have one prime factor of respective size 2/3
// since k.u = 1 mod (p-1)(q-1)
// =================================================================== //
    void keys_gen(random_generator& g, long psize, long qsize, Rep& m, Rep& k, Rep& u) const ;
    
// =================================================================== //
// log[10]
// =================================================================== //
    long log(const Rep& n, const long) const ;

// =================================================================== //
// Text conversions
// =================================================================== // 
    ostream& ecriture_str(ostream&, const Rep&) const ;
    ostream& ecriture_Int(ostream&, const Rep&, const long) const ;

public:
// =================================================================== //
// Breaking codes : finding u knowing only m an k ...
// =================================================================== //
    Rep& point_break(Rep& u) ;

};

#include "givintrsa.inl"

#endif
