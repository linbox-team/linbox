// =================================================================== //
// Givaro : Prime numbers
//              Modular powering,
//              Fermat numbers,
//              Primality tests, Factorization :
//                      (There are parameters to fix)
// Time-stamp: <30 Aug 00 16:58:00 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef _GIVARO_INTEGERS_PRIME_H_
#define _GIVARO_INTEGERS_PRIME_H_

#include "givinteger.h"


// =================================================================== //
// Fermat numbers
// =================================================================== //
class FermatDom : public IntegerDom {
public:
    FermatDom() : IntegerDom() {}
    Rep& fermat (Rep&, const long)  const ;
    int pepin (const long) const ;
};


// =================================================================== //
// Primality tests and factorization algorithms
// =================================================================== //

// Those macros are parameters to fix

// primes known
// first array
#define LOGMAX 3512
#define TABMAX 32768
// second array
#define LOGMAX2 3031
#define TABMAX2 65536
// Bounds between big and small
#define BOUNDARY_isprime TABMAX    
#define BOUNDARY_2_isprime TABMAX2    

// =================================================================== //
// Primality tests 
// =================================================================== //
class IntPrimeDom : public IntegerDom {
public:
    IntPrimeDom() :  IntegerDom() {}

    int isprime(const Rep& n, int r=1) const 
        {
/*
  return probab_prime(n);
*/
//             return ((n)<BOUNDARY_isprime ?  isprime_Tabule(n) : 
//                     (n)<BOUNDARY_2_isprime ? isprime_Tabule2(n) : 
//                     probab_prime(n));
            long l;
            return (islt(n,BOUNDARY_isprime) ?  isprime_Tabule(access(l,n)): 
                    islt(n,BOUNDARY_2_isprime) ? isprime_Tabule2(access(l,n)): 
                    local_prime(n,r));
        }

    template<class RandIter>
    int Miller(RandIter& g, const Rep& n) const  ;
    template<class RandIter>
    Rep& test_Lehmann(RandIter& g, Rep&, const Rep& n) const  ;
    template<class RandIter>
    int Lehmann(RandIter& g, const Rep& n)  const ;
    int isprime_Tabule(const int n) const ;
    int isprime_Tabule2(const int n) const ;
    


    Rep& nextprime(Rep&, const Rep&) const ;
    Rep& prevprime(Rep&, const Rep&) const ;
    Rep& nextprimein(Rep&) const ;
    Rep& prevprimein(Rep&) const ;


// Using Integer
    int local_prime(const Rep& n, int r=1) const { return probab_prime(n,r); }
    int& access(int& r, const Rep& a) const { return r=Integer2long(a); }
    long& access(long& r, const Rep& a) const { return r=Integer2long(a); }
    double& access(double& r, const Rep& a) const { return r=Integer2double(a); }

private:
    static int IP[LOGMAX+5];  // -- table for Tabule
    static const int * TP;    // -- shifted table 
    static int IP2[LOGMAX2+5]; // -- table for Tabule2
    static const int * TP2;    // -- shifted table 
/*
  static int Tabule2(const Integer& p) ;
  static int Tabule(const Integer& p) ;
  static int _memTab2[LOGMAX2+5];   // -- table for Tabule2
  static const int* _Tab2; // -- shifted _memTabule2
  static int _memTab[];    // -- table for Tabule
  static const int* _Tab;  // -- shifted _memTabule
*/
};

#endif _GIVARO_INTEGERS_PRIME_H_
