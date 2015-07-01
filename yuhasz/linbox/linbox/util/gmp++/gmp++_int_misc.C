// ==========================================================================
// $Source$
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id$
// ==========================================================================
// Description: 

#include <iostream>
#include "gmp++_int.h"

//-------------------------------------------fact (unsigned long l)
Integer fact ( unsigned long l) 
{
  Integer Res ;
  mpz_fac_ui( (mpz_ptr)&(Res.gmp_rep), l ) ;
  return Res ;
}

//-------------------------------------------square root
Integer sqrt(const Integer &a)
{
  Integer q;
  mpz_sqrt( (mpz_ptr)&(q.gmp_rep),
              (mpz_ptr)&(a.gmp_rep)) ;
  return q;
}

Integer sqrt(const Integer &a, Integer& r)
{
  Integer q;
  mpz_sqrtrem( (mpz_ptr)&(q.gmp_rep),
              (mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(a.gmp_rep)) ;
  return q;
}

// base p logarithm of a
long logp(const Integer& a, const Integer& p) {
    std::list< Integer > pows;
    Integer puiss = p, sq;
    do {
        pows.push_back( puiss );
    } while ( (puiss *= puiss) <= a );
    puiss = pows.back(); pows.pop_back();
    long res = (1 << pows.size());
    while (! pows.empty() ) {
        if ((sq = puiss * pows.back()) <= a) {
            puiss = sq;
            pows.pop_back();
            res += (1 << pows.size());
        } else
            pows.pop_back();
    }
    return res;
}
    

//------------------------------------------GMP isprime
//     If this function returns 0, OP is definitely not prime.  If it
//     returns 1, then OP is `probably' prime.  The probability of a
//     false positive is (1/4)^r.  A reasonable value of r is 25.

Integer& nextprime(Integer& r, const Integer &p)
{
  mpz_nextprime ((mpz_ptr)&(r.gmp_rep), (mpz_ptr)&(p.gmp_rep)) ;
  return r;
}
int probab_prime(const Integer &p)
{
  return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),1) ;
}

int probab_prime(const Integer &p, int r)
{
  return mpz_probab_prime_p ((mpz_ptr)&(p.gmp_rep),r) ;
}

// ==========================================================================
// Computes and returns the Jacobi and Legendre symbols (u/v) of the integers u and v.  
// The algorithm used is Gmp's.
int jacobi(const Integer& u, const Integer& v)
{
  return mpz_jacobi ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}

int legendre(const Integer& u, const Integer& v)
{
  return mpz_legendre ((mpz_ptr)&(u.gmp_rep),(mpz_ptr)&(v.gmp_rep)) ;
}



//--------------------------------------------Integer::operator <<   // shift left
Integer Integer::operator << (int l) const 
{ return this->operator<<( (unsigned long)l ); }
Integer Integer::operator << (unsigned int l) const 
{ return this->operator<<( (unsigned long)l ); }
Integer Integer::operator << (long l) const 
{ return this->operator<<( (unsigned long)l ); }

Integer Integer::operator << (unsigned long l) const 
{ 
	Integer tmp;
	mpz_mul_2exp((mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return tmp; 
}


//--------------------------------------------Integer::operator >>   // shift right
Integer Integer::operator >> (int l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (long l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (unsigned int l) const
{ return this->operator>>( (unsigned long)l ); }

Integer Integer::operator >> (unsigned long l) const
{ 
	Integer tmp;
	mpz_tdiv_q_2exp( (mpz_ptr)&(tmp.gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return tmp; 
}

//--------------------------------------------Integer::operator <<=   // shift left
Integer& Integer::operator <<= (int l) 
{ return this->operator<<= ( (unsigned long)l ); }
Integer& Integer::operator <<=  (unsigned int l) 
{ return this->operator<<= ( (unsigned long)l ); }
Integer& Integer::operator <<= (long l) 
{ return this->operator<<= ( (unsigned long)l ); }

Integer& Integer::operator <<= (unsigned long l)
{ 
	mpz_mul_2exp((mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return *this; 
}


//--------------------------------------------Integer::operator >>=   // shift right
Integer& Integer::operator >>= (int l)
{ return this->operator>>= ( (unsigned long)l ); }
Integer& Integer::operator >>= (long l) 
{ return this->operator>>= ( (unsigned long)l ); }
Integer& Integer::operator >>= (unsigned int l) 
{ return this->operator>>= ( (unsigned long)l ); }

Integer& Integer::operator >>= (unsigned long l) 
{ 
	mpz_tdiv_q_2exp( (mpz_ptr)&(gmp_rep), (mpz_srcptr)&(gmp_rep), l );
	return *this; 
}

//------------------------------------------- convert method
//------------------------------------------- casting method
long Integer2long  ( const Integer& n)
{
  return mpz_get_si ( (mpz_srcptr)&n.gmp_rep);
}
double Integer2double( const Integer& n)
{
  return mpz_get_d( (mpz_srcptr)&n.gmp_rep);
}

Integer::operator int() const {
	return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}
Integer::operator unsigned int() const {
	return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}
Integer::operator long() const {
	return mpz_get_si ( (mpz_srcptr)&gmp_rep);
}
Integer::operator unsigned long() const {
	return mpz_get_ui ( (mpz_srcptr)&gmp_rep);
}
#ifdef __USE_GMPPLUSPLUS_64__
Integer::operator unsigned long long() const {
	unsigned long low = (unsigned long)(*this);
	Integer rem;
	mpz_tdiv_q_2exp( (mpz_ptr)&(rem.gmp_rep), (mpz_srcptr)&(gmp_rep), CHAR_BIT*sizeof(unsigned long int) );
	unsigned long high = (unsigned long)(rem);
	unsigned long long tmp = high;
	tmp <<= CHAR_BIT*sizeof(unsigned long int) ;
	return tmp += low;
}
Integer::operator long long() const {
	unsigned long long tmp = (unsigned long long)(*this); 
	return (long long)tmp;
}
#endif

Integer::operator double() const {
	return mpz_get_d ( (mpz_srcptr)&gmp_rep);
}
Integer::operator float() const {
	return (float)mpz_get_d ( (mpz_srcptr)&gmp_rep);
}

