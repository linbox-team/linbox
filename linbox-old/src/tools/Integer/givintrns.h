#ifndef _ARITHMODU_INTRNS_H
#define _ARITHMODU_INTRNS_H
// ==========================================================================
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// Time-stamp: <06 Apr 00 20:06:10 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
// Description:
//  Modular arithmetic for GIVARO. Here is defined arithmetic functions
//  on rns representation with Givaro Integers.

#include "giverror.h"
#include "givinteger.h"

  // ---------------------------------------------  class RNSsystem
  // Structure which manages list of primes in order to do 


template< template<class> class Container>
class IntRNSsystem : public IntegerDom {
public:
//     typedef element    Ring;
//     typedef element   Modulo;
    typedef element   external;
    typedef Container< element > array;
    

        // Default Cstor, Dstor/Cstor of recopy: 
    IntRNSsystem();
    ~IntRNSsystem(); 
    IntRNSsystem(const IntRNSsystem& R); 

        // -- Cstor with given primes 
    IntRNSsystem( const array& primes );

    template<class TT>
    IntRNSsystem( const Container< TT > & primes );

        // -- Computation of a mixed-radix representation of the residus.
//     void RnsToMixedRadix(array&  mixrad, const array&  residu) const; 
    template<class TT>
    void RnsToMixedRadix(array&  mixrad, const Container<TT>&  residu) const; 

        // -- Convert a mixed radix representation to an external
    void MixedRadixToRing( external& res,  const array& mixrad ) const;

        // -- Convert an Ring element to a its representation
        // with the "this" rns system.
    void RingToRns( array& residu, const external& a ) const;

        // -- Fast conversion: requires pre-computation (first time it was called)
    void fastRingToRns( array& residu, const external& a ) const;

        // -- Convert a representation to an external element
    template<class TT>
    void RnsToRing( external& a, const Container<TT>& residu ) const;

        // -- Fast conversion: requires pre-computation (first time it was called)
    void fastRnsToRing( external& a, const array& residu ) const;

        // ------------- Access methods
 
        // -- Returns the number of primes of this ctxt
    int NumOfPrimes() const { return _primes.size(); } 

        // -- Returns a array to the begin of the array of primes
    const array& Primes() const;
        // -- Returns the ith primes of the rns system
    const element ith(const size_t i) const;

        // -- Returns a array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
    const array& Reciprocals() const;
    const element reciprocal(const size_t i) const;
    const element product() const;

protected:
        // -- Compute some fields of the structure :
    void ComputeCk();

        // -- Compute product of primes
    void ComputeProd();

        // -- Compute the Qk for Ring -> RNS, allocate U
    void ComputeQk();

    array  _primes; 	// - array of the relatively primes numbers
    element _prod;      // - product of primes
    array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes 

        // -- for fast conversion
    size_t _sizek;
    size_t _log2k;
    array  _qk;	// - cf algo Aho, Hopcroft & Ullman
    array  _u;	// - cf algo Aho, Hopcroft & Ullman
};


#include "givintrns_cstor.inl"
#include "givintrns_convert.inl"

#endif
