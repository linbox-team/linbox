  // -- Computation of a mixed-radix representation of the residu.
template<template<class> class Container> template<class TT>
inline void IntRNSsystem< Container >::RnsToMixedRadix
  (IntRNSsystem< Container >::array& mixrad, const Container<TT>& residu) const
{
  long i,j;
  size_t size = _primes.size();
  if (mixrad.size() < size) mixrad.resize( size );
  
  // -- Computation of  Ck 
  if (_ck.size()==0) ((IntRNSsystem*)this)->ComputeCk();

  // -- size-1 steps
  element tmp;
  mixrad[0] = residu[0];
  for (i=1; i < size; i++)
  {  // - computes pp_i = r_0 + r_1*p_0 + ... + r_{i-1} \prod_{j<i-2} p_j [p_i]
     // Horner scheme
     tmp = mixrad[i-1];
     for (j= i-2; j>=0; j--) {
// cerr << tmp << " * " << _primes[j] << " + " << mixrad[j] << " mod " << _primes[i] << " = ";
         modin( addin( mulin(tmp, _primes[j]), mixrad[j]), _primes[i]);
// cerr << tmp << ";#Horner scheme" << endl;
     }
     // - m_i = (r_i - pp_i)*ck_i, ck is reciprocals
// cerr << "(" << residu[i] << " - " << tmp << ") * " << _ck[i] << " mod " << _primes[i] << " = ";
     mod(mixrad[i],mulin(sub(tmp,  residu[i], tmp),_ck[i]) , _primes[i] );
// cerr << mixrad[i] << ";#mixrad" << endl;
  }
  
}
 
  

  // -- Convert a mixed radix representation to an Integer
template<template<class> class Container>
inline void IntRNSsystem< Container >::MixedRadixToRing( element& res, const IntRNSsystem< Container >::array& mixrad ) const 
{
  size_t size = _primes.size();
  if (size != mixrad.size()) 
    throw GivError("[IntRNSsystem::MixedRadixToRing]: bad size of input array");
  res = mixrad[size-1];
  for (int i=size-2; i>=0; --i) {
      addin( mulin(res, _primes[i]), mixrad[i]);
//     res *= _primes[i];
//     res += mixrad[i];
  }
}


  // Convert an integer to a RNS representation (which is given by this)
template<template<class> class Container>
inline void IntRNSsystem< Container >::RingToRns( IntRNSsystem< Container >::array& rns , const external& a) const
{
  size_t size = _primes.size();
  if (rns.size() != size) rns.resize(size);
  // -- may be faster using the recursive 
  // tree algorithm a mod p_1...p_k/2, and a mod p_k/2+1...p_k
  for (int i=0; i<size; i++) 
      mod( rns[i], a, _primes[i]);
//     rns[i] = mod(a, _primes[i]);
}

  // Convert to an Integer:
template<template<class> class Container> template<class TT>
inline void IntRNSsystem< Container >::RnsToRing( external& I, const Container<TT>& rns) const 
{
  // - Computation of a mixed radix representation of this
  IntRNSsystem< Container >::array mixrad(_primes.size());
  RnsToMixedRadix( mixrad , rns );

  // - Convert mixrad to an integer
  MixedRadixToRing( I, mixrad ) ;
  return;
}

