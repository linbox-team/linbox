// -- free memory allocated in array !
template<template<class> class Container>
inline IntRNSsystem< Container >::~IntRNSsystem()
{}

template<template<class> class Container>
inline IntRNSsystem< Container >::IntRNSsystem ()
        : _primes(0), _ck(0), _prod(one)
{}

template<template<class> class Container>
inline IntRNSsystem< Container >::IntRNSsystem (const IntRNSsystem< Container >& R)
 : _primes(R._primes), 
   _ck(R._primes), _prod(R._prod)
{}


  // -- Array of primes are given
template<template<class> class Container>
inline IntRNSsystem< Container >::IntRNSsystem( const IntRNSsystem< Container >::array& inprimes) 
 : _primes(inprimes),
   _ck(0), _prod(one)
{
   GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem< Container >::IntRNSsystem] bad size of array");
}

  // -- Array of primes are given
template<template<class> class Container> template <class TT>
inline IntRNSsystem< Container >::IntRNSsystem( const Container<TT>& inprimes) 
 : _ck(0), _prod(one)
{
   GIVARO_ASSERT( inprimes.size()>0, "[IntRNSsystem< Container >::IntRNSsystem] bad size of array");
   _primes.resize(inprimes.size());
   typename Container<TT>::const_iterator np = inprimes.begin();
   for(typename array::iterator pi = _primes.begin(); pi != _primes.end(); ++pi, ++np)
       *pi = element( *np );

}

  // -- Product of primes
template<template<class> class Container>
inline void IntRNSsystem< Container >::ComputeProd()
{
    if (isone(_prod))
        for (typename array::const_iterator pi = _primes.begin();pi != _primes.end(); ++pi)
            mulin(_prod, *pi);
}

  // -- Computes Ck , Ck = (\prod_{i=0}^{k-1} primes[i])^(-1) % primes[k],
  // for k=1..(_sz-1)
template<template<class> class Container>
inline void IntRNSsystem< Container >::ComputeCk()
{
  if (_ck.size() !=0) return; // -- already computed

  // - reallocation of a new array :
  size_t size = _primes.size();
  _ck.resize(size);
//  _ck[0] = Neutral::zero; // -- undefined and never used
  _ck[0] = 0UL; // -- undefined and never used

  for (size_t k=1; k < size; ++k)
  {
    element prod = _primes[0];
    for (size_t i= 1; i < k; ++i)
        modin( mulin( prod, _primes[i]), _primes[k]);
    
    element g,u;
    gcd(g,u, _ck[k],_primes[k],prod); // _ck[k] * prod = g mod _primes[k]
  }
}


template<template<class> class Container>
inline const IntRNSsystem< Container >::array& IntRNSsystem< Container >::Primes() const
{ 
  return _primes; 
}


template<template<class> class Container>
inline const IntRNSsystem< Container >::element IntRNSsystem< Container >::ith(const size_t i) const
{
  return _primes[i];
}


template<template<class> class Container>
inline const IntRNSsystem< Container >::element IntRNSsystem< Container >::product() const
{
    ((IntRNSsystem< Container >*)this)->ComputeProd();
    return _prod;
}

template<template<class> class Container>
inline const IntRNSsystem< Container >::array& IntRNSsystem< Container >::Reciprocals() const
{
  if (_ck.size() ==0) ((IntRNSsystem< Container >*)this)->ComputeCk();
  return _ck;
}


template<template<class> class Container>
inline const IntRNSsystem< Container >::element IntRNSsystem< Container >::reciprocal(const size_t i) const
{
  if (_ck.size() ==0) ((IntRNSsystem< Container >*)this)->ComputeCk();
  return _ck[i];
}
