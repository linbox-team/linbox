#include <linbox/field/unparametric.h>
#include <linbox/randiter/unparametric.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include <time.h>

namespace LinBox
{
  template<>
  class UnparametricRandIter<NTL::ZZ_pE>
    {
    public:
      typedef NTL::ZZ_pE Element;
      UnparametricRandIter<NTL::ZZ_pE>(const UnparametricField<NTL::ZZ_pE>& F =UnparametricField<NTL::ZZ_pE>(), 
				       const size_t& size = 0,
				       const size_t& seed = 0
				       )
	: _size(size), _seed(seed)
	{
	  if(_seed == 0)
	    NTL::SetSeed(NTL::to_ZZ(time(0)));
	  else
	    NTL::SetSeed(NTL::to_ZZ(_seed));
	}
      
      UnparametricRandIter<NTL::ZZ_pE>(const UnparametricRandIter<NTL::ZZ_pE>& R)
	: _size(R._size), _seed(R._seed) 
	
	{
	  if(_seed == 0)
	    NTL::SetSeed(NTL::to_ZZ(time(0)));
	  else
	    NTL::SetSeed(NTL::to_ZZ(_seed));
	}
      Element& random (Element& x)
	{
	   NTL::random(x);
	   return x;
	}
    protected:
      size_t _size;
      size_t _seed;
    };
}

namespace LinBox
{
  template<>
    NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::init (NTL::ZZ_pE &x, const integer &y) const
    {
      x=NTL::to_ZZ_pE(static_cast<long>(y));
      return x;
    }
  
  template<>
    bool UnparametricField<NTL::ZZ_pE>::isZero (const NTL::ZZ_pE& a) const
    {
      return NTL::IsZero(a);
    }
  
  template<>
    bool UnparametricField<NTL::ZZ_pE>::isOne (const NTL::ZZ_pE& a) const
    {
      return NTL::IsOne(a);
    }
  
  
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::convert(integer& c, const NTL::ZZ_pE& e) const
    {
      //fix me. I am not sure what I need to return
      return c;
    }
  
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::characteristic (integer &c) const
    {
      return c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
      //NTL::ZZ_p::modulus();
    }
  
  template<>
    integer& UnparametricField<NTL::ZZ_pE>::cardinality(integer& c) const
    {
      c=static_cast<integer>(to_long(NTL::ZZ_p::modulus()));
      c=pow(c,NTL::ZZ_pE::degree());
      return c;
    }
  NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::inv(NTL::ZZ_pE& x, const NTL::ZZ_pE& y) const
    {
      x=NTL::to_ZZ_pE(1)/y;
      return x;
    }
   NTL::ZZ_pE& UnparametricField<NTL::ZZ_pE>::invin(NTL::ZZ_pE& x) const
     {
       x=NTL::to_ZZ_pE(1)/x;
       return x;
     }
}
