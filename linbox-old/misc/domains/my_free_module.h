#include "LB_vector.h"
#ifndef _MY_FREE__
#define _MY_FREE__

//#include "interfaces.h"
template<class Ring, const long dim>
class my_free_module //: public FreeModule_a<Ring>
{public:
  typedef Ring entry_domain;
  typedef typename entry_domain::element entry;
  class element : public SD_vector<entry_domain, zero(), SD_dense_vector_tag> 
  {public:
    element() : SD_vector<entry_domain, entrySD_dense_vector_tag> (dim) {}
  };
  
 public: //??? protected?
  entry_domain R;

 public:
  my_free_module(entry_domain r) { R = r; }

  entry& dotprod(entry& z, const element& u, const element& v, const string& s = "sparsedense") const
  { z = R.zero();
    if (s == "densedense")
    { typename element::iterator p, q;
      for (p = u.begin(), q = v.begin(); p < u.end(); p++, q++)
      {
        entry r;
        R.mul(r, *p, *q);
        R.add(z, z, r);
      }
    }
    else if (s == "sparsesparse")
    {
    }
    else // presumably (s == "sparsedense")
    {
      typename element::sparse_iterator p; // should be const
      for (p = u.s_begin(); p < u.s_end(); p++)
      { 
	entry r;
	R.mul(r, *p, v[p.index()]);
        R.add(z, z, r);
      }
    }
    return z;
  }
}; // class my_free_module

#endif
