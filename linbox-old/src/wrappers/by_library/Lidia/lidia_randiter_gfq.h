/* File: src/wrapper/by_scope/field/LIDIA_randiter.h
 * Author: Pascal Giorgi for the LinBox group
 */

#ifndef _LIDIA_RANDITER_GFQ_
#define _LIDIA_RANDITER_GFQ_


namespace LinBox
{
using namespace LiDIA;

  


 template<class field> class lidia_randIter_gfq
    {
    public:
      

      typedef gf_element element;

       lidia_randIter_gfq(const field& F, 
		     const integer& size = 0, 
		     const integer& seed = 0)
      : _size(size), _seed(seed) , GF(F)
    { 
      if (_seed == integer(0)) _seed = time(NULL);
      
      integer cardinality ;
      F.cardinality(cardinality);
      if ( (cardinality != integer(-1)) &&  (_size > cardinality) )
	_size = cardinality;


#ifdef TRACE
      cout << "created random generator with size " << _size 
	   << " and seed " << _seed << endl;
#endif // TRACE

      // Seed random number generator
      srand(_seed);

    }


 lidia_randIter_gfq(const lidia_randIter_gfq& R)
      : _size(R._size), _seed(R._seed) {}


 ~lidia_randIter_gfq(void) {}


 lidia_randIter_gfq& operator=(const lidia_randIter_gfq& R)
    {
      if (this != &R) // guard against self-assignment
      {
	_size = R._size;
	_seed = R._seed;
      }

      return *this;
    }


 element& operator() (void)
    {
     element e(GF);
     e.randomize();
     return *(new element(e));
    }

 lidia_randIter_gfq(void) : _size(0), _seed(0) { time(NULL); }


private:

    /// Sampling size
    integer _size;
    
    /// Seed
    integer _seed;

    field GF;

  }; // class lidia_randIter_gfq

} // namespace LinBox

#endif // _LIDIA_RANDITER_GFQ_
