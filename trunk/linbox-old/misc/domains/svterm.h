// class svterm wraps value,index pairs.  
// For possible use in linbox vectors. -bds 6/00
#ifndef __SVTERM__
#define __SVTERM__
#include <pair.h>
#include "Integers.h"
//typedef long Integer;

template < class entry >
class svterm : public pair< entry, Integer >
{public:
  svterm(): pair<entry, Integer>(){}
  svterm(entry v, Integer i): pair<entry, Integer>(v,i){}
  inline const entry& value()const{ return first; }
  inline entry& value(){ return first; }
  inline const Integer& index()const{ return second; }
  inline Integer& index(){ return second; }
};
#endif

/********* some test code ***********
#include <iostream>
main()
{
  //svterm<int> t(2,Integer(INIT_VAL, 3));
  //svterm<int> t(2,Integer(3));
  svterm<int> t;
  t.value() = 5;
  t.index() = 3;
  cout << "( " << t.value() << ", " << t.index() << " )" << endl;

  t.index() = 4;
  cout << "( " << t.value() << ", " << t.index() << " )" << endl;

  const svterm<int> u =* new svterm<int>(5,Integer(INIT_VAL, 5)) ;
  cout << "( " << u.value() << ", " << u.index() << " )" << endl;
  /* error (as it should be): */ 
  // u.index() = 3;

}
*********/
