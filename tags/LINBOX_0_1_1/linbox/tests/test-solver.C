#include <iostream>
#include <linbox/field/ntl-zz_p.h>
#include <vector>
#include <iterator>
#include <linbox/blackbox/diagonal.h>
#include <algorithm>
#include <linbox/solutions/solver.h>
#include <linbox/field/vector-domain.h>
#include <linbox/solutions/minpoly.h>
#include <linbox/util/commentator.h>
int main()
{
  typedef LinBox::UnparametricField<NTL::zz_p> Field;
  const int prime=1073741789;
  NTL::zz_p::init(prime);
  typedef std::vector<NTL::zz_p> Vector;
  
  Field field;

  Vector a(10);
 
  LinBox::VectorDomain<Field> _VD(field);
  int iteration =1000;
  LinBox::commentator.start ("Testing solver in diagonal matrix", "test-solver",iteration);
  for(int i=0;i<iteration;++i)
    {
      for(int i=0;i<a.size();++i)
	a[i]=NTL::to_zz_p(rand()%4+1);
      LinBox::Diagonal<Field, Vector> BB(field,a);
      Vector x(10),y(10);
      for(int i=0;i<y.size();++i)
	y[i]=NTL::to_zz_p(rand()%(prime-1)+1);
      
      LinBox::solver(BB,x,y,field);
      BB.applyIn(x);
      if(!_VD.areEqual(x,y))
	{
	  std::cout<<"Failed\n";
	  return 0;
	}
      LinBox::commentator.progress();
    }
  LinBox::commentator.stop ("passed",  0, "test-solver");
 
  
  return 0;
}
	

