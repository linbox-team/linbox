/*file companion.h
 *author: Zhendong Wan for LinBox group.
 */

#include <vector>
#include <algorithm>

namespace LinBox
{
  template<class elt>
    class ComMatrix
    {
    public:
      typedef std::vector<elt> Vector;
      
      ComMatrix(const Vector& v):_vector(v){}

      Vector& apply(Vector& y, const Vector& x) const
	{
	  if(x.size()!=coldim())
	    return y;
	  
	  y.resize(rowdim());
	  
	  elt last=x[x.size()-1];

	  y[0]=_vector[0]*last*(-1);
	  
	  if(rowdim()==1)
	    return y;

	  for(int i=1;i<x.size();++i)
	    y[i]=x[i-1]-_vector[i]*last;
	  
	  return y;
	}


      Vector& applyTranspose(Vector& y, const Vector& x) const
	{
	  if(x.size()!=rowdim())
	    return y;
	  
	  y.resize(coldim());
	  
	  std::copy(x.begin()+1,x.end(),y.begin());
	  
	  int last=x.size()-1;

	  Vector::const_iterator p1,p2;
	  
	  for(p1=x.begin(),p2=_vector.begin();p1!=x.end();++p1,++p2)
	    y[last]-=(*p1)*(*p2);
	  
	  return y;
	}

      size_t rowdim(void) const
	{return _vector.size();}

      size_t coldim(void) const
	{return _vector.size();}

    private:
      Vector _vector;
    };
}
	    
      
