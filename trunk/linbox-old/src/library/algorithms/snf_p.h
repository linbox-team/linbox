#ifndef _SNF_P_
#define _SNF_P_

#include "LinBox/dense_matrix.h"
#include "LinBox/field_padic.h"

//#include <fstream>

namespace LinBox
{
  /// SNF_p is a function object to compute Smith forms over the integers.
  class SNF_p
    {
    public:

      /// Construct using the same parameters to the padic field constructor.
      SNF_p(int prime,int accurate):Pri(prime),Acc(accurate){}
     
      /** Smith Normal form for a dense integer matrix
       * @return  result	vector representing the diagonal of the Smith Form.
       * @param   matrix	dense matrix over the Field_Padic approximate representation
	*                       of the padic integers.
	*/
      std::vector<int>& operator()(std::vector<int>& result, Dense_Matrix<Field_Padic>& A)
	{
	  typedef  Field_Padic field_type;
	  typedef Dense_Matrix<field_type> matrix;
	  typedef matrix::pointer pointer;

	  Field_Padic field(Pri, Acc);
	  //matrix A;
	  //A.read(file, field);
	  //A.print(std::cout,field);

	  int size=A.row_dim();
	  if(A.col_dim()<A.row_dim())
	    size=A.col_dim();

	  result=std::vector<int> (size, Acc);
	  
	  for(int i=0; i<size;++i)
	    {

	      int row_index=i;
	      int col_index=i;
	      int exp=Acc;
	      int j;
	      pointer p;

	      for(j=i;j<A.row_dim();++j)
		for(p=A[j].begin()+i;p!=A[j].end();++p)
		  if((p->exp <exp) && (!field.isZero(*p)))
		  {
		    row_index=j;
		    col_index=p-A[j].begin();
		    exp=p->exp;
		  }	      	 
	      
	      Field_Padic::element x;

	      if(row_index!=i)//exchange ith row with row_indexth row.	       
		  std::swap(A[i],A[row_index]);
	
	      if(col_index!=i)//exchange ith column with col_indexth column.
		for(j=i;j<A.row_dim();++j)
		  std::swap(A[j][i],A[j][col_index]);

	      result[i]=exp;
	     
	      field.assign(x,A[i][i]);
	      if(field.isZero(x)|| (exp>=Acc))
		return result;
	      x.exp=0;
	      field.invin(x);
	   
	      for(p=A[i].begin()+i;p!=A[i].end();++p)		
		field.mulin(*p,x);


	      for(j=i+1;j<A.row_dim();++j)
		{
		  field.assign(x,A[j][i]);
		  field.dec_exp(x, exp);
		  field.negin(x);
		  if(x.exp<Acc)
		    {
		      pointer p1=A[i].begin()+i;
		      for(p=A[j].begin()+i; p!=A[j].end();++p)
			{field.axpyin(*p, x, *p1); ++p1;}
		    }
		}
	    } 
	  return result;
	}	  
    private:
      int Pri;
      int Acc;
    };
}

#endif
