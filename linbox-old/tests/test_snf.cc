#include "LinBox/snf_p.h"
#include "LinBox/field_padic.h"
#include <fstream>

using namespace LinBox;
void main()
{
  LinBox::SNF_P snfp(3, 2);
  
  std::ifstream file("matrix.txt");
  
  std::vector<int> result;
  
  Field_Padic field(3, 2);
    
  if(!file.fail())  
  {
    typedef Dense_Matrix<int> matrix;
    matrix A;
    
    int rows, cols;
    file>>rows;
    file>>cols;
    A.resize(rows);
    
    for(std::vector<std::vector<int> >::iterator pr=A.begin();pr!= A.end();++pr)
      {
	pr->resize(cols);
	for(matrix::pointer p =pr->begin(); p!=pr->end();++p)
	  {
	    file.ignore(1);
	    file>>*p;
	  }
      }
    snfp(result,A,field);

    for(int i=0; i<result.size();++i)
      std::cout<<result[i]<<" ";
  }
  
}
