#include "LinBox/snf_p.h"
#include "LinBox/field_padic.h"
#include <fstream>

void main()
{
  LinBox::SNF_p snfp(3, 2);
  std::ifstream file("matrix.txt");
  std::vector<int> a;
  Field_Padic field(3, 2);
  if(!file.fail())  
  {
     Dense_Matrix<Field_Padic> A;
     A.read(file, field);
    snfp(a,A);
  }
  for(int i=0; i<a.size();++i)
    std::cout<<a[i]<<" ";
}
