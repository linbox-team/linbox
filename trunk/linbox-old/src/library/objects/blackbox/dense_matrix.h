#ifndef _DENSE_MATRIX_
#define _DENSE_MATRIX_

#include <iostream>
#include <vector>
#include <fstream>

namespace LinBox
{
  /** Dense matrix class extends vector of rows representation, where a row is a vector of 
element.  Provides blackbox members and read() and write().
  */
// FIXME add the blackbox interface
  template <class Field>
    class Dense_Matrix :public  std::vector<std::vector<typename Field::element> >
    {
    public:
      typedef typename Field::element element;
      typedef std::vector<element>::iterator pointer;
     
      std::vector<element>& row (int i) {return this->operator[](i);}
      
      size_t rowdim (void) {return this->size();}
      
      size_t coldim (void) 
	{
	  if(this->empty())
	    return 0; 
	  else 
	    return (*this)[0].size();}
      Dense_Matrix (){}
      
      Dense_Matrix (std::ifstream& file, Field field =Field())
	{read(file, field);}
      //read from file row by row.each data separate by space.
      //the file's first row, it is the number of rows and columns.
      void read (std::ifstream& file, Field field =Field())
	{
	  int rows, cols;
	  file>>rows;
	  file>>cols;
	  this->resize(rows);

	  for(int i=0; i<rows; ++i)
	    {
	      (*this)[i].resize(cols);
	      for(pointer p =(*this)[i].begin(); p!=(*this)[i].end();++p)
		{
		  file.ignore(1);
		  field.read(file, *p);}}
	}

      std::vector<element>& apply(std::vector<element>& y, const std::vector<element>& x,field =Field())
	{if(x.size()!=coldim())
	  return y;
	else
	  {y.resize(rowdim());
	  for(int i=0;i<rowdim();i++)
	    {
	      y[i]=0; int j=0;
	      for(pointer p=(*this)[i].begin();p!=(*this)[i].end();p++)
		{
		  field.axpyin(y[i],*p,x[j]);
		  j++;
		}
	    }
	  }
	return y;}

      std::ostream& write(std::ostream& os =std::cout,Field field =Field())
	{
	  for(int i=0;i<this->rowdim();++i)
	    {
	      for(pointer p=(*this)[i].begin();p!=(*this)[i].end();++p)
		{field.write(os,*p);
		os<<" ";}
	      os<<"\n";
	    }	
	}
     

    };
}

#endif
