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
  template <class element>
    class Dense_Matrix :public  std::vector<std::vector<element> >
    {
    public:
      typedef std::vector<element>::iterator pointer;
     
      std::vector<element>& row (int i) {return this->operator[](i);}
      
      size_t rowdim (void) const {return this->size();}
      
      size_t coldim (void) const 
	{
	  if(this->empty())
	    return 0; 
	  else 
	    return (*this)[0].size();}

      Dense_Matrix (){}

      //reserve space for a row*col dimension matrix.
      Dense_Matrix (int row, int col)
	{
	  this->resize(row);
	  for(std::vector<std::vector<element> >::iterator p=this->begin();p!=this->end();p++)
	    p->resize(col);
	}

      template<class Field>
      Dense_Matrix (std::ifstream& file, Field field =Field())
	{read(file, field);}
      //read from file row by row.each data separate by space.
      //the file's first row, it is the number of rows and columns.

      template<class Field>
      void read (std::ifstream& file, Field field =Field())
	{
	  int rows, cols;
	  file>>rows;
	  file>>cols;
	  this->resize(rows);

	  for(std::vector<std::vector<element> >::iterator pr=this->begin();pr!=this->end();++pr)
	    {
	      pr->resize(cols);
	      for(pointer p =pr->begin(); p!=pr->end();++p)
		{
		  file.ignore(1);
		  field.read(file, *p);
		}
	    }
	}

      template<class Field>
      std::vector<element>& apply(std::vector<element>& y, const std::vector<element>& x,Field field =Field()) const
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
     
      template<class Field> 
      std::vector<element>& applyTranspose(std::vector<element>& y, const std::vector<element>& x,Field field =Field()) const
	{if(x.size()!=rowdim())
	  return y;
	else
	  {y.resize()(coldim());
	  for(int i=0;i<coldim();i++)
	    {y[i]=0;
	    for(j=0;j<rowdim();j++)
	      field.axpyin(y[i],(*this)[j][i],x[j]);
	    }
	  }
	}
	  

      template<class Field>
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
