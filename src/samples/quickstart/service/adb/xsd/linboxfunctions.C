// linboxfunctions.C
// Matthew Fendt, Summer '07 Summer Research
// A C++ file that executes LinBox code that will be made into a shared library
// so that the LinBox code can be executed in Java

// @@@@@@@@@@@@@ Currently supported methods:
// 1) determinant (of int size)
// 2) rank
// 3) trace

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
using namespace std;

#include "linbox/blackbox/dense.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/solutions/trace.h"



// @@@@@@@@@@@@@@@@@@@@@@@@ Customize these variables

static char dfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Det_Response.txt";

static char rfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Rank_Response.txt";

static char vfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Val_Response.txt";

static char tfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Trace_Response.txt";

static int MAX_ANSWER_CHARACTERS = 500000;
//-----------------------------------------------------------------------------

bool det(std::istream& matrix_in, std::ostream& det_out)
{
	typedef LinBox::PID_integer Integers;		
	Integers ZZ;
	
	LinBox::DenseMatrix<Integers> A(ZZ);
	A.read(matrix_in);
	
	Integers::Element det_A;
	LinBox::det (det_A, A);
	
	ZZ.write(det_out, det_A) << std::endl;
	
	return true;
}

char* detFiles(char* matfile)
{
  // Open the input and output streams
  std::ifstream input (matfile);
  std::ofstream output (dfile);

  
  // If there is a problem, return false, otherwise return true
  if (!input || !output || !det(input, output))
    return "Error in detFiles";

  output.close();
  
  char line[MAX_ANSWER_CHARACTERS];
  char* result;
  
  ifstream f(dfile);
  
  if (f.is_open())
    {
      f.getline(line, MAX_ANSWER_CHARACTERS);
      result = line;
      f.close();
      
      remove(dfile);
      return result;
    }
	else
	  {
	    remove(dfile);
	    return "Error in opening file";
	  }
}

bool rank(std::istream& matrix_in, std::ostream& rank_out)
{
  typedef LinBox::PID_integer Integers;		
  Integers ZZ;
  
  LinBox::DenseMatrix<Integers> A(ZZ);
  A.read(matrix_in);
  
  unsigned long rank_A;
  LinBox::rank (rank_A, A);
  
  rank_out << rank_A << std::endl;
  
  return true;
}

char* rankFiles(char* matfile)
{	
  // Open the input and output streams
  std::ifstream input (matfile);
  std::ofstream output (rfile);
  
  // If there is a problem, return false, otherwise return true
  // if (!input || !output || !rank(input, output))
  // return "Error in rankFiles";
  if (!input)
    return "Error in rankFiles input";

  if (!output)
    return "Error in rankFiles output";

  if (!(rank(input,output)))
      return "Error in rankFiles rank";


  output.close();
  
  char line[MAX_ANSWER_CHARACTERS];
  char* result;

  ifstream f(rfile);
  if (f.is_open())
    {
      f.getline(line, MAX_ANSWER_CHARACTERS);
      result = line;
      f.close();
      remove(rfile);
      return result;
    }
  else
    {
      remove(rfile);
      return "Error in opening file";
    }
}



bool val(std::istream& matrix_in, std::ostream& val_out)
{
  typedef LinBox::PID_integer Integers;
  Integers ZZ;
  
  LinBox::MatrixStream< Integers > ms(ZZ, matrix_in);
  typedef LinBox::SparseMatrix<Integers> Blackbox;
  Blackbox A(ms);
  
  Integers::Element val_A;
  LinBox::valence(val_A, A);
  
  val_out << val_A << std::endl;
}



char* valFiles(char* matfile)
{	
  // Open the input and output streams
  std::ifstream input (matfile);
  std::ofstream output (vfile);
  
  // If there is a problem, return false, otherwise return true
  if (!input || !output || !val(input, output))
    return "Error in valFiles";
  
  output.close();
  
  char line[MAX_ANSWER_CHARACTERS];
  char* result;
  
  ifstream f(vfile);
  if (f.is_open())
    {
      f.getline(line, MAX_ANSWER_CHARACTERS);
      result = line;
      f.close();
      remove(vfile);
      return result;
    }
  else
    {
      remove(vfile);
      return "Error in opening file";
    }
}




bool trace(std::istream& matrix_in, std::ostream& trace_out)
{
  typedef LinBox::PID_integer Integers;		
  Integers ZZ;

  LinBox::DenseMatrix<Integers> A(ZZ);
  A.read(matrix_in);
  
  Integers::Element trace_A;
  LinBox::trace(trace_A, A);

  ZZ.write(trace_out, trace_A) << std::endl;

  return true;
}


char* traceFiles(char* matfile)
{	
  // Open the input and output streams
  std::ifstream input (matfile);
  std::ofstream output (tfile);
  
  // If there is a problem, return false, otherwise return true
  if (!input || !output || !trace(input, output))
    return "Error in traceFiles";

  output.close();
  
  char line[MAX_ANSWER_CHARACTERS];
  char* result;
  
  ifstream f(tfile);
  if (f.is_open())
    {
      f.getline(line, MAX_ANSWER_CHARACTERS);
      result = line;
      f.close();
      remove(tfile);
      return result;
    }
  else
    {
      remove(tfile);
      return "Error in opening file";
    }
}
