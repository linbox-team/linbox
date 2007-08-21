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

// The directory of the matrix file
static char i[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/";

// The directory where the answer will be generated
static char o[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/";
	
static char* dfile = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Det_Response.txt";

static char* rfile = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Rank_Response.txt";

static char* vfile = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Val_Response.txt";

static char* tfile = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Trace_Response.txt";

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

long long detFiles(char* matfile)
{
	char* outfile = "Det_Response.txt";

	// Add on the inputs
	strcat(i, matfile);
	strcat(o, outfile);
	
	// Open the input and output streams
	std::ifstream input (i);
	std::ofstream output (o);


	// If there is a problem, return false, otherwise return true
	if (!input || !output || !det(input, output))
		return -666;

	output.close();

	char line[200];
	long long result;

	ifstream f(dfile);

	if (f.is_open())
		{
		  f.getline(line, 200);
		  result = atoll(line);
		  f.close();
		  
		  remove(dfile);
		  return result;
		}
	else
		{
		  remove(dfile);
		  return -666;
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

int rankFiles(char* matfile)
{	
	char* outfile = "Rank_Response.txt";

	// Add on the inputs
	strcat(i, matfile);
	strcat(o, outfile);
	
	// Open the input and output streams
	std::ifstream input (i);
	std::ofstream output (o);

	// If there is a problem, return false, otherwise return true
	if (!input || !output || !rank(input, output))
		return -666;

	output.close();

	char line[200];
	int result;

	ifstream f(rfile);
	if (f.is_open())
		{
		  f.getline(line, 200);
		  result = atoi(line);
		  f.close();
		  remove(rfile);
		  return result;
		}
	else
		{
		  remove(rfile);
		  return -666;
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



int valFiles(char* matfile)
{	
  char* outfile = "Val_Response.txt";
  
  // Add on the inputs
  strcat(i, matfile);
  strcat(o, outfile);
  
  // Open the input and output streams
  std::ifstream input (i);
  std::ofstream output (o);
  
  // If there is a problem, return false, otherwise return true
  if (!input || !output || !val(input, output))
    return -666;
  
  output.close();
  
  char line[200];
  int result;
  
  ifstream f(vfile);
  if (f.is_open())
    {
      f.getline(line, 200);
      result = atoi(line);
      f.close();
      remove(vfile);
      return result;
    }
  else
    {
      remove(vfile);
      return -666;
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


int traceFiles(char* matfile)
{	
  char* outfile = "Trace_Response.txt";
  
  // Add on the inputs
  strcat(i, matfile);
  strcat(o, outfile);
  
  // Open the input and output streams
  std::ifstream input (i);
  std::ofstream output (o);
  
  // If there is a problem, return false, otherwise return true
  if (!input || !output || !trace(input, output))
    return -666;
  
  output.close();
  
  char line[200];
  int result;
  
  ifstream f(tfile);
  if (f.is_open())
    {
      f.getline(line, 200);
      result = atoi(line);
      f.close();
      remove(tfile);
      return result;
    }
  else
    {
      remove(tfile);
      return -666;
    }
}
