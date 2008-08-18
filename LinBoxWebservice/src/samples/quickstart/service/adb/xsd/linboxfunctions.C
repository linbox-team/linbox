// linboxfunctions.C
// Matthew Fendt, Summer '07 Summer Research
// A C++ file that executes LinBox code that will be made into a shared library
// so that the LinBox code can be executed in Java

// @@@@@@@@@@@@@ Currently supported methods:
// 1) determinant
// 2) rank
// 3) trace
// 4) valence
// 5) Smith normal form


#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <list>
#include <algorithm>
#include <vector>
#include <cstdlib>
using namespace std;

#include "linbox/blackbox/dense.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/solutions/trace.h"
#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/smith-form-adaptive.h"

// @@@@@@@@@@@@@@@@@@@@@@@@ Customize these variables

static char dfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Det_Response.txt";

static char rfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Rank_Response.txt";

static char vfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Val_Response.txt";

static char tfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/Trace_Response.txt";

static char sfile[100] = "/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/SNF_Response.txt";

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
  // Open the output stream
  std::ofstream output (dfile);

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !det(iss, output))
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
  // Open the output stream
  std::ofstream output (rfile);

  string s(matfile);
  istringstream iss(s);

  
  // If there is a problem, return false, otherwise return true
   if (!iss || !output || !rank(iss, output))
     return "Error in rankFiles";

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

double estimateRankTime(char* matfile)
{
  return 2.3;
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
  // Open the output stream
  std::ofstream output (vfile);

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !val(iss, output))
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
  // Open the output stream
  std::ofstream output (tfile);

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !trace(iss, output))
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



bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out)
{
  typedef LinBox::PID_integer Integers;
  //    typedef Integers::Element integer;
  Integers Z;
  typedef LinBox::DenseMatrix<Integers> Matrix;

  Matrix A(Z);
  A.read(matrix_in);

  vector<Integers::Element> ans(min(A.rowdim(), A.coldim()));
  LinBox::SmithFormAdaptive::smithForm(ans, A);

  snf_out << "(";
  for (vector<Integers::Element>::iterator p = ans.begin();
       p != ans.end();
       ++p)
    Z.write(snf_out, *p) << " ";

  snf_out << ")" << endl;







  //  snf_out << "snf test";




  /*
  typedef LinBox::PID_integer Integers;		
  typedef Integers::Element integer;
  Integers ZZ;

  LinBox::DenseMatrix<Integers> A(ZZ);
  A.read(matrix_in);
  
  typedef std::list<std::pair<integer, size_t> > SNF;
  SNF snf_A;
  LinBox::smithForm(snf_A, A);

  snf_out << "(";
  SNF::iterator p = snf_A.begin();

  if (snf_A.size() > 0) 
    ZZ.write(snf_out << ", [", p->first) << ", " << p->second << "]";
  for (++p; p != snf_A.end(); ++p)
    ZZ.write(snf_out << ", [", p->first) << ", " << p->second << "]";
  snf_out << ")" << endl;
  */
  return true;
}


char* smithNormalFormFiles(char* matfile)
{	
  // Open the output stream
  std::ofstream output (sfile);

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !smithNormalForm(iss, output))
    return "Error in smithNormalFormFiles";

  output.close();
  
  char line[MAX_ANSWER_CHARACTERS];
  char* result;
  
  ifstream f(sfile);
  if (f.is_open())
    {
      f.getline(line, MAX_ANSWER_CHARACTERS);
      result = line;
      f.close();
      remove(sfile);
      return result;
    }
  else
    {
      remove(sfile);
      return "Error in opening file";
    }
}
