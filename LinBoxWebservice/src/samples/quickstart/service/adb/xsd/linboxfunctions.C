// linboxfunctions.C
// Matthew Fendt, Summer '07 Summer Research
// A C++ file that executes LinBox code that will be made into a shared library
// so that the LinBox code can be executed in Java

// @@@@@@@@@@@@@ Currently supported methods:
// 1) determinant
// 2) rank
// 3) trace
// 4) valence
// 5) Smith normal form -- problem with LinBox right now?


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
//#include "linbox/pipes.h" // iostream form of det, rank, etc.
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/solutions/trace.h"
#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/util/error.h"
#include "linbox/util/matrix-stream.h"

static char dfile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Det_Response.txt";
static char rfile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Rank_Response.txt";
static char tfile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Trace_Response.txt";
static char vfile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Val_Response.txt";
static char snffile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/SNF_Response.txt";
//-----------------------------------------------------------------------------

bool det(std::istream& matrix_in, std::ostream& det_out)
{
  try {
    typedef LinBox::PID_integer Integers;		
    Integers ZZ;

    // LinBox::MatrixStream<Integers> ms(ZZ, matrix_in);    
    LinBox::DenseMatrix<Integers> A(ZZ);
       A.read(matrix_in);
    
    Integers::Element det_A;
    
    LinBox::det (det_A, A);
  
    ZZ.write(det_out, det_A) << std::endl;
  } 
  
  catch (...) { return false; }
  
  return true;
}

char* detFiles(char* mat)
{
  // std:ofstream output(dfile);
  ofstream output(dfile);

  string s(mat);
  istringstream iss(s);

  if (!iss || !output || !det(iss, output))
    return const_cast<char*>("Error in computing determinant");
  output.close();

  char line[5000];
  char* result;

  ifstream f(dfile);
  if (f.is_open())
    {
      f.getline(line, 5000);
      result = line;
      f.close();
      remove(dfile);
      return result;
    }
  else
    {
      remove(dfile);
      return const_cast<char*>("Error in opening file");
    }
}

bool rank(std::istream& matrix_in, std::ostream& rank_out)
{
  try {
    typedef LinBox::PID_integer Integers;		
    Integers ZZ;
    LinBox::MatrixStream<Integers> ms(ZZ, matrix_in);
    
    LinBox::DenseMatrix<Integers> A(ms);
    //A.read(matrix_in);
    
    unsigned long rank_A;

    LinBox::rank (rank_A, A);
    rank_out << rank_A << std::endl;  
  } 
  
  catch (...) { return false; }
  
  return true;

}

const char* rankFiles(char* mat)
{	
  // std:ofstream output(rfile);
  ofstream output(rfile);

  string s(mat);
  istringstream iss(s);

  if (!iss || !output || !rank(iss, output))
    return const_cast<char*>("Error in computing rank");
  output.close();

  char line[5000];
  char* result;

  ifstream f(rfile);
  if (f.is_open())
    {
      f.getline(line, 5000);
      result = line;
      f.close();
      remove(rfile);
      return result;
    }
  else
    {
      remove(rfile);
      return const_cast<char*>("Error in opening file");
    }
}

// @@@@@@@@@@ Just a stub for now.  The idea is that LinBox will be able to 
// give us updated time estimates as the computation is taking place.
int estimateRankTime(char* matfile)
{
  return rand();
}

bool val(std::istream& matrix_in, std::ostream& val_out)
{
  try {
    typedef LinBox::PID_integer Integers;		
    Integers ZZ;
    
    LinBox::MatrixStream<Integers> ms(ZZ, matrix_in);
    typedef LinBox::SparseMatrix<Integers> Blackbox;
    Blackbox A(ms);

    Integers::Element val_A;
   
    LinBox::valence(val_A, A);
    
    val_out << val_A << std::endl;
  } 

  catch (...) 
    { return false; }

  return true;
}



const char* valFiles(char* mat)
{	
 std:ofstream output(vfile);

  string s(mat);
  istringstream iss(s);

  if (!iss || !output || !val(iss, output))
    return const_cast<char*>("Error in computing valence");
  output.close();

  char line[5000];
  char* result;

  ifstream f(vfile);
  if (f.is_open())
    {
      f.getline(line, 5000);
      result = line;
      f.close();
      remove(vfile);
      return result;
    }
  else
    {
      remove(vfile);
      return const_cast<char*>("Error in opening file");
    }
}



bool trace(std::istream& matrix_in, std::ostream& trace_out)
{
  try {
    typedef LinBox::PID_integer Integers;		
    Integers ZZ;

    //LinBox::MatrixStream<Integers> ms(ZZ, matrix_in);    
        LinBox::DenseMatrix<Integers> A(ZZ);
        A.read(matrix_in);
    
    Integers::Element trace_A;
    
    
    LinBox::trace(trace_A, A);
    ZZ.write(trace_out, trace_A) << std::endl;
  } 

  catch (...) { return false; }

  return true;
}


const char* traceFiles(char* mat)
{	
 std:ofstream output(tfile);

  string s(mat);
  istringstream iss(s);

  if (!iss || !output || !trace(iss, output))
    return const_cast<char*>("Error in computing trace");
  output.close();

  char line[5000];
  char* result;

  ifstream f(tfile);
  if (f.is_open())
    {
      f.getline(line, 5000);
      result = line;
      f.close();
      remove(tfile);
      return result;
    }
  else
    {
      remove(tfile);
      return const_cast<char*>("Error in opening file");
    }
}



bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out)
{
  /*
  typedef LinBox::PID_integer Integers;
  //    typedef Integers::Element integer;
  Integers ZZ;
  typedef LinBox::DenseMatrix<Integers> Matrix;

  Matrix A(ZZ);
  A.read(matrix_in);

  vector<Integers::Element> ans(min(A.rowdim(), A.coldim()));
  //LinBox::SmithFormAdaptive::smithForm(ans, A);

  snf_out << "(";
  for (vector<Integers::Element>::iterator p = ans.begin();
       p != ans.end();
       ++p)
       ZZ.write(snf_out, *p) << " ";

  snf_out << ")" << endl;
  */

	snf_out << "SNF stub" << endl;

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


const char* smithNormalFormFiles(char* mat)
{	
 std:ofstream output(snffile);

  string s(mat);
  istringstream iss(s);

  if (!iss || !output || !smithNormalForm(iss, output))
    return const_cast<char*>("Error in snfFiles");
  output.close();

  char line[5000];
  char* result;

  ifstream f(snffile);
  if (f.is_open())
    {
      f.getline(line, 5000);
      result = line;
      f.close();
      remove(snffile);
      return result;
    }
  else
    {
      remove(snffile);
      return const_cast<char*>("Error in opening file");
    }
}
