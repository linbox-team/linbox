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
//#include "linbox/pipes.h" // iostream form of det, rank, etc.
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/valence.h"
#include "linbox/solutions/trace.h"
#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/util/error.h"

static char dfile[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Det_Response.txt";
static char dfile2[100]= "/home/fendt/apache-tomcat-6.0.18/webapps/axis2/WEB-INF/services/Det_Response2.txt";
//-----------------------------------------------------------------------------

bool det(std::istream& matrix_in, std::ostream& det_out)
{
  typedef LinBox::PID_integer Integers;		
  Integers ZZ;
  
  LinBox::DenseMatrix<Integers> A(ZZ);
  A.read(matrix_in);
  
  Integers::Element det_A;
  try {
    LinBox::det (det_A, A);
  } catch (LinBox::LinboxError e) { return false; }
  
  ZZ.write(det_out, det_A) << std::endl;
  
  return true;
}

const char* detFiles(char* matfile)
{
  // Open the output stream
  ostringstream output;

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !det(iss, output))
    return const_cast<char*>("Error in detFiles");


  // FOR DEBUGGING
  ofstream output2(dfile);  
  ofstream output3(dfile2);
  output2 << "Part 1 " << output.str();
  output3 << "Part 2 " << output.str().c_str();
  output2.close();
  output3.close();
  

  return output.str().c_str();


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

const char* rankFiles(char* matfile)
{	
  ostringstream output;

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !rank(iss, output))
    return const_cast<char*>("Error in rankFiles");

  return output.str().c_str();
}

int estimateRankTime(char* matfile)
{
  return rand();
}

bool val(std::istream& matrix_in, std::ostream& val_out)
{
  typedef LinBox::PID_integer Integers;		
  Integers ZZ;
  
  LinBox::MatrixStream<Integers> ms(ZZ, matrix_in);
  typedef LinBox::SparseMatrix<Integers> Blackbox;
  Blackbox A(ms);

  Integers::Element val_A;
  LinBox::valence(val_A, A);

  val_out << val_A << std::endl;

  return true;
}



const char* valFiles(char* matfile)
{	
  ostringstream output;

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !val(iss, output))
    return const_cast<char*>("Error in valFiles");

  return output.str().c_str();

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


const char* traceFiles(char* matfile)
{	
  ostringstream output;

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !trace(iss, output))
    return const_cast<char*>("Error in traceFiles");

  return output.str().c_str();
}



bool smithNormalForm(std::istream& matrix_in, std::ostream& snf_out)
{
#if 0
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

#endif

	string s;
  matrix_in >> s;
    snf_out << "SNF stub: " << s <<  endl;

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


const char* smithNormalFormFiles(char* matfile)
{	
  ostringstream output;

  string s(matfile);
  istringstream iss(s);

  // If there is a problem, return false, otherwise return true
  if (!iss || !output || !smithNormalForm(iss, output))
    return const_cast<char*>("Error in snfFiles");

  return output.str().c_str();


    //  return const_cast<char*>("A test of SNF");





  /*
 std:ofstream output(sfile);

  string s(matfile);
  istringstream iss(s);

  if (!iss || !output || !smithNormalForm(iss, output))
    return const_cast<char*>("Error in snfFiles");

  output.close();

  char line[500000];
  char* result;

  ifstream f(sfile);
  if (f.is_open())
    {
      f.getline(line, 500000);
      result = line;
      f.close();
      remove(sfile);
      return result;
    }
  else
    {
      remove(sfile);
      return const_cast<char*>("Error in opening file");
    }
  */

}
