// example code generator gen.h bds 2021Apr

/* 
These calls from solutions/ are handled (with and without M)
charpoly(f, A, M)
det(d, A, M)
isPositiveDefinite(A, M) // over Z only
minpoly(f, A, M)
rank(r, A, M)
smithForm(S, A, M) // over Z only.  Aspirational: Z/nZ, F[x], F[x]/<f>
trace(t, A)
valence(v, A, M)

These calls from solutions/ have yet to be handled 
rowEchelon(E, A, M)
rowEchelon(E, T, A, M)
colEchelon(E, A, M)
colEchelon(E, T, A, M)
getEntry(e, A, i, j)
hadamard-bound.h ??
isPositiveSemidefinite(A, M) // over Z only
solve(x, A, b, M)

These calls (and solutions) are aspirational
frobeniusForm(F, A, M) 
nullspace(N,A,M)
over ID:
solve(xNum,xDen,A,B,M) // A_mxn * Num_nxk = den*B_mxk
  aka ratSolve
signature (s, A, M) 
solve(X, A, B, M) // matrices

To add a solution call go to solution "settings" section.
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <time.h>
using namespace std;

// line containers //
// These containers are used so that program code lines can be created in 
// a convenient order, then produced in the proper order.
vector<string> header;
vector<string> usage;
vector<string> definitions;
vector<string> inputs;
vector<string> computations;
vector<string> outputs;
vector<string> tail;

string empty = "";
string indent = "   ";

// have(), alias(), computation() are defined first (to give the general idea)
// The other support functions are defined after main().

bool have(string& s) { return s != empty; }

// add a 'using command
void alias(string nickname, string name) { 
   if (have(nickname)) 
      header.push_back( (string)"using " + nickname + " = " + name + ";" );
}

// add a line for calling the solution
void computation(string line, string prefix = empty) { computations.push_back(prefix + line); }

// support function declarations //

// Make the generator's command line arglist into one string.
string argString (int ac, char* av[], string sep = "-") ;

// Add an #include line with an optional using command following it.
void include(string path, string nickname = empty, string name = empty) ; 

// define a variable called <name>  
void define(string name, string type, string cstorArgs = empty, string comment = empty) ; 

void input (string name, size_t i) ; // assume it can read itself for now.

/*
 * domainName empty: use cout.operator<<()
 * domainName "self": use obj.write()
 * domainName any other: use domain.write();
*/
void output (string label, string objName, string domainName = empty) ; 

// Produce sequence of lines from one of the vectors above.
template<class Container>
ostream& progLines(ostream& os, string indentation, Container& C) ;

void setUsage(size_t p, size_t e, string indent) ;

void camelCaps2SepWords(string & solutionFile, string solutionName, char sep = '-') ;

int main(int ac, char* av[]) {
   if (ac < 3) {
      cout << "Usage: " << av[0] << " prob rep [[m] [p [e]]]" << endl;
      cout << "<prob> is rank, det, charpoly, minpoly, smith-form, trace, or valence." << endl;
      cout << "<rep> is d (dense internal rep), s (sparse internal rep)," << endl; 
      cout << " or t (Toeplitz as an illustration of a pure blackbox)" << endl;
      cout << "Note: For s and d, both dense and sparse file reps are accepted." << endl;

      cout << "Optional m specifies a method (a for Auto, b for Blackbox, ...)." << endl;
      // Assuming sparse elim is covered by auto.
      // If not, add option for this second char.

      cout << "Computation is over ZZ unless prime p is given," << endl;
      cout << " in which case e is optional (default 1) and computation is over GF(p^e)." << endl;

      return 0;
   } 
   /** Establish top few lines **/
   time_t genTime; time (& genTime);
   // det-d-3-10.C built by: "gen det d 3 10" at <date>
   header.push_back( (string)"// " + argString(ac, av, "-") + " built by: \"" + av[0] + " " + argString(ac,av," ") + "\" at " + ctime(&genTime) );
   header.push_back( "// Copyright (C) Givaro Team. Your license is the GNU Lesser GPL" );

   include("fstream");
   include("linbox/linbox-config.h");

   /** av[1]: Establish solution file **/
   string solutionFile; 
   string solutionName = av[1]; // solution function name in camelCaps.
   camelCaps2SepWords(solutionFile, solutionName, '-');
   include((string)"linbox/solutions/" + solutionFile + ".h");

   /** av[2]: Establish primary matrix A **/
   // We do this after establishing the ring 

   /** av[3] non-numeric: Establish method **/
   string method = "M", method_type = "LinBox::Method::"; 
   int p_index = 4; // av[p_index], if it exists, is p.
   if (ac > 3) 
      switch (av[3][0]) {
      case ('a'): define(method, method_type + "Auto"); break;
      case ('b'): define(method, method_type + "Blackbox"); break;
      case ('e'): define(method, method_type + "Elimination"); break;
      case ('s'): define(method, method_type + "SparseElimination"); break;
      default:    method = empty; p_index = 3; 
      }
   else 
      method = empty;
//  include("linbox/solutions/methods.h"); // included with solutionFile

   /** the rest, av[3,4] or av[4,5]: Establish usage message, 
     * define p and e (if necessary), and establish the ground ring.  
     **/
   size_t p = 0, e = 1;
   if (ac > p_index) p = atoi(av[p_index]); // the characteristic
   if (ac > p_index+1) e = atoi(av[p_index+1]); // the exponent if nonprime Galois field.

   setUsage(p, e, indent);

   // Get p,e from run time command line args with default to the gen time args.
   if (p > 0) {
      define ("p", "LinBox::integer", 
         (string)"= argc > 2 ? atoi(argv[2]) : " + av[p_index]);
      if (e > 1) 
         define ("e", "size_t", 
            (string)"= argc > 3 ? atoi(argv[3]) : " + av[p_index+1]);
   }

   string ring, ringType, ringClass, ringCstorArgs, ringPath;
   if (e > 1) {
      ringCstorArgs = "p, e";
      cerr << av[0] << ": Non prime Galois Field not yet supported" << endl;
      return -1;
   } else if (p == 0) {
      ring = "ZZ"; 
      ringType = "IntDom"; 
      ringClass = "Givaro::ZRing<Givaro::Integer>";
      ringPath = "givaro/givintfactor.h"; 
      ringCstorArgs = empty;
   } else { // p > 0, e = 1
      ring = "F";
      ringType = "Field";
      ringCstorArgs = "p";
      ringPath = "givaro/modular.h"; 
      // potential special field rep for tiny primes
      // just illustrative for now
      if (p == 2) 
         ringClass = "Givaro::Modular<int16_t>";
      // potential selection of field rep by prime range
      // just illustrative for now
      size_t thresh1 = 4096, thresh2 = thresh1*thresh1;
      if (2 <= p  and p < thresh1) 
         ringClass = "Givaro::Modular<float>"; // 2^12
      if (thresh1 <= p  and p < thresh2) 
         ringClass = "Givaro::Modular<float,double>"; // 2^24
      size_t thresh3 = thresh2*1024;
      if (thresh2 <= p  and p < thresh3) 
         ringClass = "Givaro::Modular<uint64_t>"; // 2^32
      if (thresh3 <= p) {
         cerr << av[0] << ": characteristic > " << thresh3 << " not supported currently."; return -1;
      }
   } // end of else p > 0, e = 1

   include(ringPath, ringType, ringClass);
   define (ring, ringType, ringCstorArgs);


   /** Establish A (there is always input of a primary matrix named A) **/
   string matPath, matClass, matComment;
   switch (av[2][0]) {
   case 'd':
      matPath = "linbox/matrix/dense-matrix.h"; 
      matClass = "LinBox::DenseMatrix<" + ringType + ">";
      break;
   case 's':
      matPath = "linbox/matrix/sparse-matrix.h"; 
      matClass = "LinBox::SparseMatrix<" + ringType + ">";
      break;
   case 't':
      matComment = "//Using Toeplitz as illustrative blackbox type.\n";
      matPath = "linbox/blackbox/toeplitz.h";
      matClass = "LinBox::Toeplitz<" + ringType + ">";
      break;
   default: 
      cerr << "unknown matrix representation selector " << av[2] << endl; 
      return -1;
   }
   //include (matPath);
   //define ("A", matClass, ring, matComment);
   include (matPath, "Matrix", matClass);
   define ("A", "Matrix", ring, matComment);
   input("A", 1);

   /*** solution settings section ***/
   /**: establish other vars, inputs, solution call, outputs **/

   string LB = "LinBox::";// prefix for solution calls.
   // mTail is end of arg list with or without a method.
   string mTail = ", M);"; if (not have (method)) mTail = ");"; 

   if (solutionName == "charpoly" or solutionName == "minpoly") { 
      include ( "linbox/ring/polynomial-ring.h" );
      string polyRing = ring + 'x';
      string polyType = (string)"LinBox::DensePolynomial<" + ringType + ">";
      string polyRingType = polyType + "::Domain_t";
      //define(polyRing, polyRingType, ring + ", 'x'", "// just for io");
      // place this def in the output are to be clear it is only needed here.
      outputs.push_back(polyRingType + " " + polyRing + "(" + ring + ", 'x');");

      if (solutionName == "charpoly") { 
         define("c", polyType, ring);
         computation(LB + "charpoly(c, A" + mTail);
         output("charpoly = ", "c", polyRing);
      }
      if (solutionName == "minpoly") { 
         define("m", polyType, ring);
         computation(LB + "minpoly(m, A" + mTail);
         output("minpoly = ", "m", polyRing);
      }
   }
   if (solutionName == "det") {
      define("d", ringType + "::Element");
      computation(LB + "det(d, A" + mTail);
      output("det = ", "d", ring);
   }

   if (solutionName == "isPositiveDefinite") {
      if (p > 0) {
         cerr << "isPositiveDefinite applies to integer matrix only." << endl; 
         return -1; 
      }
      define("b", "bool");
      computation ((string)"b = " + LB + "isPositiveDefinite(A" + mTail);
      outputs.push_back("std::cout << \"isPositiveDefinite: \" << (b ? \"true\" : \"false\") << std::endl;");
   }

   if (solutionName == "rank") { 
   	//workaround for glitch in rank.h
      if (p == 0) solutionName = "integral_rank";

      define("r", "size_t");
      computation(LB + "rank(r, A" + mTail);
      output("rank = ", "r");
   }

   if (solutionName == "smithForm") { // smithForm(S, A, M)
      if (p > 0) {
         cerr << "Currently SmithForm applies to integer matrix only." << endl; 
         return -1; 
      }
   // illustrating that a solution/ might be skipped (discouraged)
   // (solutions/smith-form seems to have problems...)
      define("S", "LinBox::BlasVector<" + ringType + ">", ring);
      switch (av[2][0]) {
      case('d'):
         include("linbox/algorithms/smith-form-adaptive.h");
         computation(LB + "smithForm(S, A" + mTail);
      case('s'):
         include("linbox/algorithms/smith-form-valence.h");
         computation(LB + "smithForm(S, A" + mTail);
         break;
      default:
         break;
      }
      // does S need to be sized?
      computation(LB + "smithForm(S, A" + mTail);
      output("Smith Invariants = ", "S", "self"); 
      // change to compressed form...
   }

   if (solutionName == "trace") { 
      if (have(method)) {
         // in exception to solution convention about methods
         cerr << "trace accepts no method." << endl; 
         return -1; 
      }
      computation(LB + "trace(t, A);");
      define("t", ringType + "::Element");
      output("trace = ", "t", ring);
   }

   if (solutionName == "valence") { 
      define("v", ringType + "::Element");
      computation(LB + "valence(v, A" + mTail);
      output("valence = ", "v", ring);
   }

   /*** output the program ***/

   string programFile = argString(ac, av, "-") + ".C";
   ofstream pout(programFile);

   /** header **/
   progLines(pout, empty, header) << endl;

   pout << "int main(int argc, char* argv[]) {" << endl;
   progLines(pout, indent, usage) << endl;

   progLines(pout, indent, definitions) << endl;
   progLines(pout, indent, inputs) << endl;
   progLines(pout, indent, computations) << endl;
   progLines(pout, indent, outputs) << endl;
   progLines(pout, indent, tail);

   pout << "}" << endl;
   pout.close();
   
} // main

// Make the generator's command line arglist into one string
string argString (int ac, char* av[], string sep) {
// A typical sep is " ", "-", or "_".
   string ret = av[1];
   for (int i = 2; i < ac; ++i) 
      ret += sep + av[i];
   return ret;
}
// Add a #include line with optional using command following
void include(string path, string nickname, string name) { 
   if (not have(path)) return;
   header.push_back( (string)"#include <" + path + ">");
   if (have(nickname)) alias(nickname, name);
}

// define var called <name>  
void define(string name, string type, string cstorArgs, string comment) { 
   if (not have(name)) 
      return;
   if (not have(type))
      cerr << "Defining " << name << ", must have a type." << endl;

   if (have(comment)) definitions.push_back(comment);
   // two kinds constructor call: "T x(args)" and "T x = stuff"
   if (cstorArgs.size() >= 1 and cstorArgs[0] != '=')
      definitions.push_back(type + " " + name + "(" + cstorArgs + ");");
   else
      definitions.push_back(type + " " + name + cstorArgs + ";");
}

void input (string name, size_t i) { // assume it can read itself for now.
   if (not have(name)) return;
   stringstream s; s << "std::ifstream matStream(argv[" << i << "]);";
   inputs.push_back( s.str() );
   inputs.push_back( name + ".read(matStream);" );
   inputs.push_back( "matStream.close();" );
}

/*
 * domainName empty: use cout.operator<<()
 * domainName "self": use obj.write()
 * domainName any other: use domain.write();
*/
void output (string label, string objName, string domainName) { 
   if (not have(objName)) return;
   // label should end with space or newline
   string labeledOs = (string)"std::cout << \"" + label + '\"';
   if (domainName == "self") // object writes
      outputs.push_back( objName + ".write(" + labeledOs + ") << \"\\n\";" );
   else if (domainName != empty) // domain of object writes
      outputs.push_back( domainName + ".write(" + labeledOs + ", " + objName +") << \"\\n\";" );
   else // domainName = empty means object type is known to stdio
      outputs.push_back( labeledOs + " << " + objName + " << \"\\n\";" );
}

template<class Container>
ostream& progLines(ostream& os, string indentation, Container& C) {
   for (typename Container::iterator p  = C.begin(); p != C.end(); ++p) 
      os << indentation << *p << "\n";
   return os;
}

void camelCaps2SepWords(string & solutionFile, string solutionName, char sep) {
   solutionFile.resize(0);
   bool prevIsCap = false;
   for (int i = 0; i < solutionName.size(); ++i) {
      char x = solutionName[i];
      if (x < 'A' or 'Z' < x) { 
         solutionFile.push_back(x); 
         prevIsCap = false;
      } else {
         if (not prevIsCap) solutionFile.push_back(sep);
         solutionFile.push_back(x + 'a' - 'A'); 
         prevIsCap = true;
      }
   }
}

void setUsage(size_t p, size_t e, string indent) {
   int ac = 3; // ac is maximal runtime arg count 
   if (p == 0) ac = 2;
   if (p > 0 and e == 1) ac = 3;
   if (p > 0 and e > 1) ac = 4;
   string firstLine = "std::cout << \"Usage: \" << argv[0] << \" matfile ";
   switch (ac) { // set condition and prepare first line of usage statement
      case (2): usage.push_back("if (argc != 2) {"); 
                break;
      case (3): usage.push_back("if (argc < 2 or argc > 3) {"); 
                firstLine += "[p]";
                break;
      case (4): usage.push_back("if (argc < 2 or argc > 4) {");
                firstLine += "[p [e]]";
                break;
      default:  ;
   }
   usage.push_back(indent + firstLine + "\\n\";");

   usage.push_back( indent +  // 2nd line
      "std::cout << \"where matfile contains the matrix in any supported file format.\\n\";"); 
   if (ac > 2) {
      stringstream s; s << "std::cout << \"Optional p (default " << p << 
         ") is the characteristic in range [\" << Field::minCardinality() << \"..\" << Field::maxCardinality() << \").\\n\";";
      usage.push_back( indent + s.str() ); // 3rd line
   }
   if (ac > 3) {
      stringstream s; s << "std::cout << \"Optional e (default " << e << ") is the exponent for GF(p^e).\\n\";";
      usage.push_back( indent + s.str() ); // 4rd line
   }
   usage.push_back( indent + "return -1;" );
   usage.push_back("}");
}

/*
types: s, d, e, m, v, p, vp
ground domain Domain  
index size_t
Domain::Element
MatrixDomain ?
MatrixDomain::Element
Matrix< Ring >
VectorDomain ?
VectorDomain::Element Vector< Ring >
PolynomialDomain ?
PolynomialDomain::Element Polynomial< Ring >
Vector< Polynomial< Ring > >
types have name, path, cstorArgs
special includes: io, prob, ground ring (domain)
var domain includes: at most one per domain
var def, input, output: one per arg
only ground ring is typedef'ed to one of categories Ring/ID/PIR, 
domain_t domain_n domain_cstorArgs
for def
domain classes deliver: include, domdef,  from path, category, name
each domain element also delivers: eltdef, eltinput, eltoutput
size_t name = "size_t"; path = ""; 
def(var) = name + var + ['(' + cstorArgs + ')] +";\n"
input(var) = cin >> var;
output(var) = cout << "<var> = " << <var> << end; 
*/

