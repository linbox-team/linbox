/* -*- mode:C++ -*- */
/* File: test-rational-solver.C
 * Author: Zhendong Wan, modified by David Pritchard
 */

#include <linbox/field/ntl-ZZ.h>
#include <linbox/field/modular-int.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/diagonal.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/matrix/dense.h>
#include <linbox/algorithms/diophantine-solver.h>
#include <iostream>
#include <fstream>
#include <linbox/randiter/random-prime.h>

#include <linbox/field/unparametric.h>
#include <linbox/field/PID-integer.h>
#include <linbox/field/PID-double.h>
#include <linbox/field/ntl.h>
// #include <linbox/field/givaro.h>
// #include <linbox/field/lidia.h>

#include <linbox/field/archetype.h>
#include <linbox/field/givaro.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/integer.h>
#include <vector>

#include <linbox/field/archetype.h>
#include <linbox/vector/vector-domain.h>

using namespace std;
using namespace LinBox;

enum ArgumentType {
  TYPE_NONE, TYPE_INT, TYPE_INTEGER, TYPE_DOUBLE, TYPE_BOOL
};

struct Argument 
{
	char             c;
	char            *helpString;
	ArgumentType     type;
	void            *data;
};

/* Display a help message on command usage */

void printHelpMessage (const char *program, Argument *args) 
{
	int i, l;

	// Skip past libtool prefix in program name
	if (!strncmp (program, "lt-", strlen ("lt-")))
		program += strlen ("lt-");

	std::cout << "Usage: " << program << " [options] [<report file>]" << std::endl;
	std::cout << std::endl;
	std::cout << "Where [options] are the following:" << std::endl;

	for (i = 0; args[i].c != '\0'; i++) {
		switch (args[i].type) {
		case TYPE_NONE:
		  cout << "  -" << args[i].c << "        ";
		  break;		  
		case TYPE_INT:
		case TYPE_INTEGER:
		case TYPE_DOUBLE:
		  cout << "  -" << args[i].c << ' ' << args[i].c << "      ";
		  break;
		case TYPE_BOOL:
		  cout << "  -" << args[i].c << " {YN+-} ";
		  break;
		}
		std::cout << args[i].helpString;
		for (int j=strlen(args[i].helpString); j<54; j++) std::cout<<' ';
		std::cout << " (default ";
		switch (args[i].type) {
		case TYPE_NONE:
		  cout << "off";
		  break;		  
		case TYPE_INT:
		  cout << *(int *) args[i].data;
		  break;
		case TYPE_INTEGER:
		  cout << *(Integer *) args[i].data;
		  break;
		case TYPE_DOUBLE:
		  cout << *(double *) args[i].data;
		  break;
		case TYPE_BOOL:
		  cout << ((*(bool *)args[i].data)?"ON":"OFF");
		  break;
		}
		std::cout << ")" << std::endl;
	}

	std::cout << "  -h or -?  Display this message" << std::endl;
	std::cout << std::endl;
	std::cout << "If <report file> is not given, then no detailed reporting is done. This is" << std::endl;
	std::cout << "suitable if you wish only to determine whether the tests succeeded." << std::endl;
	std::cout << std::endl;
	std::cout << "[1] N.B. This program does not verify the primality of Q, and does not use a" << std::endl;
	std::cout << "    field extension in the event that Q=p^n, n > 1" << std::endl;
	std::cout << std::endl;
}

/* Find an argument in the argument list for a character */

Argument *findArgument (Argument *args, char c) 
{
	int i;

	for (i = 0; args[i].c != '\0' && args[i].c != c; i++);

	if (args[i].c != '\0')
		return &(args[i]);
	else
		return (Argument *) 0;
}

/* Parse command line arguments */

void parseArguments (int argc, char **argv, Argument *args)
{
	int i;
	Argument *current;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			if (argv[i][1] == 'h' || argv[i][1] == '?') {
				printHelpMessage (argv[0], args);
				exit (1);
			}
			else if ((current = findArgument (args, argv[i][1])) != (Argument *) 0) {
				switch (current->type) {
				case TYPE_NONE:
					*(bool *) current->data = true;
					break;

				case TYPE_INT:
					*(int *) current->data = atoi (argv[i+1]);
					i++;
					break;

				case TYPE_INTEGER:
				  *(Integer *) current->data = Integer(argv[i+1]);
					i++;
					break;

				case TYPE_DOUBLE:
					*(double *) current->data = atof (argv[i+1]);
					i++;
					break;

				case TYPE_BOOL:
				  if (argc == i+1 || strlen(argv[i+1]) > 1) {
				        *(bool *) current->data = true;
				        break;
				  }
					*(bool *) current->data = (argv[i+1][0] == '+' || argv[i+1][0] == 'Y' || argv[i+1][0] == 'y') ;
					i++;
					break;
				}
			} else {
				std::cerr << "ERROR: Bad argument " << argv[i] << std::endl;
				break;
			}
		} 
	}
}

/*********************************************************************************************/
//end of header stuff

//here Domain is a ring supporting the gcd, eg NTL_ZZ
template<class Domain>
bool reduceIn(Domain& D, std::pair<typename Domain::Element, 
	      typename Domain::Element> &frac)
{
  if (D.isZero(frac.second))
    {
      if (!D.isZero(frac.first)) D.init(frac.first, 1);
      // cerr << "division by zero has occured\n";
      return false;
    }
  if (D.isZero(frac.first))
    {
      D.init(frac.second, 1);
      return true;
    }
  typename Domain::Element gcd;
  D.gcd(gcd, frac.first, frac.second);
  D.divin(frac.first, gcd);
  D.divin(frac.second, gcd);
  return true;
};

template<class Domain> 
class VectorFraction
{
public:
  typedef typename Domain::Element Element;
  typedef typename std::pair<Element, Element> Fraction;
  typedef typename std::vector<Fraction> FVector;
  typedef typename Vector<Domain>::Dense Vector;

  //fields
  Vector numer;
  Element denom;

  //constructor
  VectorFraction(Domain& D, FVector &frac, bool reduceInPlace = true)
  {
    typename FVector::iterator i;

    if (reduceInPlace)
      {
	for (i=frac.begin(); i!=frac.end(); i++)
	  if (!reduceIn(D, *i)) cout << "error, div by zero\n";
      }
    
    D.init(denom, 1);
    for (i=frac.begin(); i!=frac.end(); i++)
      D.lcmin(denom, i->second);
    
    numer = Vector(frac.size());
    typename Vector::iterator j;

    for (i=frac.begin(), j=numer.begin(); i!=frac.end(); i++, j++)
      {
	D.mul(*j, denom, i->first);
	D.divin(*j, i->second);
      }
    
  }
};

#define random() ((double)rand() / ((double)(RAND_MAX)+1))


int n = 12;  
int defaultPrime = 101; 
int numPrimes = 5;
bool useRandom = true;
bool useDeterm = true;
bool printStuff = false;

int useFiles = false;
integer eBound = 10000; 
double singularProportion = 0.1;
bool inconsistent = false;

int useTimer = true;     
int entrySeed = 12345;  

bool testPidDouble = false;


static Argument args[] = {
  { 'n', "Set dimension of test matrix to nxn",              TYPE_INT,     &n },
  { 'q', "Solve first over the field Z/qZ",                  TYPE_INT,     &defaultPrime },
  { 'm', "Try solving with up to m primes",                  TYPE_INT,     &numPrimes },
  { 'r', "Set random solving on/off",                        TYPE_BOOL,    &useRandom },
  { 'd', "Set deterministic solving on/off",                 TYPE_BOOL,    &useDeterm },
  { 'p', "Print lots of detail?",                            TYPE_BOOL,    &printStuff },
  { 'f', "Read space-separated data from files td-{A, b}.txt?",   TYPE_BOOL,    &useFiles},
  { 'b', "(If f=OFF) Entry bound is (-b, b]",                     TYPE_INTEGER, &eBound},
  { 'z', "(If f=OFF) Make nullity roughly z*n, z in [0, 1]",      TYPE_DOUBLE,     &singularProportion},
  { 'i', "(If f=OFF) Force inconsistent system",                  TYPE_BOOL,    &inconsistent},
  { 't', "(If f=OFF) Randomize with timer?",                 TYPE_BOOL,    &useTimer},
  { 'w', "(If f=OFF, t=OFF) Randomize with seed w",          TYPE_INT,     &entrySeed},
  { 'e', "Test PID_double",          TYPE_BOOL,     &testPidDouble}
};

int trialCount=0;
integer* Aentries;
integer* bentries;

template <class Ring, class Field>
int test()
{
  trialCount++;
  if (trialCount == 1) {
    cout << "n: "<<n<<"\n";
    if (useTimer) 
      entrySeed = static_cast<unsigned>(time(NULL));
    if (!useFiles)
      cout << "Seed: " << entrySeed <<"\n";
    srand(entrySeed);
    
    Aentries = new integer[n*n];
    bentries = new integer[n];

    if (useFiles) {
      ifstream in;
      in.open("td-b.txt");
      for (int i=0; i<n; i++)
	in >> bentries[i];
      in.close();
      in.open("td-A.txt");
      for (int i=0; i<n*n; i++)
	in >> Aentries[i];
      in.close();
    }
    else {
      bool* auxRow = new bool[n];
      int auxRows = 0;
      for (int i=0; i<n; i++) {
	auxRow[i] = (random() <= singularProportion);
	if (auxRow[i]) auxRows++;
      }
      cout << "nullity >= " << auxRows << endl;

      if (inconsistent && auxRows == 0) 
	{ auxRows++; auxRow[(int)(random()*n)] = true; }

      if (auxRows > 0 && auxRows < n)
	eBound /= (n-auxRows);

      for (int i=0; i<n; i++) 	
	bentries[i] = random()*(double)(2*eBound+1)-(double)eBound;

      for (int i=0; i<n*n; i++) 
	Aentries[i] = random()*(double)(2*eBound+1)-(double)eBound;

      //make singular rows
      for (int i=0; i<n; i++)
	if (auxRow[i]) {
	  for (int j=0; j<n; j++) Aentries[n*i+j] = 0;
	  bentries[i] = 0;
	  for (int k=0; k<n; k++) 
	    if (!auxRow[k]) {
	      int m = (int)(random()*2)*2 - 1;
	      for (int j=0; j<n; j++)
		Aentries[n*i+j] += Aentries[n*k+j]*m;
	      bentries[i] += bentries[k]*m;
	    }
	  if (inconsistent) 
	    bentries[i] += (int)(random()*2)*2 - 1;
	}
    } 
  }

  typedef typename Ring::Element RingElement;
  
  Ring R;
  VectorDomain<Ring> VD (R);
  typedef Vector<Ring>::Dense Vector;
  Vector b(n);
  typedef DenseMatrix<Ring> Matrix; 
  Matrix A(R, n, n);
  MatrixDomain<Ring> MD(R);

  typename Vector::iterator bi=b.begin();
  for (int i=0; bi!=b.end(); bi++, i++)
    R.init(*bi, bentries[i]);
  if (trialCount==1 && printStuff) {cout << "b:\n"; VD.write(cout, b);}
	
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      R.init(A[i][j], Aentries[i*n+j]);
  if (trialCount==1 && printStuff) {cout << "\nA:\n"; A.write(cout);}

  Field F(2); //some givaro fields crash on small primes
  cout << "Testing with Z of type '";
  R.write(cout);
  cout<<"' and Z/pZ of type '";
  F.write(cout)<<"'"<<endl;

  typedef RationalSolver<Ring, Field, class RandomPrime, DixonTraits> QSolver; 


  typedef std::vector<std::pair<RingElement, RingElement> > FractionVector;
  FractionVector x(n);

  for (int i=0; i<2; i++)
    {
      if (i==0 && !useDeterm) continue;
      if (i==1 && !useRandom) continue;
      //clear x
      {
	FractionVector::iterator i;
	for (i=x.begin(); i!=x.end(); i++) {
	  R.init(i->first, 0);
	  R.init(i->second, 0);
	}
      }

      QSolver rsolver(defaultPrime, R, LinBox::RandomPrime(14)); //to avoid Log16 overflows
      DiophantineSolver<QSolver> zsolver(rsolver);
      SolverReturnStatus s;
      
      if (i==0) {
	cout << "Solving deterministically.\n";
	s = zsolver.solve(x, A, b, numPrimes);
      }
      else {
	cout << "Solving randomly.\n";
	s = zsolver.randomSolve(x, A, b, numPrimes);
      }
      cout << "solverReturnStatus: " << solverReturnString[(int)s] << "\n";

      if (s == OK)
	{
	  if (printStuff) {
	    /*    cout << "Solution:[";
	    {
	      FractionVector::iterator i;
	      for (i=x.begin(); i!=x.end(); i++)
		  cout << i->first << "/" << i->second << " ";
	      cout << "]\n";
	      }*/
	  }
	  
	  bool badSolution = false;
	  {
	    FractionVector::iterator i;
	    for (i=x.begin(); i!=x.end(); i++)
	      badSolution |= R.isZero(i->second);
	  }
	  
	  if (badSolution)
	    {
	      cout << "SOLUTION HAS ZERO DENOMINATOR.\n";
	      continue;
	    }
	  VectorFraction<Ring> red(R, x);
	  
	  if (printStuff)
	    {
	      cout << "Reduced solution: ";
	      VD.write(cout, red.numer);
	      cout << "/"<< red.denom << "\n";
	    }
	  Vector LHS(n), RHS(b);
	  // check that Ax = b, if it thought it was okay
	  MD.vectorMul(LHS, A, red.numer);
	  VD.mulin(RHS, red.denom);
	  /* if (printStuff)
	    { 
	      cout << "Ax * denom: ";
	      VD.write(cout, LHS) << "\n";
	  
	      cout << "b * denom: ";
	      VD.write(cout, RHS) << "\n";
	      } */
	  cout << "Ax=b : " << (VD.areEqual(LHS, RHS)?"Yes":"NO!") << "\n";
	  
	  // check certificate of minimality
	}
      else if (s==INCONSISTENT)
	{
	  // check certificate of inconsistency
	}
    }

};

template <class Field>
int fieldTest()
{
  return
    0
    +test<NTL_ZZ, Field>()   
    +test<PID_integer, Field>() 
    +(testPidDouble?test<PID_double, Field>():0) 
	// */
  ;
};

int main (int argc, char **argv)
{
  parseArguments (argc, argv, args);

  fieldTest<GivaroGfq>();  

  fieldTest<Modular<integer> >(); 
  fieldTest<GivaroZpz<Log16> >(); 
  fieldTest<NTL_zz_p>();          
  
  //fieldTest<GivaroMontg>(); 
  fieldTest<GivaroZpz<Std16> >(); 
  fieldTest<GivaroZpz<Std32> >(); 
  fieldTest<GivaroZpz<Std64> >(); 

  fieldTest<NTL_ZZ_p>(); 

  fieldTest<Modular<int> >();   
  fieldTest<Modular<double> >(); 
  

  // */
  // this takes a long time to compile with all fields 
  // so comment out unused ones when debugging
  return 0;
}

