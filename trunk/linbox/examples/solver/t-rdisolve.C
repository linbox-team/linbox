/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* File: t-rdisolve.C - demo and testing of rational/diophantine system solver
 * Author: David Pritchard
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

#include <linbox/../tests/test-common.C>
using namespace std;
using namespace LinBox;

#define random() ((double)rand() / ((double)(RAND_MAX)+1))

int n = 30;  
int c = 30;
int defaultPrime = 101; 
int numPrimes = 5;
bool useRandom = false;
bool useDeterm = false;
bool useDiophantine = true;
bool useCertDio = true;
bool printStuff = true;

int useFiles = false;
integer eBoundCmd = 100; 
double singularProportion = 0;
bool inconsistent = false;

int useTimer = true;     
int entrySeed = 12345;  
int trials = 1;

bool testPidDouble = false;

static Argument args[] = {
	{ 'n', 0, "Row dimension of test matrix",                        TYPE_INT,     &n },
	{ 'c', 0, "Column dimension of test matrix (c<=0 => c=n)",       TYPE_INT,     &c },
	{ 'q', 0, "Solve first over the field Z/qZ",                     TYPE_INT,     &defaultPrime },
	{ 'm', 0, "Try solving with up to m primes",                     TYPE_INT,     &numPrimes },
	{ 'r', 0, "Set random solving on/off",                           TYPE_BOOL,    &useRandom },
	{ 'd', 0, "Set deterministic solving on/off",                    TYPE_BOOL,    &useDeterm },
	{ 'z', 0, "Set diophantine solving on/off",                      TYPE_BOOL,    &useDiophantine },
	{ 'v', 0, "Set certified diophantine solving on/off",            TYPE_BOOL,    &useCertDio},
	{ 'p', 0, "Print lots of detail?",                               TYPE_BOOL,    &printStuff },
	{ 'f', 0, "Read space-separated data from files td-{A, b}.txt?", TYPE_BOOL,    &useFiles},
	{ 'b', 0, "(If f=OFF) Entry bound is (-b, b]",                   TYPE_INTEGER, &eBoundCmd},
	{ 'x', 0, "(If f=OFF) Make roughly x*n dependant rows",          TYPE_DOUBLE,  &singularProportion},
	{ 'i', 0, "(If f=OFF) Force inconsistent system",                TYPE_BOOL,    &inconsistent},
	{ 't', 0, "(If f=OFF) Randomize with timer?",                    TYPE_BOOL,    &useTimer},
	{ 'w', 0, "(If f=OFF, t=OFF) Randomize with seed w",             TYPE_INT,     &entrySeed},
	{ 'e', 0, "Test PID_double",                                     TYPE_BOOL,    &testPidDouble},
	{ 'g', 0, "Repeat trials g times",                               TYPE_INT,     &trials}
};

int trialCount=0;
integer* Aentries;
integer* bentries;

template <class Ring, class Field>
int test() {
	trialCount++;

	typedef typename Ring::Element RingElement;
  
	Ring R;
	VectorDomain<Ring> VD (R);
	typedef Vector<Ring>::Dense Vector;
	Vector b(n);
	typedef DenseMatrix<Ring> Matrix; 
	Matrix A(R, n, c);
	MatrixDomain<Ring> MD(R);

	typename Vector::iterator bi=b.begin();
	for (int i=0; bi!=b.end(); bi++, i++)
		R.init(*bi, bentries[i]);
	if (trialCount==1 && printStuff) {cout << "b:\n"; VD.write(cout, b);}
	
	for (int i=0; i<n; i++)
		for (int j=0; j<c; j++)
			R.init(A[i][j], Aentries[i*c+j]);
	if (trialCount==1 && printStuff) {cout << "\nA:\n"; A.write(cout);}

	Field F(2); //some givaro fields crash on small primes
	cout << "Testing with Z of type '";
	R.write(cout);
	cout<<"' and Z/pZ of type '";
	F.write(cout)<<"'"<<endl;

	typedef RationalSolver<Ring, Field, class RandomPrime, DixonTraits> QSolver; 
	typedef DiophantineSolver<QSolver> ZSolver; 

	typedef std::vector<std::pair<RingElement, RingElement> > FractionVector;
	FractionVector x(n);
	int result=0;

	for (int iteration=0; iteration<4; iteration++) {
		if (iteration==0 && !useDeterm) continue;
		if (iteration==1 && !useRandom) continue;
		if (iteration==2 && !useDiophantine) continue;
		if (iteration==3 && !useCertDio) continue;

		//clear x				
		for (FractionVector::iterator i=x.begin(); i!=x.end(); i++) {
			R.init(i->first, 0);
			R.init(i->second, 0);
		}

		QSolver rsolver(defaultPrime, R, LinBox::RandomPrime(14)); //to avoid Log16 overflows
		ZSolver zsolver(rsolver);
		SolverReturnStatus s;
      
		if (iteration==0) {
			cout << "Solving deterministically.\n";
			s = zsolver.solve(x, A, b, numPrimes);
		}
		else if (iteration==1) {
			cout << "Solving randomly.\n";
			s = zsolver.randomSolve(x, A, b, numPrimes);
		}
		else if (iteration==2) {
			cout << "Solving diophantically.\n";
			s = zsolver.diophantineSolve(x, A, b, false, numPrimes);
		}
		else {
			cout << "Solving diophantically with certification of min-denominator.\n";
			s = zsolver.diophantineSolve(x, A, b, true, numPrimes);
		}
		cout << "solverReturnStatus: " << solverReturnString[(int)s] << "\n";

		if (s == SS_OK)	{
			if (printStuff) {
// 				cout << "Solution:[";
// 				for ( FractionVector::iterator i=x.begin(); i!=x.end(); i++)
// 					cout << i->first << "/" << i->second << " ";
// 				cout << "]\n";
			}
	  
			bool errorInSolution = false;
			VectorFraction<Ring> red(R, x, errorInSolution);
			if (errorInSolution) {
				cout << "SOLUTION HAS ZERO DENOMINATOR.\n";
				continue;
			}
	  
			if (printStuff) {
				cout << "Reduced solution: ";
				red.write(cout) << "\n";
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
			
			if (iteration==3) {
				// check certificate of minimality z
				// should satisfy that zA is integral, and den(z.b) == den(y)
				
				typedef typename Ring::Element Integer;
				Integer dp, tmp, denzb;
				VectorFraction<Ring> z(zsolver.lastCertificate);
				R.init(dp, 0);
				typename Vector::iterator zi = z.numer.begin();
				typename Vector::iterator bi = b.begin();
				for (; bi != b.end(); bi++, zi++)
					R.addin(dp, R.mul(tmp, *bi, *zi));
				
				R.gcd(denzb, dp, z.denom);
				R.div(denzb, z.denom, denzb);
				
				bool error;
				VectorFraction<Ring> tmpvf(R, x, error);
				bool certified = (!error) && R.areEqual(denzb, tmpvf.denom);
				if (!certified)
					cout << "ERROR Failed den(z.b) == den(y)" << endl;

				bool certified2 = true;
				Integer* nza = new Integer[c]; //z.numer * A
				for (int i=0; i<c; i++) R.init(nza[i], 0);
				for (int i=0; i<n; i++)
					for (int j=0; j<c; j++)
						R.addin(nza[j], R.mul(tmp, z.numer[i], A[i][j]));
				for (int i=0; i<c; i++)
					certified2 &= R.isDivisor(nza[i], z.denom);
				if (!certified2)
					cout << "ERROR Failed zA integral" << endl;

				if (certified && certified2) 
					cout << "Solution is certified correctly." << endl;
			}
		}
		else if (s==SS_INCONSISTENT) {
			// check certificate of inconsistency
		}
	}
	return result;

};

template <class Field>
int fieldTest()
{
	return
		0
		+test<NTL_ZZ, Field>()   
		//+test<PID_integer, Field>() 
		+(testPidDouble?test<PID_double, Field>():0) 
		// */
		;
};

void genTestData() {
	bool* auxRow = new bool[n];
	int auxRows = 0;
	for (int i=0; i<n; i++) {
		auxRow[i] = (random() <= singularProportion);
		if (auxRow[i]) auxRows++;
	}
	cout << "at least " << auxRows << " dependent rows" << endl;

	if (inconsistent && auxRows == 0) 
		{ auxRows++; auxRow[(int)(random()*n)] = true; }

	integer eBound(eBoundCmd);
	if (auxRows > 0 && auxRows < n)
		eBound /= (n-auxRows);
	if (eBound == 0) {
		cout << "WARNING; 'b' dropped to 0. Changed to 1, try increasing 'b'." << endl;
		eBound = 1;
	}

	for (int i=0; i<n; i++) 	
		bentries[i] = random()*(double)(2*eBound+1)-(double)eBound;

	for (int i=0; i<n*c; i++) 
		Aentries[i] = random()*(double)(2*eBound+1)-(double)eBound;

	//make singular rows
	for (int i=0; i<n; i++)
		if (auxRow[i]) {
			for (int j=0; j<c; j++) Aentries[c*i+j] = 0;
			bentries[i] = 0;
			for (int k=0; k<n; k++) 
				if (!auxRow[k]) {
					int m = (int)(random()*2)*2 - 1;
					for (int j=0; j<c; j++)
						Aentries[c*i+j] += Aentries[c*k+j]*m;
					bentries[i] += bentries[k]*m;
				}
			if (inconsistent) 
				bentries[i] += (int)(random()*2)*2 - 1;
		}
}

int main (int argc, char **argv)
{
	parseArguments (argc, argv, args, true);

	if (useTimer) {
		entrySeed = static_cast<unsigned>(time(NULL));
		useTimer = false;
	}
	
	writeCommandString (cout, args, argv[0]);

	if (c <= 0) c += n;
	if (c <= 0) {
		cout << "WARNING, c <= -n; resulting column dimension changed from nonpositive value to 1" << endl;
		c = 1;
	}
	cout << "Matrix is n="<<n<<" rows by c="<<c<<" columns\n";

	cout << "Seed: " << entrySeed <<"\n";
	srand(entrySeed);
    
	Aentries = new integer[n*c];
	bentries = new integer[n];

	if (useFiles) {
		ifstream in, in2;
		in.open("td-b.txt");
		for (int i=0; i<n; i++)
			in >> bentries[i];
		in.close();
		in2.open("td-A.txt");
		for (int i=0; i<n*c; i++)
			in2 >> Aentries[i];
		in2.close();
	}

	for (int j=0; j < trials; j++) {
		if (!useFiles) genTestData();
		fieldTest<Modular<integer> >(); 
		fieldTest<GivaroZpz<Log16> >(); 
		fieldTest<NTL_zz_p>();          
		/*
		  fieldTest<GivaroGfq>();  
		  fieldTest<GivaroZpz<Std16> >(); 
		  fieldTest<GivaroZpz<Std32> >(); 
		  fieldTest<GivaroZpz<Std64> >(); 
		  
		  //fieldTest<GivaroMontg>();    -- appears to be broken
		  //fieldTest<NTL_ZZ_p>();       -- appears to be broken
		  
		  fieldTest<Modular<int> >();   
		  fieldTest<Modular<double> >(); 

	  // */
		// this takes a long time to compile with all fields 
		// so comment out unused ones when debugging
		cout << "finished trial " << (j+1) << " of " << trials << endl;
	}
	return 0;
}

// TODO: fix norm problems, add max-norm
// TODO: come up with better test data, so can have a big singular matrix of all 0..9

//more space; ./t-rdisolve -p n -n 10 -c 11 -x 0.05 -b 10 -g 30 -m 2 -q 7 > toto; more space; egrep -i "warn|erro|dam|no\!|fail" toto; more space; tail toto
