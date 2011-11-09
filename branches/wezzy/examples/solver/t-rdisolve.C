/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * examples/solver/t-rdisolve.C
 *
 * Copyright (C) 2004, 2005, 2010  D. Pritchard, P. Giorgi
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see
 *   <http://www.gnu.org/licenses/>.
 */

/* linbox/examples/solver/t-rdisolve.C
 * demo, testing, time-comparison of certified rational/diophantine system solver
 *
 * Written by David Pritchard  <daveagp@mit.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#define LIFTING_PROGRESS
#define RSTIMING


#include <linbox/integer.h>
#include <linbox/field/ntl-ZZ.h>
#include <linbox/field/modular-int.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/diagonal.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/algorithms/vector-fraction.h>
#include <linbox/matrix/dense.h>
#include <linbox/algorithms/diophantine-solver.h>
#include <iostream>
#include <fstream>
#include <linbox/randiter/random-prime.h>

#include <linbox/field/unparametric.h>
#include <linbox/field/PID-integer.h>
#include <linbox/field/PID-double.h>
#include <linbox/field/ntl.h>

#include <linbox/field/archetype.h>
//#include <linbox/field/givaro.h>
#include <linbox/vector/vector-domain.h>

#include <vector>

#include <linbox/field/archetype.h>
#include <linbox/vector/vector-domain.h>

#include <linbox/../tests/test-common.C>
using namespace std;
using namespace LinBox;

#define random_01() ((double)rand() / ((double)(RAND_MAX)+1))

int  n               = 5;
int  c               = 5;
int  defaultPrime    = 0;
int  primeBits       = 14;         // note: should be <= 15 to use GivaroZpz<Log16>
int  numPrimes       = 1;
bool useDeterm       = true;
bool useRandom       = false;
bool useDiophantine  = false;
int  printStuff      = 0;
int  showTiming      = 0;

bool    useFiles           = false;
bool    sparseMatrix       = false;
integer eBoundCmd          = 1000;
double  singularProportion = 0;
bool    inconsistent       = false;

int useTimer  = true;
int entrySeed = 12345;
int trials    = 1;

int destroyColumns = 0;

bool testPidDouble = false;

int levelAsInt = (int)SL_CERTIFIED;

static Argument args[] = {
	{ 'n', 0, "Row dimension of test matrix",                        TYPE_INT,     &n },
	{ 'c', 0, "Column dimension of test matrix (c<=0 => c=n)",       TYPE_INT,     &c },
	{ 'm', 0, "Try solving with up to m primes",                     TYPE_INT,     &numPrimes },
	{ 'q', 0, "Solve first over the field Z/qZ (0: pick randomly)",  TYPE_INT,     &defaultPrime },
	{ 'g', 0, "Subsequently generate primes that are g bits long",   TYPE_INT,     &primeBits },
	{ 'r', 0, "Set random solving on/off",                           TYPE_BOOL,    &useRandom },
	{ 'd', 0, "Set deterministic solving on/off",                    TYPE_BOOL,    &useDeterm },
	{ 'z', 0, "Set diophantine solving on/off",                      TYPE_BOOL,    &useDiophantine },
	{ 'p', 0, "Print lots of detail, tree levels (0,1,2,3)",           TYPE_INT,     &printStuff },
	{ 'f', 0, "Read space-separated data from files td-{A, b}.txt?", TYPE_BOOL,    &useFiles},
	{ 's', 0, "(If f=0N)  Say td-A.txt is in sparse format",         TYPE_BOOL,    &sparseMatrix},
	{ 'b', 0, "(If f=OFF) Entry bound is (-b, b]",                   TYPE_INTEGER, &eBoundCmd},
	{ 'x', 0, "(If f=OFF) Make roughly x*n dependant rows",          TYPE_DOUBLE,  &singularProportion},
	{ 'i', 0, "(If f=OFF) Force inconsistent system",                TYPE_BOOL,    &inconsistent},
	{ 't', 0, "(If f=OFF) Randomize with timer?",                    TYPE_BOOL,    &useTimer},
	{ 'w', 0, "(If f=OFF, t=OFF) Randomize with seed w",             TYPE_INT,     &entrySeed},
	{ 'e', 0, "Test PID_double",                                     TYPE_BOOL,    &testPidDouble},
	{ 'k', 0, "Repeat trials k times",                               TYPE_INT,     &trials},
	{ 'l', 0, "Level: 0=Monte Carlo, 1=Las Vegas, 2=Certified",      TYPE_INT,     &levelAsInt},
	{ 'o', 0, "Set o columns to zero at random",                     TYPE_INT,     &destroyColumns}
}; // 7 more options (yu vs jah) and the whole alphabet is covered

int trialCount=0;
integer* Aentries;
integer* bentries;

template <class Ring, class Field>
int test()
{
	trialCount++;

	typedef typename Ring::Element RingElement;

	Ring R;
	VectorDomain<Ring> VD (R);
	typedef typename Vector<Ring>::Dense Vector;
	typedef DenseMatrix<Ring> Matrix;
	Matrix A(R, n, c);
	MatrixDomain<Ring> MD(R);
	typedef typename Ring::Element Integer;
	Vector b(n);

	if (sparseMatrix) {
		// reading A from td-A.txt file
		ifstream inA, inb;
		inA.open("td-A.txt");
		A.read(inA);
		cout << "Matrix is sparse with n="<<A.rowdim()<<" rows by c="<<A.coldim()<<" columns\n";
		n = (int) A.rowdim();
		c = (int) A.coldim();
		inA.close();
		// reading b from td-b.txt
		b.resize(n);
		inb.open("td-b.txt");
		for (int i=0; i<n; i++)
			inb >> b[i];
		inb.close();
	}
	else {
		// reading A from Aentrie vector
		for (int i=0; i<n; i++)
			for (int j=0; j<c; j++)
				R.init(A[i][j], Aentries[i*c+j]);
		// reading b from bentry vector
		typename Vector::iterator bi=b.begin();
		for (int i=0; bi!=b.end(); bi++, i++)
			R.init(*bi, bentries[i]);
	}


	if (trialCount==1 && (printStuff>2)) {cout << "b:\n"; VD.write(cout, b);}

	if (trialCount==1 && (printStuff>2)) {cout << "\nA:\n"; A.write(cout);}

	Field F(defaultPrime>0 ? defaultPrime : 2);
	cout << "Testing with Z of type '";
	R.write(cout);
	cout<<"' and Z/pZ of type '";
	F.write(cout)<<"'"<<endl;

	typedef RationalSolver<Ring, Field, class RandomPrime, DixonTraits> QSolver;
	typedef DiophantineSolver<QSolver> ZSolver;

	//typedef std::vector<std::pair<RingElement, RingElement> > FractionVector;
	typedef VectorFraction<Ring> FractionVector;
	FractionVector x(R,c);
	int result=0;
	SolverLevel level = (SolverLevel)levelAsInt;

	for (int iteration=0; iteration<3; iteration++) {
		if (iteration==0 && !useDeterm) continue;
		if (iteration==1 && !useRandom) continue;
		if (iteration==2 && !useDiophantine) continue;

		// no more cleaning
		//clear x
		//for (FractionVector::Dense::iterator i=x.begin(); i!=x.end(); i++) {
		//	R.init(i->first, 0);
		//	R.init(i->second, 0);
		//}

		QSolver* rsolver;
		if (defaultPrime == 0)
			rsolver = new QSolver(R, LinBox::RandomPrime(primeBits));
		else
			rsolver = new QSolver(defaultPrime, R, LinBox::RandomPrime(primeBits));

		ZSolver zsolver(*rsolver);
		SolverReturnStatus s;

		if (iteration==0) {
			cout << "Solving deterministically.\n";
			s = zsolver.solve(x.numer, x.denom, A, b, numPrimes, level);
		}
		else if (iteration==1) {
			cout << "Solving randomly.\n";
			s = zsolver.randomSolve(x.numer, x.denom, A, b, numPrimes, level);
		}
		else {
			cout << "Solving diophantically.\n";
			s = zsolver.diophantineSolve(x.numer, x.denom, A, b, numPrimes, level);
		}
		cout << "solverReturnStatus: " << solverReturnString[(int)s] << "\n";

#ifdef RSTIMING
		rsolver->reportTimes(cout);
#endif
		if (s == SS_OK)	{
			VectorFraction<Ring> red(x);

			if (printStuff > 0) {
				if (useDiophantine){
					cout<<"Number of system solved   : "<<zsolver.numSolutionsNeeded<<endl;
					cout<<"Number of system failed   : "<<zsolver.numFailedCallsToSolver<<endl;
					cout<<"Number of system revelant : "<<zsolver.numRevelantSolutions<<endl;

				}
				cout<<"Reduced solution: ";
				integer tmp;
				size_t maxbits=0;
				for (int i=0;i<n;++i){
					R.convert(tmp,x.numer[i]);
					maxbits=(maxbits > tmp.bitsize() ? maxbits: tmp.bitsize());
				}
				R.convert(tmp,x.denom);
				cout<<"numerators hold over "<<maxbits<<" bits and denominators hold over "<<tmp.bitsize()<<" bits\n";
			}
			if (printStuff > 1) {
				if (useFiles){
					ofstream out("td-x.txt");
					out<< "Reduced solution: ";
					red.write(out) << "\n";
					out.close();
				}
				else{
					cout << "Reduced solution: ";
					red.write(cout) << "\n";
				}
			}
			Vector LHS(n), RHS(b);
			// check that Ax = b, if it thought it was okay
			MD.vectorMul(LHS, A, red.numer);
			VD.mulin(RHS, red.denom);
			if (VD.areEqual(LHS, RHS))
				cout << "Ax=b : Yes" << endl;
			else {
				cout << "Ax=b : No" << endl;
				if (level >= SL_LASVEGAS)
					cout << "ERROR: Las Vegas or Certified solver should never return wrong answer" << endl;
			}

			if (iteration==2 && level == SL_CERTIFIED) {
				// check certificate of minimality z
				// should satisfy that zA is integral, and den(z.b) == den(y)

				Integer dp, tmp, denzb;
				VectorFraction<Ring> z(zsolver.lastCertificate);
				R.init(dp, 0);
				typename Vector::iterator zi = z.numer.begin();
				typename Vector::iterator bi = b.begin();
				for (; bi != b.end(); bi++, zi++)
					R.addin(dp, R.mul(tmp, *bi, *zi));

				R.gcd(denzb, dp, z.denom);
				R.div(denzb, z.denom, denzb);

				VectorFraction<Ring> tmpvf(x);
				bool certified = R.areEqual(denzb, tmpvf.denom);
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
					cout << "Solution is certified correctly as having minimal denominator." << endl;
			}
		}
		else if (s==SS_INCONSISTENT && level == SL_CERTIFIED) {
			cout << "About to check certificate of inconsistency";
			VectorFraction<Ring> cert(zsolver.lastCertificate);
			if (printStuff > 1) {
				cout << ": ";
				cert.write(cout);
			}
			cout << endl;
			std::vector<Integer> certA(c);
			if (R.isZero(cert.denom))
				cout << "ERROR: Zero denom in inc-certificate. May not have been generated." << endl;

			Integer certb, tmp;

			for (int i=0; i<c; i++) R.init(certA[i], 0);
			R.init(certb, 0);

			for (int i=0; i<n; i++)
				for (int j=0; j<c; j++)
					R.addin(certA[j], R.mul(tmp, cert.numer[i], A[i][j]));
			for (int i=0; i<n; i++)
				R.addin(certb, R.mul(tmp, cert.numer[i], b[i]));

			bool certifies1 = true; //check certificate
			if (R.isZero(certb)) {
				cout << "ERROR: Product of certificate . b is zero!" << endl;
				certifies1 = false;
			}
			bool certifies2 = true;
			for (size_t i=0; certifies2 && i<A.rowdim(); i++)
				if (!certifies2) {
					certifies2 = false;
					cout << "ERROR: entry " << i << " of certificate . A is nonzero" << endl;
				}
			if (certifies1 && certifies2)
				cout << "System is certified correctly as inconsistent." << endl;
		}
		delete rsolver;
	}
	return result;

};

template <class Field>
int fieldTest()
{
	return
	0
	//+test<NTL_ZZ, Field>()
	+test<PID_integer, Field>()
	//+(testPidDouble?test<PID_double, Field>():0)
	// */
	;
};

void testAllFields()
{
	//fieldTest<GivaroZpz<Log16> >();
	//fieldTest<NTL_zz_p>();

	//fieldTest<GivaroZpz<Std16> >();

	//fieldTest<Modular<int> >();
	fieldTest<Modular<double> >();


	//fieldTest<GivaroZpz<Std32> >();           //broken?
	//fieldTest<GivaroZpz<Std64> >();           //broken?
	//fieldTest<GivaroGfq>();                   //broken?

	//fieldTest<GivaroMontg>();    // appears to be broken in current build
	//fieldTest<NTL_ZZ_p>();       // appears to be broken in current build
	//fieldTest<Modular<integer> >();

	// */
	// this takes a long time to compile with all fields
	// so comment out unused ones when debugging
}

void genTestData()
{
	bool* auxRow = new bool[n];
	int auxRows = 0;
	for (int i=0; i<n; i++) {
		auxRow[i] = (random_01() <= singularProportion);
		if (auxRow[i]) auxRows++;
	}
	cout << "at least " << auxRows << " dependent rows" << endl;

	if (inconsistent && auxRows == 0)
	{ auxRows++; auxRow[(int)(random_01()*n)] = true; }

	integer eBound(eBoundCmd);
	if (auxRows > 0 && auxRows < n)
		eBound /= (n-auxRows);
	if (eBound == 0) {
		cout << "WARNING: 'b' dropped to 0. Changed to 1, try increasing 'b'." << endl;
		eBound = 1;
	}

	PID_integer Z;
	PID_integer::RandIter ri(Z, 2*eBound,entrySeed); //for some reason this iterator tends to give numbers with
	//large common factors, so we perturb the data a bit
	bool notRandomEnough = (eBound >> 64) > 0;
	double bigStuff = ((long long)1)<<25;
	for (int i=0; i<n; i++) {
		ri.random(bentries[i]);
		bentries[i] -= eBound-1;
		if (notRandomEnough) bentries[i] += static_cast<int>((random_01()-0.5)*bigStuff);
	}

	for (int i=0; i<n*c; i++) {
		ri.random(Aentries[i]);
		Aentries[i] -= eBound-1;
		if (notRandomEnough) Aentries[i] += static_cast<int>((random_01()-0.5)*bigStuff);
	}

	int whichInconsistent = (int)(random_01()*auxRows);
	//make singular rows
	for (int i=0; i<n; i++)
		if (auxRow[i]) {
			for (int j=0; j<c; j++) Aentries[c*i+j] = 0;
			bentries[i] = 0;
			for (int k=0; k<n; k++)
				if (!auxRow[k]) {
					int m = (int)(random_01()*2)*2 - 1;
					for (int j=0; j<c; j++)
						Aentries[c*i+j] += Aentries[c*k+j]*m;
					bentries[i] += bentries[k]*m;
				}
			if (inconsistent) {
				if (whichInconsistent == 0)
					bentries[i] += (int)(random_01()*2)*2 - 1;
				whichInconsistent--;
			}
		}
	trialCount = 0; //so new data get printed

	int columnsToDestroy = destroyColumns;
	if (columnsToDestroy > c) {
		cout << "WARNING, o > c. Lowering o." << endl;
		columnsToDestroy = c;
	}

	for (int i=0; i<c; i++) {
		if (random_01()*(c-i-1) < columnsToDestroy) {
			for (int j=0; j<n; j++)
				Aentries[c*j+i] = 0;
			columnsToDestroy--;
		}
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
	if (!sparseMatrix)
		cout << "Matrix is dense with n="<<n<<" rows by c="<<c<<" columns\n";

	cout << "Seed: " << entrySeed <<"\n";
	srand(entrySeed);

	Aentries = new integer[n*c];
	bentries = new integer[n];

	if (useFiles && !sparseMatrix) {

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
		testAllFields();
		cout << "finished trial " << (j+1) << " of " << trials << endl;
	}

	delete[] Aentries;
	delete[] bentries;

	return 0;
}

//! @todo come up with better test data, so can have a big singular matrix of all 0..9
//! @todo change "probability of dependence" to "set X dependent rows"
//! @bug fixme : seems to not work for n >= 10000
