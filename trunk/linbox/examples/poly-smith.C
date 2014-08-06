#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <linbox/linbox-config.h>

#include <linbox/util/commentator.h>
#include <linbox/field/modular.h>
#include <linbox/vector/stream.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/matrix/dense-matrix.h>

#include <linbox/algorithms/blackbox-block-container.h>
#include <linbox/algorithms/block-coppersmith-domain.h>

#include <givaro/givtimer.h>
#include <givaro/givzpz.h>
#include <givaro/givpoly1.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/algorithms/smith-form-direct.h>
#include <linbox/algorithms/smith-form-kannan-bachem.h>

using namespace LinBox;
using namespace std;

/// Prints out a matrix in maple format
template<class MatrixDom, class Matrix>
void printMatrix(MatrixDom MD, Matrix A)
{
	typename MatrixDom::Field::Element tmp;
	
	cout << "<";
	for (int i = 0; i < A.rowdim(); i++)
	{
		if (i != 0)
			cout << ",";
		
		MD.field().write(cout << "<", A.getEntry(tmp, i, 0));
		for (int j = 1; j < A.coldim(); j++)
		{
			MD.field().write(cout << "|", A.getEntry(tmp, i, j));
		}
		cout << ">";
	}
	cout << ">" << endl;
}

/// Prints a vector as a list
template<class Field, class Vector>
void printVector(Field F, Vector S)
{	
	F.write(cout << "[", S[0]);
	for (size_t i = 1; i < S.size(); i++)
		F.write(cout << ",", S[i]);
	cout << "]" << endl;
}

/// Makes vector polynomials monic
template<class PolyDom, class Vector>
void normalize(PolyDom &PD, Vector &S)
{
	typedef typename PolyDom::Type_t Scalar;
	
	for (size_t i = 0; i < S.size(); i++)
	{
		Scalar lcoef;
		PD.leadcoef(lcoef, S[i]);
		PD.divin(S[i], lcoef);
	}
}

int main(int argc, char **argv)
{
	typedef Givaro::ZpzDom<Givaro::Std32> BaseDom;
	typedef Givaro::Poly1Dom<BaseDom, Givaro::Dense> PolyDom;
	typedef GivaroPoly<PolyDom> Field;
	typedef MatrixDomain<Field> MatrixDom;
	typedef BlasMatrix<Field> Matrix;
	
	if (argc < 4)
	{
		cout << "Genertes a random GF(p)^(n-by-n) with polynomial coefficients with degree=d" << endl;
		cout << "Use: ./polysmith <p> <n> <d>" << endl;
		return 1;
	}
	
	int p = atoi(argv[1]);
	int n = atoi(argv[2]);
	int d = atoi(argv[3]);
	
	BaseDom BD(p);
	PolyDom PD(BD, "x");
	Field F(PD);
	MatrixDom MD(F);
	
	Matrix M(F, n, n);
	
	GivaroPolyRandIter<Field> Rand(F);
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Field::Element tmp;
			Rand.random(tmp, Givaro::Degree(d));
			M.setEntry(i, j, tmp);
		}
	}
	
	printMatrix(MD, M);
	
	// Text Book SF Domain
	// Comment out this block if things are getting too slow
	SmithFormDirectDomain<MatrixDom> SFD(MD);
	BlasVector<Field> S1(F, n, F.zero);
	
	Givaro::Timer t1;
	t1.start();
	SFD.solve(S1, M);
	t1.stop();
	normalize(PD, S1);
	printVector(F, S1);
	
	// Kannan-Bachem w/ Chou-Collins Improvement Domain 
	SmithFormKannanBachemDomain<MatrixDom> SFKB(MD);
	BlasVector<Field> S2(F, n, F.zero);
	
	Givaro::Timer t2;
	t2.start();
	SFKB.solve(S2, M);
	t2.stop();
	normalize(PD, S2);
	printVector(F, S2);
	
	cout << "SFD: " << n << " " << p << " " << d << " " << t1.usertime() << endl;
	cout << "SFK: " << n << " " << p << " " << d << " " << t2.usertime() << endl;
	
	return 0;
}
