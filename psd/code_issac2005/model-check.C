/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */

/* lieatlas/codes/model-check.C */
#include "linbox-config.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "linbox/field/gmp-rational.h"
#include "linbox/field/modular-int.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/vector/vector-domain.h"
#include <linbox/solutions/minpoly.h>
#include <linbox/solutions/rank.h>
#include "read.h"
#include "makeop.h"

using namespace std;
using namespace LinBox;

template<class BL>
bool ordercheck(const BL& op);

template<class BB>
ostream& matrix_write(ostream& out, BB* lbb, typename BB::Field& K);
	
int main(int argc, char* argv[]) {   
	//int prime = 1000000007;
	if (argc < 2 || argc > 4) {	
		cerr << endl << 
			"Usage: a.out <n> [checks]" 
			 << endl << endl << 
			"The E8 model Q<n>Cpt is checked (if file exists) to catch errors in"
			 << endl << 
			"the model and rootData files and/or in reading them." << endl;
		cerr << 
			"It is assumed that, eg., model 2 is in file `Q2Cpt' in the current directory"
			 << endl << 
			"or in the subdirectory `models' and that the root data is either in `e8data.m'"
			 << endl << 
			"in the current directory or in subdir `rootData'.  The checks are:" << endl;
		cerr << endl <<
			"1.  Operator products should be as expected from the Dynkin diagram:"
			 << endl << endl <<
			"            2" << endl <<
			"            |" << endl <<
			"    1 - 3 - 4 - 5 - 6 - 7 - 8" << endl << endl <<
			"	X(i)X(j) has order 1 if i = j, 2 if i adjacent to j, 3 otherwise."
			 << endl << endl <<
			"2.  The operator constructed from the zero facet should be F."
			 << endl << endl <<
			"4.  The operator constructed from the facet f := (0, 1, 2, 3, 4, 5, 6, 23)"
			 << endl <<
			"    Should be the zero matrix."
			 << endl << endl <<
			"8.  The operator constructed from the facet f/40 should be a symmetric"
			 << endl <<
			"    non-singular matrix  (also pos def - but we don't check that here)."
			 << endl;
		return 0;
	}

	int mn = atoi(argv[1]);
	int checks = argc == 3 ? atoi(argv[2]) : 15/* all checks */;

	bool pass = true;
    typedef GMPRationalField Rationals;
    typedef vector<Rationals::Element> Q_Vector;
	typedef SparseMatrix<Rationals> Blackbox;
    typedef Modular<int> Field;
	typedef vector<Field::Element> K_Vector;
	typedef SparseMatrix<Field> FBlackbox;
	typedef LieMatrix<Rationals, Blackbox> RBlackbox;
	typedef LieMatrix<Field, FBlackbox> LBlackbox;
	int prime = 9973, rk = 8;
	Field K(prime); Rationals Q; 
	VectorDomain<Rationals> QVD(Q);
	VectorDomain<Field> KVD (K);
	vector<Q_Vector> root(rk);
	Q_Vector rho(rk);
	vector<Blackbox*> op(rk);
	Blackbox* F; LBlackbox* lbb;
	Q_Vector nu(rk); Q_Vector r, coeff2; vector<int> i_l;
    integer c; Rationals::Element one; Q. init (one, 1);
	Q_Vector::iterator r_p, coeff2_p;

	cout << "Read rootData\n";
	string rootData("e8data.m");
	read_rootData(root, rho, rootData, Q, rk);

	ostringstream smn; smn << mn;
	string modelf("Q" + smn.str() + "Cpt");
	read_model(op, F, modelf, Q, rk);

	int dim = F -> rowdim();
	K_Vector w1(dim), w2(dim), v(dim);
	cout << "Read model " << smn.str() << " of dimension " << dim << endl;
    // initialize v at random.
    Field::RandIter RI(K);
    for (int i = 0; i < dim; ++i) RI.random(v[i]);
    /* model check
	 * . check operator orders
	 * . check that at lambda = zero vector, the produced operator is F,
	 * and is symmetric.
	 * . check that at special facet the produced operator is the zero matrix
	 * . check that at desired facet the produced operator is symmetric mod
	 * a prime and full rank.
	 * . In each case above check that a product of 120 ops was made,
	 * Show the parameters.

	 * This checks model reading and consistency with it's own F.
	 * This checks rootData reading and usage as well.
	*/

	if(checks & 1){  // test orders
		cout << "Table of orders of operator products X_i X_j" << endl; 
		if (ordercheck(op)) cout << "orders check OK" << endl;
		else { cout << "orders check FAIL" << endl; pass = false; }
	}

    if(checks & 14) // test at least one facet
	    cout << "Operators will be constructed from facets and tested mod " 
			 << prime << endl;


	if(checks & 2){ // test nu = zero vector 
		cout << endl << "For facet ";
		QVD.write(cout, nu) << ", ";
		compute_r (r, i_l, nu, root, rho, Q);

		coeff2. resize (r. size());
		for (r_p = r. begin(), coeff2_p = coeff2. begin(); 
			 r_p != r. end(); ++ r_p, ++ coeff2_p) { 
			 Q. add (*coeff2_p, *r_p, one); 
			 Q. invin (*coeff2_p);
		} 
		RBlackbox* Lbb = new RBlackbox(op, i_l, r, coeff2, F, Q); 

		if (r.size() != 120) { 
			pass=false; cout << "Error in coefficients construction" << endl; }

		FBlackbox* FP;
		MatrixMod::mod (FP, *F, K); 
		MatrixMod::mod (lbb, *Lbb, K);

		if (lbb->coldim() <= 28){
			cout << endl;
			matrix_write(cout, lbb, K);
		}

		// mul lbb and F by random vector v and compare results
		lbb->apply(w1, v); 
		FP -> apply(w2, v);
		if ( KVD.areEqual(w1, w2)) {
			cout << "correct: operator equals F" << endl;
		} 
		else {  
			cout << "wrong: operator should equal F" << endl; pass = false;
		}
		delete lbb;
	}

	if(checks & 4){ //test at nu = 0 1 2 3 4 5 6 23;
		for (int i = 0; i < rk; ++i) Q.init(nu[i] , i); Q.init(nu[rk-1],23);
		cout << endl << "For facet ";
		QVD.write(cout, nu) << ", ";
		compute_r (r, i_l, nu, root, rho, Q);
		coeff2. resize (r. size());
		for (r_p = r. begin(), coeff2_p = coeff2. begin(); 
			 r_p != r. end(); ++ r_p, ++ coeff2_p) { 
			 Q. add (*coeff2_p, *r_p, one); 
			 Q. invin (*coeff2_p);
		} 
		RBlackbox* Lbb = new RBlackbox(op, i_l, r, coeff2, F, Q); 

		if (r.size() != 120) { 
			pass=false; cout << "Error in coefficients construction" << endl; }

		MatrixMod::mod (lbb, *Lbb, K);
		if (lbb->coldim() <= 28){
			cout << endl;
			matrix_write(cout, lbb, K);
		}
		// mul lbb by random vector and check for zero.
		lbb->apply(w1, v);
		// cout << endl;
		// KVD.write(cout, v) << " is v" << endl;
		// KVD.write(cout, w1) << " is w1"  << endl;
		if ( KVD.isZero(w1)) {
			cout << "correct: operator is the zero matrix" << endl;
		} else {  
			cout << "wrong: operator should be zero matrix" << endl; pass = false;
		}
		delete lbb;
	}

	if(checks & 8){ //test at nu = (1/40)*[0 1 2 3 4 5 6 23].
		for (int i = 0; i < rk; ++i) Q.init(nu[i] , i); Q.init(nu[rk-1],23);
		Rationals::Element s; Q.init(s, 40); Q.invin(s); // s = 1/40.
		QVD.mulin(nu, s);// nu = nu/40;

		cout << endl << "For facet ";
		QVD.write(cout, nu) << ", ";
		compute_r (r, i_l, nu, root, rho, Q);
		coeff2. resize (r. size());
		for (r_p = r. begin(), coeff2_p = coeff2. begin(); 
			 r_p != r. end(); ++ r_p, ++ coeff2_p) { 
			 Q. add (*coeff2_p, *r_p, one); 
			 Q. invin (*coeff2_p);
		} 
		RBlackbox* Lbb = new RBlackbox(op, i_l, r, coeff2, F, Q); 

		if (r.size() != 120) { 
			pass=false; cout << "Error in coefficients construction" << endl; }
	
		MatrixMod::mod(lbb, *Lbb, K);
		if (lbb->coldim() <= 28){
			cout << endl;
			matrix_write(cout, lbb, K);
		}

		// mul lbb by random vector and lbb^T by v and compare
		lbb->apply(w1, v); lbb->applyTranspose(w2, v);
		if ( KVD.areEqual(w1, w2)) {
			cout << "correct: symmetric operator" << endl;
		} else {  
			cout << "wrong: operator is UNsymmetric" << endl; pass = false;
		}

		cout << "Computing rank:";
		unsigned long d;
		rank (d, *lbb, K);
		cout << "rank: " << d << endl;
		if ( d == lbb->rowdim() ) { 
			cout << "\rOperator of full rank" << endl;
		} else {  
			cout << "\rOperator NOT full rank" << endl; pass = false;
			matrix_write(cout, lbb, K) << endl;
		}
		delete lbb;
	}

	cout << "    model " << mn;
	if (pass && checks == 15) cout << " passes all tests " << endl;
	else if (pass && checks < 15) {
		cout << " passes tests"; 
		if (checks & 1) cout << " - orders";
		if (checks & 2) cout << " - zero-facet->F" ;
		if (checks & 4) cout << " - 40test-facet->zero";
		if (checks & 8) cout << " - test-facet->sym-rank";
		cout << endl;
	}
	else cout << " FAILS one or more tests" << endl;
	return 0;
}

template<class BB>
ostream& matrix_write(ostream& out, BB* lbb, typename BB::Field& K){
	typedef typename BB::Field Field;
	typedef typename Field::Element Element;
	VectorDomain<Field> V(K);
	Element zero, one;
	K.init(zero, 0); K.init(one, 1);
	int n = lbb->coldim();
	vector<Element> v(n, zero), w(n, zero);
	for( int i = 0; i < n; ++i){
		K.assign(v[i], one);
		lbb->apply(w, v);
		K.assign(v[i], zero);
		V.write(out, w) << endl;
	}
	return out;
}

template<class BB>
int productorder(const BB& A, const BB& B, size_t max = 3){
	typedef typename BB::Field Field;
	const Field& F(A.field());
	typedef vector<typename Field::Element> Vector;
	VectorDomain<Field> VD(F);
	typename Field::Element one, zero;
	F.init(one, 1); F.init(zero, 0);
	int n = B.coldim(), k, kb;
	Vector e(n), w(n), w2(n);
	for (int i = 0; i < n; ++i) F.init(e[i], random()%1000000);  
	//{ F.init(e[i], random()%1000000);  F.write(cout, e[i]) << endl;}
	//{
	//	if (i > 0) e[i-1] = zero;
	//	e[i] = one;
	// for what k > 0, (AB)^k e == e ?
	w = e;
	k = 0;
	do{
		B.apply(w2, w);
		A.apply(w, w2);
		++k;
	} while (k <= max && !VD.areEqual(w, e));
	return k > max ? 0 : k;
	//	if (k > max) return 0;
	//	if (i == 0) kb = k;
	//	if (i > 0)  kb = lcm(k, kb);
	//} 
	//return kb;
}

bool falsep(int i, int j){
	cout << "wrong order at " << i << ", " << j << endl;
	return false;
}

template<class BL>
bool ordercheck(const BL& op){
    int rk = op.size();

    for (int i = 1; i <= rk; ++i){
	    for (int j = 1 ; j <= rk; ++j) {
			int o = productorder(*(op[i-1]),*(op[j-1])); 
			if (o == 0) cout << "- "; else cout << o << " "; 
			if (1 == o) continue;
			if (i == j){ 
				if ( 1 != o)
					return falsep(i, j); }
            else if (3 <= i && j == i+1 ){ 
				if ( 3 != o)
					return falsep(i, j); }
            else if (3 <= j && i == j+1 ){ 
				if ( 3 != o)
					return falsep(i, j); }
            else if (i*j == 3 ){ 
				if ( 3 != o)
					return falsep(i, j); }
            else if ( ((i == 2 && j == 4 )||(i == 4 && j == 2))  ){ 
				if ( 3 != o)
					return falsep(i, j); }
            else { 
				if ( 2 != o)
					return falsep(i, j); }
		}
		cout << endl;
	}
	return true;
}
