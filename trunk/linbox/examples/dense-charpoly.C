/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// Copyright (C) 2004 Clément Pernet
// See COPYING for license information.
/** @name examples/dense-charpoly.C
 *
 * @author Clément Pernet <clement.pernet@imag.fr>
 *
 * @memo 
 * Small program that computes timings for dense-charpoly computation of dense matrices
 *
 * @doc
 * Load the input matrix from a file, compute its charpoly over a prime finite field.
 */
//@{
#include "linbox-config.h"
#define DISP 0
#define DEBUG 0
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
//#include "Matio.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
#include <linbox/field/givaro-zpz.h>
#include "linbox/blackbox/sparse.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"

using namespace LinBox;
using namespace std;
//typedef  GivaroZpz<Std32> Field;
typedef Modular<double> Field;
typedef Field::Element Element;
typedef BlasMatrix<Element> Matrix;
typedef vector<Field::Element> Polynomial;

template<class T, template <class T> class Container>
std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
	for(typename Container<T>::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o << std::endl;
}
template<class Field>
void 
print_poly(const Field & F, vector<typename Field::Element> & P){
	int i = 0;
	integer tmp;
	typename vector<typename Field::Element>::iterator it = P.begin();
	F.convert(tmp,*it++);
	while (tmp==0){ 
		F.convert(tmp,*it++); 
		i++;
	}
	cerr<<tmp;
	if (i) cerr<<"X^"<<i;
	i++;
	for (;it!=P.end(); it++,i++){
		F.convert(tmp,*it);
		if ( tmp!=0 ) cerr<<" + "<<tmp<<"X^"<<i; 
	}
	cerr<<endl;
}

template<class Field>
vector<typename Field::Element> &
mulpoly(const Field & F, 
	vector<typename Field::Element> &res, 
	const vector<typename Field::Element> & P1,
	const vector<typename Field::Element> & P2){
	size_t i,j;
	for (i=0;i<res.size();i++)
		F.init(res[i],0.0);
	for ( i=0;i<P1.size();i++)
		for ( j=0;j<P2.size();j++)
			F.axpyin(res[i+j],P1[i],P2[j]);
	return res;
		
}

/// dense-charpoly matrix-file field-characteristic number-of-computations
int main (int argc, char **argv)
{
	
	if (argc != 4){
		cerr << " Usage: dense-charpoly A p i" <<endl
		     << " p: the characteristic of the field"<<endl
		     << " A: a square matrix"<<endl
		     << " i: the number of iterations"<<endl;
		return -1;
	}

	Field F (atoi(argv[1]));
	ifstream input (argv[2]);
	int nbi = atoi(argv[3]);
	
	if (!input) {
		cerr << "Error: Cannot load matrix " << argv[1] << endl;
		return -1;
	}
	cerr<<"Loading Matrix ....";
	
	SparseMatrix<Field> Ad(F);
	Ad.read (input);
	size_t n=Ad.coldim();
	if ( Ad.coldim() != Ad.rowdim() ) {
		cerr << "A is " << Ad.rowdim() << " by " << Ad.coldim() << endl;
		cerr << "Error: A is not a square matrix" << endl;
		return -1;
	}
	cerr<<"..";
	BlasMatrix<Element> A(Ad);
	cerr<<"Done"<<endl;
	list<Polynomial> P;
	list<Polynomial>::iterator it;
	BlasMatrixDomain<Field> BMD(F);
	
	Timer tim, tot;
	tot.clear();
	for (int i=0; i<nbi;++i){
		tim.clear();
	        tim.start();
	        BMD.charpoly( P, A);
	        tim.stop();
	        tot+=tim;
	}


	Polynomial prod(n+1);
	Polynomial tmp(n+1);
 	it=P.begin();
	F.init(prod[0],1.0);
	
	list<pair<Polynomial,int> > factorList;
	list<pair<Polynomial,int> >::iterator factorList_it;
	size_t factornum=0;
	while(it!=P.end()){
#if DEBUG
		print_poly(F, *it);
#endif
		factornum++;
		factorList_it=factorList.begin();
		while(factorList_it!=factorList.end()){
			
			if (*it == factorList_it->first){
				factorList_it->second++;
				break;
			}
		  factorList_it++;
		}
		if (factorList_it==factorList.end()){
			pair<Polynomial,int> pa(*it,1);
			factorList.push_front(pa);
		}
		tmp = prod;
		mulpoly(F, prod, tmp, *it);
		it++;
	}
	
#if DISP
	cerr<<"----------------Factors-------------------"<<endl;
	factorList_it=factorList.begin();
	while(factorList_it!=factorList.end()){
		cerr<<factorList_it->second<<" ";
		print_poly(F, factorList_it->first);
		factorList_it++;  
	}
	cerr<<"----------------CharPoly(A)-------------------"<<endl;
	print_poly(F, prod );

	cerr<<"----------------Time--------------------------"<<endl;
#endif	
	
	
	cout<<factornum<<" "<< ((double) tot.usertime()) / nbi<<endl;
	
	return 0;
}
//@}
