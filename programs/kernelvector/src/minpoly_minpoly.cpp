/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
/*  Copyright (c) 2015 HPAC
 *  Written by 	Pascal Giorgi <Pascal.Giorgi@lirmm.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */


#define PROFILE_PMBASIS
#define LOW_MEMORY_PMBASIS
#define TRACK_MEMORY_MATPOL

// #define PROFILE_PMBASIS 1
//#define MEM_PMBASIS 1

#include <linbox/linbox-config.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/randiter/random-fftprime.h>
#include <linbox/randiter/random-prime.h> 
#include <linbox/ring/modular.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>
#include <linbox/algorithms/polynomial-matrix/order-basis.h>
#include <givaro/givtimer.h>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

using namespace LinBox;
using namespace std;

template<typename Field, typename Matrix>
Matrix* loadSerie(const Field& F, ifstream& in, bool trans){
	size_t M,N,D,row,col;
	in>>M>>N>>D;
	row=(trans?N:M);
	col=(trans?M:N);
	Matrix* S= new Matrix(F,row+col,col,D);
	for(size_t k=0;k<D;k++) // loop on the degree		
		for(size_t i=0;i<M;i++) // loop on the row of coefficient matrix
			for(size_t j=0;j<N;j++) // loop on the column of coefficient matrix
				if (!trans)
					in>>S->ref(i,j,k);
				else
					in>>S->ref(j,i,k);
	for(size_t i=0;i<col;i++)
		F.assign(S->ref(row+i,i,0),F.one);

	return S;
}

template<typename Field, typename Matrix>
void writeMatrixLinGen(const Field& F, size_t r, Matrix &Pol, vector<size_t> &shift, ofstream& out, bool trans){
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain, Field> PMatrix;
	size_t m=shift.size();
	size_t d=Pol.size();
	// convert to polynomial of matrices
	PMatrix Sigma(F,m,m,d);
	Sigma.copy(Pol);
	
	// take the r rows which have lowest shift (compute associated permutation)
	std::vector<size_t> Perm(m);
	for (size_t i=0;i<m;++i)
		Perm[i]=i;
	for (size_t i=0;i<m;++i) {
		size_t idx_min=i;
		for (size_t j=i+1;j<m;++j)
			if (shift[j]< shift[idx_min])
				idx_min=j;
		std::swap(shift[i],shift[idx_min]);
		Perm[i]=idx_min;
	}
	BlasPermutation<size_t> BPerm(Perm);
	BlasMatrixDomain<Field> BMD(F);;
	// Apply BPerm to  Sigma 
	for (size_t i=0;i<d;++i)
		BMD.mulin_right(BPerm,Sigma[i]);

	// Compute the reverse polynomial of Sigma according to row shift 
	size_t max= *std::max_element(shift.begin(),shift.begin()+r);
	PMatrix lingen(F,r,r,max+1);
	for (size_t i=0;i<r;i++)
		for (size_t j=0;j<=shift[i];j++)
			for (size_t k=0;k<r;k++)
				F.assign(lingen[shift[i]-j].refEntry(i,k), Sigma[j].getEntry(i,k));
#ifdef __VERBOSE
	std::cout<<"LinGen:="<<lingen<<"\n";
	std::cout<<"f:=LinearAlgebra:-Determinant(LinGen): f:=f/lcoeff(f) mod "<<F.cardinality()<<";\n";
#endif
	// output the lingen
	out<<r<<" "<<r<<" "<<max+1<<endl;
	if (!trans) {
		for(size_t k=0;k<max+1;k++){
			for (size_t i=0;i<r;i++)
				for(size_t j=0;j<r;j++)
					F.write(out,lingen[k].getEntry(i,j))<<" ";
			out<<endl;
		}
	}
	else {
		for(size_t k=0;k<max+1;k++){
			for(size_t j=0;j<r;j++)
				for (size_t i=0;i<r;i++)		
					F.write(out,lingen[k].getEntry(i,j))<<" ";
			out<<endl;
		}
	}
}

template<typename Field>
void launch_computation(integer &p, string in, string out, bool trans){
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain, Field> PMatrix;
	Field F(p);
	ifstream input(in);
	ofstream output(out);
	// load the sequence
    std::cerr << "[MINP] load Serie modulo " << p << " ..." << std::endl;
    Givaro::Timer Chrono, GlobTime; Chrono.clear();Chrono.start();
 	PMatrix *Serie=loadSerie<Field,PMatrix>(F,input,trans);
    Chrono.stop();
    std::cerr << "[MINP] Serie loaded in " << Chrono << std::endl;
    GlobTime = Chrono;

    size_t m,n,d;
	m=Serie->rowdim();
	n=Serie->coldim();
	d=Serie->size();

    std::cerr << "[MINP] Order basis (" << m << 'x' << n << ")x" << d << " ..." << std::endl;
    Chrono.clear();Chrono.start();

	// set the shift to [ 0 .. 0 1 .. 1]
	vector<size_t> shift(m,0);
	std::fill(shift.begin()+m-n,shift.end(),1);
	
#ifndef  LOW_MEMORY_PMBASIS
	PMatrix Sigma(F,m,m,d+1);
	// compute order basis 
    Chrono.stop();
    GlobTime += Chrono;
    std::cerr << "[MINP] Sigma/shifts initialized in " << Chrono << std::endl;
    Chrono.clear();Chrono.start();
	OrderBasis<Field> SB(F);
	SB.PM_Basis(Sigma, *Serie, d, shift);
#ifdef __VERBOSE
	std::cout<<"p:="<<F.cardinality()<<";"<<std::endl;
	std::cout<<"d:="<<Serie->size()<<";"<<std::endl;
	std::cout<<"Serie:="<<*Serie<<std::endl;
	std::cout<<"Sigma:="<<Sigma<<std::endl;
#endif
	delete Serie;	
    Chrono.stop();
    std::cerr << "[MINP] Basis computed in " << Chrono << std::endl;
    GlobTime += Chrono;
	// write the lingen according to order basis
	writeMatrixLinGen(F, m-n, Sigma, shift, output,trans);

    std::cerr << "[MINP] Minpoly Done in " << GlobTime << std::endl;
	
#else

	// OPTIMIZED MEMORY VARIANT
	
	PMatrix *Sigma;
	// compute order basis 
    Chrono.stop();
    GlobTime += Chrono;
    std::cerr << "[MINP] Sigma/shifts initialized in " << Chrono << std::endl;
    Chrono.clear();Chrono.start();
	OrderBasis<Field> SB(F);
	SB.PM_Basis_low(Sigma, Serie, d, shift);
#ifdef __VERBOSE
	std::cout<<"p:="<<F.cardinality()<<";"<<std::endl;
	std::cout<<"d:="<<Serie->size()<<";"<<std::endl;
	//std::cout<<"Serie:="<<*Serie<<std::endl;
	std::cout<<"Sigma:="<<Sigma<<std::endl;
#endif
    Chrono.stop();
    std::cerr << "[MINP] Basis computed in " << Chrono << std::endl;
    GlobTime += Chrono;
	// write the lingen according to order basis
	writeMatrixLinGen(F, m-n, *Sigma, shift, output,trans);
	delete Sigma;
    std::cerr << "[MINP] Minpoly Done in " << GlobTime << std::endl;

#endif

	
	

	
}

int main(int argc, char** argv){
    std::cerr << "command-line:";
    for(int i = 0; i < argc; ++i)
       std::cerr << ' ' << argv[i];
    std::cerr << std::endl;
	
	static integer p = 0; // prime field caharacteristic
	static string input="input.txt";
	static string output="output.txt";
	static bool   right=false;

	static Argument args[] = {
		{ 'p', "-p P"  , "Set the characteristic of the prime field."        , TYPE_INTEGER, &p },
		{ 'i', "-i IN" , "Set the input file of the sequence ."              , TYPE_STR    , &input },
		{ 'o', "-o OUT", "Set the output file for matrix lingen."            , TYPE_STR    , &output },
		{ 'r', "-r "   , "Compute a right matrix lingen (default is left)."  , TYPE_BOOL   , &right },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	
	if (p == 0)
	    throw std::runtime_error("No field characteristic specified with option -p");
	    	
	typedef Givaro::Modular<double>              SmallField;	
	typedef Givaro::Modular<Givaro::Integer>     LargeField;
	typedef Givaro::Modular<RecInt::ruint128,RecInt::ruint256>  Field128;

	
	if (p.bitsize() < 26) {
		launch_computation<SmallField>(p,input,output,right);
	}
	else {
		if (p.bitsize() < 128)
			launch_computation<Field128>(p,input,output,right);
		else
			launch_computation<LargeField>(p,input,output,right);
	}
	

	return 0;
}


