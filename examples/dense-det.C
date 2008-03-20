/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/** \file dense-det.C examples/dense-det.C \ingroup examples
 *
 * @author Clément Pernet <clement.pernet@imag.fr>
 *
 * \brief Small program that computes timings for dense-charpoly computation of dense matrices
 *
 * Load the input matrix from a file, compute its charpoly over a prime finite field.
 */
// Copyright (C) 2004 Clément Pernet
// See COPYING for license information.
#include "linbox/linbox-config.h"
#define DISP 0
#define DEBUG 0
#include <fstream>
#include <iostream>
#include <vector>
//#include "Matio.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
#include <linbox/field/givaro-zpz.h>
#include "linbox/blackbox/sparse.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"

using namespace LinBox;
//typedef  GivaroZpz<Std32> Field;
typedef Modular<double> Field;
typedef Field::Element Element;
typedef BlasMatrix<Element> Matrix;
typedef std::vector<Field::Element> Polynomial;

template<class T, template <class T> class Container>
std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
	for(typename Container<T>::const_iterator refs =  C.begin();
	    refs != C.end() ;
	    ++refs )
		o << (*refs) << " " ;
	return o << std::endl;
}

/// dense-det field-characteristic matrix-file number-of-computations
int main (int argc, char **argv)
{
	
	if (argc != 4){
		std::cerr << " Usage: dense-charpoly p Afile i" << std::endl
		     << " p: the characteristic of the field" << std::endl
		     << " Afile: name of file containing a square matrix" << std::endl
		     << " i: the number of iterations" << std::endl;
		return -1;
	}

	Field F (atoi(argv[1]));
	std::ifstream input (argv[2]);
	int nbi = atoi(argv[3]);
	
	if (!input) {
		std::cerr << "Error: Cannot load matrix " << argv[1] << std::endl;
		return -1;
	}
	std::cerr<<"Loading Matrix ....";
	
	SparseMatrix<Field> Ad(F);
	Ad.read (input);
	if ( Ad.coldim() != Ad.rowdim() ) {
		std::cerr << "A is " << Ad.rowdim() << " by " << Ad.coldim() << std::endl;
		std::cerr << "Error: A is not a square matrix" << std::endl;
		return -1;
	}
	std::cerr << "..";
	BlasMatrix<Element> A(Ad);
	std::cerr << "Done" << std::endl;
	BlasMatrixDomain<Field> BMD(F);
	
	Timer tim, tot;
	tot.clear();
	Element d ;
        for (int i=0; i<nbi;++i){
		tim.clear();
	        tim.start();
	        d = BMD.detin(A);
	        tim.stop();
	        tot+=tim;
	}

	F.write( std::cout , d ) << std::endl;
	std::cerr << tot << std::endl;
	
	return 0;
}
