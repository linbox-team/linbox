/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**  * Time-stamp: <09 Mar 07 17:05:49 Jean-Guillaume.Dumas@imag.fr>
 */
/**\file examples/linbox-functions.C 
\brief Determinant and rank of dense matrix over Z, 
matrix from file, det or rank to a file. 
\ingroup examples
*/

#include <iostream>

#include "linbox/blackbox/dense.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"

void det_function(std::istream& matrix_in, std::ostream& det_out)
{
		typedef LinBox::PID_integer Integers;		
		Integers ZZ;

		LinBox::DenseMatrix<Integers> A(ZZ);
		A.read(matrix_in);

		Integers::Element det_A;
		LinBox::det (det_A, A);

		ZZ.write(det_out, det_A) << std::endl;
}

void det_function_files(const char* matfile, char* outfile)
{
		std::ifstream input (matfile);
		if (!input) throw std::bad_exception();

		std::ofstream output (outfile);
		if (!output) throw std::bad_exception();

		det_function(input, output);
		//output.close()
}

void rank_function(std::istream& matrix_in, std::ostream& rank_out)
{
		typedef LinBox::PID_integer Integers;		
		Integers ZZ;

		LinBox::DenseMatrix<Integers> A(ZZ);
		A.read(matrix_in);

		unsigned long rank_A;
		LinBox::rank (rank_A, A);

		rank_out << rank_A << std::endl;
}

void rank_function_files(const char* matfile, char* outfile)
{
		std::ifstream input (matfile);
		if (!input) throw std::bad_exception();

		std::ofstream output (outfile);
		if (!output) throw std::bad_exception();

		rank_function(input, output);
		//output.close()
}
/*
int main()
{

bool err = false; 
try
{ det_function_files("m1000", "foodet"); } 
catch (...)
{ err = true; 
			 std::cerr << "det_function error occured" << std::endl;
}
if (!err) std::cerr << "def_function good result" << std::endl;

rank_function_files("m1000", "foorank");
}
*/
