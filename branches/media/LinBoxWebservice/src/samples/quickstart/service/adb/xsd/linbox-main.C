/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**  * Time-stamp: <09 Mar 07 17:05:49 Jean-Guillaume.Dumas@imag.fr>
 */
/**\file examples/linbox-functions.C 
\brief Determinant and rank of dense matrix over Z, 
matrix from file, det or rank to a file. 
\ingroup examples
*/

#include <iostream>
#include "linboxfunctions.h"

int main()
{

bool err = false; 
try
{ detFiles("m1000", "foodet"); } 
catch (...)
{ err = true; 
			 std::cerr << "det_function error occured" << std::endl;
}
if (!err) std::cerr << "def_function good result" << std::endl;

rankFiles("m1000", "foorank");
}
