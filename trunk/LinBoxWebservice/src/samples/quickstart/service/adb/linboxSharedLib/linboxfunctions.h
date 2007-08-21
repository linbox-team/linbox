/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**  * Time-stamp: <09 Mar 07 17:05:49 Jean-Guillaume.Dumas@imag.fr>
 */
/**\file examples/linbox-functions.C 
\brief Determinant and rank of dense matrix over Z, 
matrix from file, det or rank to a file. 
\ingroup examples
*/

#include <iostream>

//#include "linbox/blackbox/dense.h"
//#include "linbox/solutions/det.h"
//#include "linbox/solutions/rank.h"

void det(std::istream& matrix_in, std::ostream& det_out);

void detFiles(const char* matfile, char* outfile);

void rank(std::istream& matrix_in, std::ostream& rank_out);

void rankFiles(const char* matfile, char* outfile);
