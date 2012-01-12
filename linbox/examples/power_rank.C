/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/*
 * examples/power_rank.C
 *
 * Copyright (C) 2005, 2010 J-G Dumas
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** \file examples/power_rank.C
 * @example  examples/power_rank.C
  \brief Rank of sparse matrix over Z or Zp.
  \ingroup examples
  */
#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/field/givaro.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/smith-form-sparseelim-local.h"

using namespace LinBox;
using namespace std;



int main (int argc, char **argv)
{
	commentator.setMaxDetailLevel (-1);
	commentator.setMaxDepth (-1);
	commentator.setReportStream (std::cerr);

	if (argc < 4 || argc > 4)
	{	cerr << "Usage: rank <matrix-file-in-supported-format> <prime> <prime-power>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 4) {
		LinBox::int64_t p = atoi(argv[2]);
		LinBox::int64_t q = atoi(argv[3]);
		typedef GivaroZpz<Std64> Field;
		Field F(q);
		MatrixStream<Field> ms( F, input );
		SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;

		// using Sparse Elimination
		PowerGaussDomain< Field > PGD( F );
		std::vector<std::pair<size_t,size_t> > local;

		Timer tq; tq.clear(); tq.start();
		PGD(local, B, q, p);
		tq.stop();


		std::cout << "Local Smith Form : (";
		for (std::vector<std::pair<size_t,size_t> >::const_iterator  p = local.begin();
		     p != local.end(); ++p)
			std::cout << p->first << " " << p->second << ", ";
		cout << ")" << endl;


		std::cerr << tq << std::endl;
	}

	return 0;
}
