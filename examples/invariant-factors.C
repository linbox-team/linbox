/* linbox/algorithms/coppersmith-invariant-factors.h
 * Copyright (C) 2015 Gavin Harrison
 *
 * Written by Gavin Harrison <gavin.har@gmail.com>
 *
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"

#include "linbox/algorithms/invariant-factors.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

typedef InvariantFactors<Field,SparseMat> FactorDomain;
typedef typename FactorDomain::PolyDom PolyDom;
typedef typename FactorDomain::PolyRing PolyRing;
typedef DenseVector<PolyRing> FactorVector;

int main(int argc, char** argv)
{
	int earlyTerm;
	int p = 99991, b = 30, r = 3, t = 5;
	std::string mFname = "jumat",oFname="";

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
		{ 't', "-t T", "First blocking factor", TYPE_INT, &t},
		{ 'r', "-r R", "R-th factor used as mod for iliopoulos", TYPE_INT, &r},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &oFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	SparseMat M(F);

	{
		std::ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}

	std::cout << "Finished reading, dimension " << M.rowdim() << "x" << M.coldim() << ", nnz " << M.size() << std::endl;

	PolyDom PD(F,"x");
	PolyRing R(PD);
	FactorVector List1(R), List2(R), factorList(R);
	FactorDomain CIF(F, R);

#if 1
	CIF.solve(factorList, M, b, t, r, earlyTerm);
#else
	// from solve in invariant-factors.h
        // Compute first b1 factors
        CIF.computeFactors(List1, M, t, earlyTerm);

	std::cout << "Finished computing early factors" << std::endl;
	{
		for (int i = List1.size()-1; i >= 0; i--) {
			std::cout << PD.degree(List1[i]) << ", ";
		}
		std::cout << std::endl;
	}

        // get r-th factor
		PolyRing::Element d; R.assign(d, List1[t - r]);

        // Compute factors mod r-th factor
        CIF.computeFactors(List2, M, d, b, earlyTerm);

	// Fill in zeros before r-th factor with r-th factor
	    for (int i = 0; i < b - r; i++) {
			if (R.isZero(List2[i])) R.assign(List2[i], d);
		}

        // Fill in remaining factors with original values
		factorList.resize(b);
        for (int i = 0; i < r; i++) 
			factorList[i] = List1[t-1-i];
        for (int i = r; i < b; i++) 
			factorList[i] = List2[b-1-i];
#endif
	std::cout << "Finished computing factors" << std::endl;
		for (size_t i = 0; i<factorList.size(); i++) 
			std::cout << PD.degree(factorList[i]) << ", ";
		std::cout << std::endl;

	if (oFname.size() > 0) 
	{
		std::ofstream out(oFname);
		for (size_t i = 0; i<factorList.size(); i++) {
			R.write(out,factorList[i]);
			out << std::endl;
		}
		out.close();
	}

	return 0;
}


