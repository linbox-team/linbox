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

#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/smith-form-iliopoulos2.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::TPL> SparseMat;

typedef InvariantFactors<Field,SparseMat> FactorDomain;
typedef typename FactorDomain::PolyDom PolyDom;
typedef typename FactorDomain::PolyRing PolyRing;
typedef DenseVector<PolyRing> FactorVector;

int main(int argc, char** argv)
{
	int earlyTerm;
	int p = 3, b = 3;
	std::string mFname,oFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
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

	std::cout << "Finished reading" << std::endl;

	PolyDom PD(F,"x");
	PolyRing R(PD);
	FactorVector factorList(R);
	FactorDomain CIF(F, R);
	IliopoulosDomain<PolyRing> ID(R);

	CIF.computeFactors(factorList, M, b, earlyTerm);
	std::cout << "Finished computing factors" << std::endl;

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


