/*
 * examples/smithvalence.C
 * Copyright (c) Linbox
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

/**\file examples/smithvalence.C
 * @example  examples/smithvalence.C
 \brief Valence of sparse matrix over Z or Zp.
 \ingroup examples
*/
#ifndef DISABLE_COMMENTATOR
#define DISABLE_COMMENTATOR
#endif

#define __VALENCE_REPORTING__ 1

#include <linbox/linbox-config.h>

#include <iostream>
#include <givaro/givintfactor.h>
#include <fflas-ffpack/paladin/parallel.h>
#include <linbox/solutions/smith-form.h>
#include <linbox/algorithms/smith-form-valence.h>
#include <vector>

using namespace LinBox;




int main (int argc, char **argv)
{

	if (argc < 2 || argc > 4) {
		std::cerr << "Usage: smithvalence <matrix-file-in-supported-format> [-ata|-aat|valence] [coprime]" << std::endl;
        std::cerr << "       Optional parameters valence and coprime are integers." << std::endl;
        std::cerr << "       Prime factors of valence will be used for local computation." << std::endl;
        std::cerr << "       coprime will be used for overall Z rank computation." << std::endl;
		return -1;
	}

    const std::string filename(argv[1]);

	std::ifstream input (filename);
	if (!input) { std::cerr << "Error opening matrix file " << filename << std::endl; return -1; }

	Givaro::ZRing<Integer> ZZ;
	MatrixStream< Givaro::ZRing<Integer> > ms( ZZ, input );
	typedef SparseMatrix<Givaro::ZRing<Integer>>  Blackbox;
	Blackbox A (ms);
	input.close();

	std::clog << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

    size_t method=0;
	Givaro::Integer val_A(0);

	if (argc >= 3) {
		if (strcmp(argv[2],"-ata") == 0) {
            method=2;
		}
		else if (strcmp(argv[2],"-aat") == 0) {
            method=1;
		}
		else {
			std::clog << "Suppose primes are contained in " << argv[2] << std::endl;
			val_A = LinBox::Integer(argv[2]);
		}
	}

    Givaro::Integer coprimeV=1;
	if (argc >= 4) {
		std::clog << "Suppose " << argv[3] << " is coprime with Smith form" << std::endl;
		coprimeV =Givaro::Integer(argv[3]);
	}

    std::vector<Givaro::Integer> SmithDiagonal;

#ifdef  __LINBOX_USE_OPENMP
	std::clog << "num procs: " << omp_get_num_procs() << std::endl;
    std::clog << "max threads: " << MAX_THREADS << std::endl;
#endif
        // Returns the Smith form as a diagonal,
        // the valence val_A,
        // the coprimeV used to compute the integral rank
	LinBox::Timer chrono; chrono.start();
    PAR_BLOCK {
        smithValence(SmithDiagonal, val_A, A, filename, coprimeV, method);
    }

	chrono.stop();

	std::clog << "Integer Smith Form :" << std::endl;
    writeCompressedSmith(std::cout, SmithDiagonal, ZZ, A.rowdim(), A.coldim()) << std::endl;

	std::clog << chrono << std::endl;
    

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
