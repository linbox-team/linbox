
/*
 * examples/power_rank.C
 *
 * Copyright (C) 2012 J-G Dumas
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

/** \file examples/power_rank.C
 * @example  examples/power_rank.C
  \brief Rank of sparse matrix over Z or Zp.
  \ingroup examples
  */
#include "linbox/linbox-config.h"

#include <iostream>

#include "linbox/field/givaro.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/smith-form-sparseelim-poweroftwo.h"

using namespace LinBox;
using namespace std;



int main (int argc, char **argv) {
    commentator().setMaxDetailLevel (-1);
    commentator().setMaxDepth (-1);
    commentator().setReportStream (std::cerr);
    
    if (argc < 3 || argc > 3)
    {	cerr << "Usage: rank <matrix-file-in-supported-format> <power of two exponent>]" << endl; return -1; }
    
    ifstream input (argv[1]);
    if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }
    
    long unsigned int r;
    
    if (argc == 3) {
        LinBox::Timer tim; 
        size_t exponent = atoi(argv[2]);
        std::vector<std::pair<size_t,size_t> > local;
        if (exponent > 63) {
            typedef LinBox::PID_integer Ring;
            Ring ZZ;
            LinBox::MatrixStream<Ring> ms( ZZ, input );
            LinBox::SparseMatrix<Ring, LinBox::Vector<Ring>::SparseSeq > A (ms);
            input.close();
            LinBox::PowerGaussDomainPowerOfTwo< Givaro::Integer > PGD;
            tim.clear(); tim.start();
            PGD(local, A, exponent);
            tim.stop();
            
        } else {
            typedef LinBox::UnparametricField<int64_t> Ring;
            Ring R;
            LinBox::MatrixStream<Ring> ms( R, input );
            LinBox::SparseMatrix<Ring, LinBox::Vector<Ring>::SparseSeq > A (ms);
            input.close();
            LinBox::PowerGaussDomainPowerOfTwo< uint64_t > PGD;
            
            tim.clear(); tim.start();
            PGD(local, A, exponent);
            tim.stop();
        }
        

        std::cout << "Local Smith Form : (";
        for (std::vector<std::pair<size_t,size_t> >::const_iterator  p = local.begin();
             p != local.end(); ++p)
            std::cout << '[' << p->second << ',' << p->first << "] ";
        cout << ')' << endl;
        
        
        std::cerr << tim << std::endl;
    }
    
    return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

