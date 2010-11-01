/** 
 * examples/solve.C
 *
 * Copyright (C) 2005, 2010 J-G Dumas, D. Saunders, P. Giorgi
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{

 	commentator.setMaxDetailLevel (-1);
 	commentator.setMaxDepth (-1);
 	commentator.setReportStream (std::cerr);


    if (argc < 2 || argc > 4) {
        cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>] [<p>]" << endl;
        return 0;
    }
    srand48( BaseTimer::seed() );

    std::ifstream input (argv[1]);
    if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }
    std::ifstream invect;

    bool createB = false;
    int ModComp = 0;
    if (argc == 2) {
        createB = true;
        ModComp = 0;
    }

    if (argc == 3) {
        invect.open (argv[2], std::ifstream::in);
        if (!invect) { 
            createB = true;
            ModComp = 2;
        } else {
            createB = false;
            ModComp = 0;
        }
    }       

    if (argc == 4) {
        ModComp = 3;
        invect.open (argv[2], std::ifstream::in);
        if (!invect) { 
            createB = true;
        } else
            createB = false;
    }
            

    if (ModComp) {
            
        typedef Modular<double> Field;
        double q = atof(argv[ModComp]);
        Field F(q);
        MatrixStream< Field > ms ( F, input );
        SparseMatrix<Field> A (ms);  // A.write(std::cout);
        cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;
            
        std::vector<Field::Element> X( A.coldim()),B(A.rowdim());
        if (createB) {
            cerr << "Creating a random {-1,1} vector " << endl;
            std::vector<Field::Element> U( A.coldim() );
            for(std::vector<Field::Element>::iterator it=U.begin();
                it != U.end(); ++it)
                if (drand48() <0.5)
                    F.init(*it,-1);
                else
                    F.init(*it,1);
            A.apply(B,U);
        } else {
            for(std::vector<Field::Element>::iterator it=B.begin();
                it != B.end(); ++it)
                invect >> *it;
        }

//         A.write(std::cout << "A: ") << std::endl;

        std::cout << "B is [";
        for(std::vector<Field::Element>::const_iterator it=B.begin();it != B.end(); ++it)
            F.write(cout, *it) << " ";
        std::cout << "]" << std::endl;
                
        Timer chrono; 

            // Sparse Elimination 
        chrono.clear();
        chrono.start();		
        solve (X, A, B, Method::SparseElimination());
        chrono.stop();
		
        std::cout << "(Sparse Gauss) Solution is [";
        for(std::vector<Field::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            F.write(cout, *it) << " ";
        std::cout << "]" << std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;
		
           // BlasElimination
        chrono.start();		
        solve (X, A, B, Method::BlasElimination());
        chrono.stop();
		
        std::cout << "(BlasElimination) Solution is [";
        for(std::vector<Field::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            F.write(cout, *it) << " ";
        std::cout << "]" << std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;
		
            // Wiedemann 
        chrono.clear();
        chrono.start();		
        solve (X, A, B, Method::Blackbox());
        chrono.stop();
		
        std::cout << "(Wiedemann) Solution is [";
        for(std::vector<Field::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            F.write(cout, *it) << " ";
        std::cout << "]" << std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<<std::endl;;
		
 //             // Lanczos
//         chrono.clear();
//         chrono.start();		
//         solve (X, A, B, Method::Lanczos());
//         chrono.stop();
		
//         std::cout << "(Lanczos) Solution is [";
//         for(std::vector<Field::Element>::const_iterator it=X.begin();it != X.end(); ++it)
//             F.write(cout, *it) << " ";
//         std::cout << "]" << std::endl;		
//         std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;


//             // Block Lanczos
//         Method::BlockLanczos MBL;
//         MBL.preconditioner(Specifier::FULL_DIAGONAL);
//         chrono.clear();
//         chrono.start();		
//         solve (X, A, B, MBL);
//         chrono.stop();
		
//         std::cout << "(Block Lanczos) Solution is [";
//         for(std::vector<Field::Element>::const_iterator it=X.begin();it != X.end(); ++it)
//             F.write(cout, *it) << " ";
//         std::cout << "]" << std::endl;		
//         std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl<< std::endl;
		
    } else { 

        PID_integer ZZ;
        MatrixStream< PID_integer > ms( ZZ, input );
        SparseMatrix<PID_integer> A (ms);
        PID_integer::Element d;
        std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

        std::vector<PID_integer::Element> X( A.coldim()),B(A.rowdim());

        if (createB) {
            cerr << "Creating a random {-1,1} vector " << endl;
            for(std::vector<PID_integer::Element>::iterator it=B.begin();
                it != B.end(); ++it)
                if (drand48() <0.5)
                    *it = -1;
                else
                    *it = 1;
        } else {
            for(std::vector<PID_integer::Element>::iterator it=B.begin();
                it != B.end(); ++it)
                invect >> *it;
        }


        std::cout << "B is [";
        for(std::vector<PID_integer::Element>::const_iterator it=B.begin();
            it != B.end(); ++it)
            ZZ.write(cout, *it) << " ";
        std::cout << "]" << std::endl;
                
	
        Timer chrono; 
		
	// Wiedemann
        chrono.start();
        solve (X, d, A, B, Method::Wiedemann());
        chrono.stop();
		
        std::cout << "(Wiedemann) Solution is [";
        for(std::vector<PID_integer::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            ZZ.write(cout, *it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d) << std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

	// BlasElimination
        chrono.start();
        solve (X, d, A, B, Method::BlasElimination());
        chrono.stop();

        std::cout << "(BlasElimination) Solution is [";
        for(std::vector<PID_integer::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            ZZ.write(cout, *it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

	// Sparse Elimination
        chrono.start();
        solve (X, d, A, B, Method::SparseElimination());
        chrono.stop();

        std::cout << "(SparseElimination) Solution is [";
        for(std::vector<PID_integer::Element>::const_iterator it=X.begin();it != X.end(); ++it)
            ZZ.write(cout, *it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;		
        std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;


//             // Lanczos
//         chrono.start();
//         solve (X, d, A, B, Method::Lanczos());
//         chrono.stop();

//         std::cout << "(Lanczos) Solution is [";
//         for(std::vector<PID_integer::Element>::const_iterator it=X.begin();it != X.end(); ++it)
//             ZZ.write(cout, *it) << " ";
//         std::cout << "] / ";
//         ZZ.write(std::cout, d) << std::endl;		
//         std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
		

//             // Block Lanczos
//         chrono.clear();
//         chrono.start();
//         solve (X, d, A, B, Method::BlockLanczos());
//         chrono.stop();

//         std::cout << "(Block Lanczos) Solution is [";
//         for(std::vector<PID_integer::Element>::const_iterator it=X.begin();it != X.end(); ++it)
//             ZZ.write(cout, *it) << " ";
//         std::cout << "] / ";
//         ZZ.write(std::cout, d) << std::endl;		
//         std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;

    }

    return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
