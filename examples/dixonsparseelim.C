/* Copyright (C) The LinBox group
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

/**\file examples/dixonsparseelim.C
 @example examples/dixonsparseelim.C
 @author Jean-Guillaume.Dumas@univ-grenoble-alpes.fr
 * \brief Dixon System Solving Lifting using sparse LU
 *        (direct match to dixendenseLU)
 * \ingroup examples
 */
#include <iostream>
#include <omp.h>
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;
typedef Givaro::ZRing<Givaro::Integer> Ints;
typedef DenseVector<Ints> ZVector;


struct FixPrime {
    typedef Givaro::Integer Prime_Type;
    const Prime_Type _myprime;
    FixPrime(const Givaro::Integer& i) : _myprime(i) {}
    inline FixPrime &operator ++ () { return *this; }
    const Prime_Type &operator * () const { return randomPrime(); }
    const Prime_Type & randomPrime() const { return _myprime; }
    void setBits(uint64_t bits) {}
    template<class _ModField> void setBitsField() { }
};    


int main (int argc, char **argv) {
        // Usage
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>]" << std::endl;
        return 0;
    }

        // File
    std::ifstream input (argv[1]);
    if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

    
    std::ifstream invect;
    bool createB = false;
    if (argc == 2) {
        createB = true;
    }
    if (argc == 3) {
        invect.open (argv[2], std::ifstream::in);
        if (!invect) { 
            createB = true;
        } else {
            createB = false;
        }
    }       
    
        // Read Integral matrix from File
    Ints ZZ;
    MatrixStream< Ints > ms( ZZ, input );
    SparseMatrix<Ints> A (ms);
    Ints::Element d;
    std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
    
    {
            // Print Matrix
        
            // Matrix Market
            // std::cout << "A is " << A << std::endl;
        
            // Maple
        if ( (A.rowdim() < 100) && (A.coldim() < 100) )
            A.write(std::cout << "Pretty A is ", Tag::FileFormat::Maple) << std::endl;
    }
    
        // Vectors
    ZVector X(ZZ, A.coldim()),B(ZZ, A.rowdim());
    
    if (createB) {
        std::cerr << "Creating a random {-1,1} vector " << std::endl;
        srand48( BaseTimer::seed() );
        for(ZVector::iterator it=B.begin();
            it != B.end(); ++it)
            if (drand48() <0.5)
                *it = -1;
            else
                *it = 1;
    } else {
        for(ZVector::iterator it=B.begin();
            it != B.end(); ++it)
            invect >> *it;
    }
    
    {
            // Print RHS
        
        std::cout << "B is [";
        for(auto it:B) ZZ.write(std::cout, it) << " ";
        std::cout << "]" << std::endl;
    }
    
    std::cout << "B is " << B.size() << "x1" << std::endl;
    
    Timer chrono; 

        // BlasElimination
    Method::SparseElimination M;
    M.singularity = Singularity::NonSingular;

        //====================================================
        // BEGIN Replacement solve with fixed prime
    Method::Dixon m(M);
    typedef Givaro::Modular<double> Field;
		// 0.7213475205 is an upper approximation of 1/(2log(2))
    size_t bitsize((size_t)( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)));
    Givaro::Integer randomPrime( *(PrimeIterator<>(bitsize)) );

    FixPrime fixedprime( randomPrime );
    DixonSolver<Ints, Field, FixPrime, Method::SparseElimination> rsolve(A.field(), fixedprime);
    std::cout << "Using: " << *fixedprime << " as the fixed p-adic." << std::endl;


    chrono.start();
    rsolve.solve(X, d, A, B);

        // END Replacement solve with fixed prime
        //====================================================
//     solve (X, d, A, B, M);
    chrono.stop();

    std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;
    
    {
            // Solution size 

        std::cout<<"Reduced solution: \n";
        size_t maxbits=0;
        for (size_t i=0;i<A.coldim();++i){
            maxbits=(maxbits > X[i].bitsize() ? maxbits: X[i].bitsize());
        }
        std::cout<<" numerators of size   "<<maxbits<<" bits" << std::endl
                 <<" denominators hold over "<<d.bitsize()<<" bits\n";	
    }
    
    
    {
			// Check Solution

        VectorDomain<Ints> VD(ZZ);
        MatrixDomain<Ints> MD(ZZ);
        ZVector LHS(ZZ, A.rowdim()), RHS(ZZ, B);
            // check that Ax = d.b
        MD.vectorMul(LHS, A, X);
        VD.mulin(RHS, d);
        if (VD.areEqual(LHS, RHS))
            std::cout << "Ax=b : Yes" << std::endl;
        else
            std::cout << "Ax=b : No" << std::endl;
    }
    
    {
            // Print Solution
        
        std::cout << "(SparseElimination) Solution is [";
        for(auto it:X) ZZ.write(std::cout, it) << " ";
        std::cout << "] / ";
        ZZ.write(std::cout, d)<< std::endl;		
    }
    
    return 0;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
