/*  Copyright (c) 2015 HPAC
 *  Written by 	Alexis Breust <alexis.breust@gmail.fr>
 * 		Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 */

#include <linbox/linbox-config.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <givaro/zring.h>
#include <linbox/vector/light_container.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/permutation.h>
#include <linbox/blackbox/submatrix.h>
#include <linbox/blackbox/diagonal.h>
#include "common_common.hpp"
#include "common_kernel.hpp"


using namespace LinBox;
typedef Givaro::ZRing<Modulus_t> Ring;

typedef SparseMatrix<Ring, SparseMatrixFormat::SparseSeq > SpMatrix;

SpMatrix& Transpose(SpMatrix& AT, const SpMatrix& A) {
    AT.resize(A.coldim(),A.rowdim());
    for(size_t i=0; i<A.rowdim(); ++i) {
        for(size_t jj=0;jj<A[i].size();++jj)
            AT.setEntry(A[i][jj].first,i,A[i][jj].second);
    }
    return AT;
} 

SpMatrix Transpose(const SpMatrix& A) {
    SpMatrix At(A.field(), A.coldim(), A.rowdim());
    return Transpose(At,A);
} 


typedef SpMatrix::Row::value_type Pair_t;

SpMatrix::Row& applyTranspose (SpMatrix::Row& y, const std::vector<size_t>& P, const SpMatrix::Row& x) {
    y.resize(0);
    for (size_t i = 0; i < x.size(); ++i) {
        size_t k = x[i].first;
        y.push_back(Pair_t(P[k],x[i].second));
    }
    return y;
}

struct ComparePair {
    bool operator()(const Pair_t& x, const Pair_t& y) {
        return x.first < y.first;
    }
};

            
    

// Y=PX
SpMatrix& applyP(SpMatrix& Y, const Permutation<Ring>& P, const SpMatrix& X) {
    Y.resize(X.rowdim(),X.coldim());
    for (size_t i = 0; i < Y.rowdim(); ++i) {
        Y[i].resize(0);
        Y[i].insert(Y[i].begin(),X[P[i]].begin(),X[P[i]].end());
    }
    return Y; 
}

// Y=XP
SpMatrix& applyP(SpMatrix& Y, const SpMatrix& X, const std::vector<size_t>& P) {
    Y.resize(X.rowdim(),X.coldim());
    for (size_t i = 0; i < Y.rowdim(); ++i) {
        applyTranspose(Y[i],P,X[i]);
        std::sort(Y[i].begin(),Y[i].end(), ComparePair());
    }
    return Y;
}

int main(int argc, char ** argv) 
{
    std::cerr << "command-line:";
    for(int i = 0; i < argc; ++i) 
       std::cerr << ' ' << argv[i];
    std::cerr << std::endl;
    
// P A Q = [ D C | Z1 ]
//         [ 0 B | Z2 ]
// It might be sufficient to solve for B x = 0 ...

    if (argc != 2) {
        std::cerr << "[DECO] usage: " << argv[0] << " filename" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ifstream input (argv[1]);
    std::string filename=removeExtension(argv[1]);
    Ring ZZ;
    Ring::Element val;
    size_t i,j; char t;
    input >> i >> j >> t;
    std::cerr << "[DECO] reading " << i << 'x' << j << " ..." << std::endl;

    assert(i<j);

    SpMatrix A (ZZ); A.resize(j,i,0);
    while( i!= 0 ) {
        ZZ.read(input >> i >> j,val);
        if (! ZZ.isZero(val)) A.setEntry(j-1,i-1,val);
    }
    std::cerr << "[DECO] Transpose Matrix read."<< std::endl;
    
    Permutation<Ring> P(ZZ, A.rowdim());
    size_t last=A.rowdim()-1;
    size_t diagonalsize=0;

    std::set<size_t> pivots;

    for(size_t ii=0;ii<=last;++ii) {
        if (A[ii].size()==0) P.permute(ii,last--);
        if (A[ii].size()==1) {
            auto search = pivots.find( A[ii][0].first );
            if(search != pivots.end()) {
                if (last>ii) P.permute(ii,last--);
            } else {
                P.permute(diagonalsize++,ii);
                pivots.insert(A[ii][0].first);
            }
        }
        
    }

    std::ofstream matD(filename+".matD.sms"), matC(filename+".matC.sms"), matB(filename+".matB.sms"), matA2(filename+".matZ.sms"), matP(filename+".matP.sms");

    P.write(matP) << std::endl;
    matP.close();
    std::cerr << "[DECO] P " << P.rowdim() << 'x' << P.coldim() << " wrote."<< std::endl;

    std::vector<size_t> Q(A.coldim()); 
    for(size_t ii=0;ii<Q.size();++ii) Q[ii]=ii;
    for(size_t ii=0;ii<diagonalsize;++ii)
    {
        Q[ii] = A[P[ii]][0].first;
        Q[Q[ii]]=ii;
    }
    

// P.write(std::cout << "P:=") << ':' << std::endl;
// std::cout << "Q:=[";for(auto vals:Q) std::cout << vals << ' '; std::cout << ']' << std::endl;

    std::cerr << "[DECO] Permutation found with " << (A.rowdim()-1-last) << " useless columns."<< std::endl;
    std::cerr << "[DECO] Permutation found with " << diagonalsize << " singletons."<< std::endl;

   
// A.write(std::cout, Tag::FileFormat::Guillaume) << std::endl;
// A.write(std::cerr << "A:=", Tag::FileFormat::Maple) << ':' << std::endl;

    SpMatrix Az (ZZ); Az.resize(A.rowdim(),A.coldim(),0);
    SpMatrix At (ZZ); At.resize(A.rowdim(),A.coldim(),0);

//     P.applyRight(Az,A);
//     Q.applyLeft(At,Az);

    applyP(Az,P,A);
// Az.write(std::cerr << "Az:=", Tag::FileFormat::Maple) << ':' << std::endl;
    applyP(At,Az,Q);
    
    std::cerr << "[DECO] Permutations applied."<< std::endl;
// At.write(std::cerr << "At:=", Tag::FileFormat::Maple) << ':' << std::endl;
    
    SpMatrix C(ZZ,At.coldim()-diagonalsize,diagonalsize);
    SpMatrix B(ZZ,At.coldim()-diagonalsize,At.coldim()-diagonalsize);
    for(size_t ii=0;ii<C.rowdim();++ii) {
        const SpMatrix::Row& Atii = At[ii+diagonalsize];
        C[ii].resize(0);
        B[ii].resize(0);
        for(size_t jj=0;jj<Atii.size();++jj) {
            const size_t y = Atii[jj].first;
            if (y < diagonalsize)
                C.setEntry(ii,y,Atii[jj].second);
            else 
                B.setEntry(ii,y-diagonalsize,Atii[jj].second);
        }
    }
    Transpose(C).write(matC, Tag::FileFormat::Guillaume) << std::endl;
    matC.close();
    std::cerr << "[DECO] C " << C.coldim() << 'x' << C.rowdim() << " wrote."<< std::endl;
    Transpose(B).write(matB, Tag::FileFormat::Guillaume) << std::endl;
    matB.close();
    std::cerr << "[DECO] B " << B.coldim() << 'x' << B.rowdim() << " wrote."<< std::endl;
    
    SpMatrix A2(ZZ,At.rowdim()-At.coldim(),At.coldim());
    for(size_t ii=0;ii<A2.rowdim();++ii) {
        SpMatrix::Row& Atii = At[ii+At.coldim()];
        A2[ii].resize(0);
        for(size_t jj=0;jj<Atii.size();++jj)
            A2.setEntry(ii,Atii[jj].first,Atii[jj].second);
    }

    std::cerr << "[DECO] Matrix cut into B, C, Z."<< std::endl;
    Transpose(A2).write(matA2, Tag::FileFormat::Guillaume) << std::endl;
    matA2.close();
    std::cerr << "[DECO] Z " << A2.coldim() << 'x' << A2.rowdim() << " wrote."<< std::endl;
    

//     C.write(std::cerr << "C:=", Tag::FileFormat::Maple) << ':' << std::endl;
//     B.write(std::cerr << "B:=", Tag::FileFormat::Maple) << ':' << std::endl;
//     A2.write(std::cerr << "Z:=", Tag::FileFormat::Maple) << ':' << std::endl;
    
    Diagonal<Ring> D(ZZ,diagonalsize);
    Ring::Element x;
    for(size_t ii=0; ii<diagonalsize;++ii) {
        At.getEntry(x,ii,ii);
        assert( ZZ.areEqual(x,At[ii][0].second) );
        D.setEntry(ii,ii,At[ii][0].second);
//         D.setEntry(ii,ii,At.getEntry(x,ii,ii));
    }
    
    std::cerr << "[DECO] Diagonal build."<< std::endl;
//     D.write(std::cerr) << std::endl;
//     Transpose(C).write(std::cerr << "C:=", Tag::FileFormat::Maple) << ':' << std::endl;
//     Transpose(B).write(std::cerr << "B:=", Tag::FileFormat::Maple) << ':' << std::endl;
//     Transpose(A2).write(std::cerr << "Z:=", Tag::FileFormat::Maple) << ':' << std::endl;

    D.write(matD) << std::endl;
    std::cerr << "[DECO] D " << D.rowdim() << 'x' << D.coldim() << " wrote."<< std::endl;

    

 
    return 0;
}

