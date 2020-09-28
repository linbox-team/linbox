/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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




#include <iostream>
#include <iomanip>
using namespace std;

#include <linbox/ring/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
#include <givaro/zring.h>
#include <recint/rint.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>
#include <linbox/algorithms/polynomial-matrix/matpoly-mult-naive.h>
#include <linbox/algorithms/polynomial-matrix/matpoly-mult-kara.h>

using namespace LinBox;

template<typename Rand, typename Matrix>
void randomMatPol(Rand& r,  Matrix& A){
	for(size_t i=0;i<A.size()*A.rowdim()*A.coldim();i++)
		r.random(A.getPointer()[i]);
}


template<typename T>
void test(T& x){
    x.ref(0,0)++;
}

template<typename PolMatMulDomain>
bool checkMatPolMul(const PolMatMulDomain& PMMD, size_t m,size_t n, size_t d, long seed, string algo){
    
    commentator().start (string("Testing polynomial matrix Multiplication with ")+algo, "testMatpolyMul", 1);
    bool finalok=true;
    bool pass= true;
  
    ostream& report = LinBox::commentator().report();
    //ostream& report =std::cout;

    typedef typename PolMatMulDomain::Field Field;
    Field F(PMMD.field());
    typename Field::RandIter G(F,seed);
    typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP;
    typedef PolynomialMatrix<Field, PMType::matfirst> PMatrix;
    typedef PolynomialMatrix<Field, PMType::matrowfirst> PMatrixP;

    
    MatrixP A1(F,m,n,d),B1(F,n,n,d), C1(F,m,n,2*d-1);
    PMatrix A2(F,m,n,d),B2(F,n,n,d), C2(F,m,n,2*d-1);
    PMatrixP A3(F,m,n,d),B3(F,n,n,d), C3(F,m,n,2*d-1);

    randomMatPol(G,A1);
    randomMatPol(G,B1);
    A2.copy(A1);  
    B2.copy(B1);
    A3.copy(A1);  
    B3.copy(B1);
    
    //A2[0].copy(A2[1]);
    //A1[0].copy(A1[1]);
    
    // std::cout<<"Problem test-polynomial-matrix: lvalue and rvalue -> todo\n";
    // std::cout<<A1<<std::endl;
    // ((A1[0]).refEntry(0,0))++;
    // std::cout<<A1<<std::endl;
    // std::cout<<A2<<std::endl;
    // ((A2[0]).refEntry(0,0))++;
    // //test(A1[0]);
    // //test(A2[0]);
    // std::cout<<A2<<std::endl;
    
    PMMD.mul(C1,A1,B1);
    PMMD.mul(C2,A2,B2);
    PMMD.mul(C3,A3,B3);
    for(size_t i=0;i<m*n;i++)
        for(size_t k=0;k<2*d-1;k++)
            if (!F.areEqual(C1.get(i,k),C2.get(i,k)) || !F.areEqual(C1.get(i,k),C3.get(i,k)))
                pass=false;

    finalok&=pass;
    report << std::endl;
    report<<"   - Polynomial Matrix (plain)    :"; 
    if (pass)
        report<<"OK"<<endl;
    else
        report<<"KO"<<endl;

    if (!pass){
        report<<"A1:="<<A1<<endl;
        report<<"B1:="<<B1<<endl;
        report<<"C1:="<<C1<<endl;
        report<<"---------------"<<endl;
        report<<"A2:="<<A2<<endl;
        report<<"B2:="<<B2<<endl;
        report<<"C2:="<<C2<<endl;
        report<<"---------------"<<endl;
        report<<"A3:="<<A3<<endl;
        report<<"B3:="<<B3<<endl;
        report<<"C3:="<<C3<<endl;
        report<<endl<<endl;
    }
    
    if (d>1){
        size_t d2=d/2;
        auto AA1=A1.at(0,d2-1); auto BB1=B1.at(d-d2,d-1); auto CC1=C1.at(0,2*d2-2);
        auto AA2=A2.at(0,d2-1); auto BB2=B2.at(d-d2,d-1); auto CC2=C2.at(0,2*d2-2);
        auto AA3=A3.at(0,d2-1); auto BB3=B3.at(d-d2,d-1); auto CC3=C3.at(0,2*d2-2);

        PMMD.mul(CC1,AA1,BB1);
        PMMD.mul(CC2,AA2,BB2);
        PMMD.mul(CC3,AA3,BB3);
  
        pass= true;
        for(size_t i=0;i<m*n;i++)
            for(size_t k=0;k<d-1;k++)
                if (!F.areEqual(CC1.get(i,k),CC2.get(i,k))  || !F.areEqual(C1.get(i,k),C3.get(i,k)))
                    pass=false;

        finalok&=pass;
        report<<"   - Polynomial Matrix (view)     :"; 
        if (pass)
            report<<"OK"<<endl;
        else
            report<<"KO"<<endl;

        if (!pass){
            report<<"AA1:="<<AA1<<endl;
            report<<"BB1:="<<BB1<<endl;
            report<<"CC1:="<<CC1<<endl;
            report<<"---------------"<<endl;
            report<<"AA2:="<<AA2<<endl;
            report<<"BB2:="<<BB2<<endl;
            report<<"CC2:="<<CC2<<endl;
            report<<"---------------"<<endl;
            report<<"AA3:="<<AA3<<endl;
            report<<"BB3:="<<BB3<<endl;
            report<<"CC3:="<<CC3<<endl;
            report<<endl<<endl;
        }
        
        A1.resize(d2);  A2.resize(d2); A3.resize(d2);
        B1.resize(d2);  B2.resize(d2); B3.resize(d2);        
        C1.resize((d2<<1)-1);  C2.resize((d2<<1)-1); C3.resize((d2<<1)-1);

        PMMD.mul(C1,A1,B1);
        PMMD.mul(C2,A2,B2);
        PMMD.mul(C3,A3,B3);
        
        pass= true;
        for(size_t i=0;i<m*n;i++)
            for(size_t k=0;k<(d2<<1)-1;k++)
                if (!F.areEqual(C1.get(i,k),C2.get(i,k)) || !F.areEqual(C1.get(i,k),C3.get(i,k)))
                    pass=false;
        
        finalok&=pass;
        
        report<<"   - Polynomial Matrix (resizing) :"; 
        if (pass)
            report<<"OK"<<endl;
        else
            report<<"KO"<<endl;
        
        if (!pass){
            report<<"A1:="<<A1<<endl;
            report<<"B1:="<<B1<<endl;
            report<<"C1:="<<C1<<endl;
            report<<"---------------"<<endl;
            report<<"A2:="<<A2<<endl;
            report<<"B2:="<<B2<<endl;
            report<<"C2:="<<C2<<endl;
            report<<"---------------"<<endl;
            report<<"A3:="<<A3<<endl;
            report<<"B3:="<<B3<<endl;
            report<<"C3:="<<C3<<endl;
            report<<endl<<endl;
        }  
    }

    commentator().stop(MSG_STATUS(pass),(const char *) 0,"testMatpolyMul");
    return pass;
}



int main(int argc, char** argv){
    static integer q= 101;
    static size_t  m = 10; // matrix dimension
    static size_t  d = 50; // polynomial size
    static long    seed = time(NULL);

    static Argument args[] = {
                              { 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
                              { 'q', "-q Q", "Set the prime field charactetistic.", TYPE_INTEGER,     &q },
                              { 'd', "-d D", "Set degree bound of test matrices to D.", TYPE_INT,     &d },
                              { 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
                              END_OF_ARGUMENTS
    };
    parseArguments (argc, argv, args);


    commentator().start ("Testing polynomial matrix", "testMatpoly", 1);
    bool pass=true;

    typedef Givaro::Modular<double> Field;
    //typedef Givaro::Modular<uint32_t> Field;

    Field F(q);
    PolynomialMatrixNaiveMulDomain<Field> PMMD_naive(F);
    PolynomialMatrixKaraDomain<Field>     PMMD_kara(F);
    PolynomialMatrixFFTMulDomain<Field>   PMMD_fft(F);    

    ostream & report = commentator().report();
    
    report<<"Polynomial matrix testing over ";F.write(report)<<std::endl;
    pass &= checkMatPolMul(PMMD_naive,m,m,d,seed, "Naive");
    pass &= checkMatPolMul(PMMD_kara,m,m,d,seed, "Karatsuba");
    pass &= checkMatPolMul(PMMD_fft,m,m,d,seed, "FFT");

    commentator().stop(MSG_STATUS(pass),(const char *) 0,"testMatpoly");
   
    return (pass? 0: -1);
} 



// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
