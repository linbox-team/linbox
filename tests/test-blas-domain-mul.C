/*
 * Copyright (C) 2020  Pascal Giorgi
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

#include <linbox/linbox-config.h>
#include <iostream>
#include <string>

#include "linbox/integer.h"

#include "givaro/modular-balanced.h"
#include "givaro/modular.h"
#include "givaro/zring.h"

#include "linbox/matrix/dense-matrix.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/randiter/random-prime.h"

#include "test-matrix-utils.h"

using namespace LinBox;
using namespace std;

template<typename Matrix1, typename Matrix2>
void transpose(Matrix1& tA, const Matrix2& A){
    for(size_t i =0;i<A.rowdim();i++)
        for(size_t j=0;j<A.coldim();j++)
            tA.setEntry(j,i,A.getEntry(i,j));
}

template <class Field>
static bool testMulAdd (const Field& F, size_t n, uint64_t b, size_t iter, size_t seed)
{

	typedef typename Field::Element     Element;
	typedef typename Field::RandIter   RandIter;
	typedef BlasMatrix<Field>          Matrix;
    typedef BlasVector<Field>          Vector;
    typedef typename Matrix::subMatrixType subMatrix;
    typedef typename Vector::subVectorType subVector;
    typedef typename Matrix::constSubMatrixType constSubMatrix;
    typedef typename Vector::constSubVectorType constSubVector;
    
    ostream &report = commentator().report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);

    std::ostringstream msg_field;
    F.write(msg_field);
    std::string msg="Testing square muladd with "+msg_field.str();
    commentator().start (msg.c_str(),"testMulAdd");

    Givaro::Integer samplesize(1); samplesize<<=b;
    RandIter G(F,seed,samplesize);

	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (size_t k=0;k<iter; ++k) {

		commentator().progress(k);

		
		Element alpha, beta,malpha,a,b;
        F.init(alpha);F.init(beta);F.init(malpha);
        F.init(a); F.init(b);

        Matrix A(F,n,n), B(F,n,n),C(F, n,n),D(F,n,n),T(F,n,n),R(F,n,n);
		Vector x(F,n),y(F,n),z(F,n),t(F,n);        

		// Create 3 random n*n matrices
		A.random(G); B.random(G);	C.random(G);

		// Create 2 random vectors
        x.random(G); y.random(G);

		// create 2 random element
		G.random(a);
		G.random(b);


        std::vector<std::pair<Element,Element>> scal={{a,b},{F.one,F.one},{F.one,F.mOne},{F.mOne,F.one},{F.mOne,F.mOne},{F.zero,F.one},{F.one,F.zero},{F.zero,F.mOne},{F.mOne,F.zero}};
        
        for(std::pair<Element,Element> ab:scal){
            F.assign(alpha,ab.first);
            F.assign(beta ,ab.second);
            F.neg(malpha,alpha);
            
            /***** CHECK 1 *****/
            // compute D = -alpha.(A*C+B*C) + alpha.(A+B)*C = 0
            BMD.mul(D,A,C);
            //cout<<"D="<<D<<endl;        
            BMD.mul(T,B,C);
            //cout<<"T="<<T<<endl;
            BMD.addin(D,T);
            //cout<<"D="<<D<<endl;
            BMD.add(T,A,B);
            //cout<<"T="<<T<<endl;
            BMD.muladd(R,malpha,D,alpha,T,C);
            if (!BMD.isZero(R)) {
                report<<"Check 1 failed (alpha=";F.write(report,alpha)<<")"<<endl;
                report<<"computing D = -alpha.(A*C+B*C) + alpha.(A+B)*C = 0"<<endl;
                report<<"A="<<A<<endl;
                report<<"B="<<B<<endl;
                report<<"C="<<C<<endl;
                report<<"D="<<R<<endl;
                ret=false;
            }

            /***** CHECK 2 *****/
            // compute z = beta.y + alpha.A*x
            // * using BlasMatrixDomain
            BMD.muladd(z,beta,y,alpha,A,x);
            // * by hand
            A.apply(t, x);
            for (size_t i=0;i<n;++i){
                F.mulin(t[i],alpha);
                F.axpyin(t[i],beta,y[i]);
            }        
            if (!(t==z)){
                report<<"Check 2 failed (alpha=";F.write(report,alpha)<<", beta="; F.write(report,beta)<<")"<<endl;
                report <<"computing z = beta.y + alpha.A*x "<<endl;
                report<<"A="<<A<<endl;
                report<<"x="<<x<<endl;
                report<<"y="<<y<<endl;
                report<<"z="<<z<<endl;
                report<<"t="<<t<<endl;
                ret=false;
            }

            /***** CHECK 3 *****/
            // using submatrices and subvectors
            size_t k1=(n+1)/2; 
            //size_t k2=n-k1;
        
            subMatrix D1(D,0,0,k1,k1);
            subMatrix T1(T,0,0,k1,k1); 
            subMatrix R1(R,0,0,k1,k1);
            constSubMatrix A1(A,0,0,k1,k1); 
            constSubMatrix B1(B,0,0,k1,k1); 
            constSubMatrix C1(C,0,0,k1,k1);
            subVector t1(t,0,1,k1);
            subVector z1(z,0,1,k1);
            constSubVector x1(x,0,1,k1);
            constSubVector y1(y,0,1,k1);

            // compute D1 = -alpha.(A1*C1+B1*C1) + alpha.(A1+B1)*C1 = 0        
            // cout<<"alpha="<<alpha<<endl;
            // cout<<"A1="<<A1<<endl;
            // cout<<"B1="<<B1<<endl;
            // cout<<"C1="<<C1<<endl;
        
            BMD.mul(D1,A1,C1);
            //cout<<"D1="<<D1<<endl;
            BMD.mul(T1,B1,C1);
            //cout<<"T1="<<T1<<endl;
            BMD.addin(D1,T1);
            //cout<<"D1="<<D1<<endl;
            BMD.add(T1,A1,B1);
            //cout<<"T1="<<T1<<endl;
            BMD.muladd(R1,malpha,D1,alpha,T1,C1);
            if (!BMD.isZero(R1)) {
                report<<"Check 3 failed (alpha=";F.write(report,alpha)<<")"<<endl;
                report<<"computing submatrix [1.."<<k1<<",1.."<<k1<<"of D = -alpha.(A*C+B*C) + alpha.(A+B)*C = 0"<<endl;
                report<<"A="<<A<<endl;
                report<<"B="<<B<<endl;
                report<<"C="<<C<<endl;
                report<<"D="<<R<<endl;
                ret=false; 
            }

            /***** CHECK 3 *****/
            // compute z1 = beta.y1 + alpha.A1*x1     
            // * using BlasMatrixDomain
            BMD.muladd(z1,beta,y1,alpha,A1,x1);
            // * by hand
            A1.apply(t1, x1);
            for (size_t i=0;i<k1;++i){
                F.mulin(t1[i],alpha);
                F.axpyin(t1[i],beta,y1[i]);
            }        
            if (!(t1==z1)){
                report<<"Check 4 failed (alpha=";F.write(report,alpha)<<", beta="; F.write(report,beta)<<")"<<endl;
                ret=false;
            }
        }
    }

	commentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAdd");
	return ret;
}


// tests MulAdd for various shapes and values of transposition.
template <class Field>
static bool testMulAddShapeTrans (const Field &F, size_t m, size_t n, size_t k, uint64_t b, size_t iter, size_t seed)
{
    std::ostringstream msg_field;
    F.write(msg_field);
    std::string msg="Testing rectangular muladd using transposition with "+msg_field.str();
    commentator().start (msg.c_str(),"testMulAddShapeTrans");

    
	typedef typename Field::Element Element;
	typedef BlasMatrix<Field>                Matrix;
    //typedef typename Matrix::subMatrixType subMatrix;
	typedef TransposedBlasMatrix<Matrix> TransposedMatrix ;
	typedef typename Field::RandIter RandIter ;
    typedef typename Field::NonZeroRandIter NonZeroRandIter ;

    Givaro::Integer samplesize(1); samplesize<<=b;
    RandIter G(F,seed,samplesize);
    NonZeroRandIter Gnz(G);

	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	// input matrix
	Matrix A(F, m,k);
	Matrix B(F, k,n);
	Matrix C(F, m,n);
	// result matrix
	Matrix D(F, m,n);
	Matrix E(F, m,n);

    for (size_t h=0;h<iter; ++h) {

		commentator().progress(h);

        // random A,B,C
        A.random(G);
        B.random(G);
        C.random(G);

        // hard tranpose A,B
        Matrix A1 (F, k,m) ;
        transpose(A1,A) ;
        Matrix B1 (F, n,k) ;
        transpose(B1,B) ;
        TransposedMatrix tA(A1); // t(tA)=A
        TransposedMatrix tB(B1); // t(tB)=B

        // random alpha, beta
        Element alpha ;    
        Element beta ;
        G.random(alpha);
        G.random(beta);
    
        // witness (D= beta.C+ alpha.A.B 
        BMD.muladd(D,beta,C,alpha,A,B);

        // E= beta.C+ alpha.A.B 
        BMD.muladd(E,beta,C,alpha,A,B);
        if (!BMD.areEqual(E,D)) {
            ret = false ;
            commentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (noTrans, noTrans) *** " << std::endl;
        }

        BMD.muladd(E,beta,C,alpha,A,tB);
        if (!BMD.areEqual(E,D))  {
            ret = false ;
            commentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (noTrans, Trans) *** " << std::endl;
        }

        BMD.muladd(E,beta,C,alpha,tA,B);
        if (!BMD.areEqual(E,D)) {
            ret = false ;
            commentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (Trans, noTrans) *** " << std::endl;
        }

        BMD.muladd(E,beta,C,alpha,tA,tB);
        if (!BMD.areEqual(E,D)) {
            ret = false ;
            commentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (Trans, Trans) *** " << std::endl;
        }
    }

	commentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAddShapeTrans");
	return ret ;
}

// tests MulAdd for various shapes and values of transposition.
template<class Field, bool LeftSide, bool UnitDiag>
static bool testTriangMulShapeTrans (const Field &F, size_t m, size_t n, uint64_t b, size_t iter, size_t seed)
{
    std::ostringstream msg_field;
    F.write(msg_field);
    std::string msg="Testing triangular muladd for shapes and transposition with "+msg_field.str();
    commentator().start (msg.c_str(),"testTriangMulShapeTrans");

	typedef BlasMatrix<Field>                                   Matrix ;
	typedef TriangularBlasMatrix<Matrix>               TriangularMatrix ;
	typedef TransposedBlasMatrix<TriangularMatrix > TransposedTriangular ;
	typedef typename Field::RandIter                            RandIter ;

    Givaro::Integer samplesize(1); samplesize<<=b;
    RandIter G(F,seed,samplesize);


	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	size_t k =(int) (LeftSide?m:n) ;
	// input matrix
	Matrix A(F, k,k); // A = L+U-I. Either L or U is unit.
	Matrix B(F, m,n);
	// result matrix
	Matrix D(F, m,n);
	Matrix E(F, m,n);

	for (size_t h=0;h<iter; ++h) {

		commentator().progress(h);

    
        // random A,B
        A.random(G);
        B.random(G);

        // hard tranpose A,B
        Matrix A1 (F,k,k) ;
        transpose(A1,A) ;

        /*  test (L+U-I) B+B = LB+UB */
        if (LeftSide)
            BMD.muladd(D,F.one,B,F.one,A,B);
        else
            BMD.muladd(D,F.one,B,F.one,B,A);

        /****  DIRECT ****/
        {
            /*  L */
            TriangularMatrix L (A, Tag::Shape::Lower, (UnitDiag?Tag::Diag::Unit:Tag::Diag::NonUnit));

            /*  U */
            TriangularMatrix U (A, Tag::Shape::Upper,(!UnitDiag?Tag::Diag::Unit:Tag::Diag::NonUnit));

            /*  make product */
            E = B ;
            // Matrix G(m,n);
            // G = E ;
            Matrix G((const Matrix&)E); //!@warning on n'oublie pas l'esperluette !!!
            if(LeftSide) {
                BMD.mulin_right(L,E) ; // B <- AB
                BMD.mulin_right(U,G) ;
            }
            else {
                BMD.mulin_left(G,L) ;  // B <- BA
                BMD.mulin_left(E,U) ;
            }
            BMD.addin(E,G);

            /*  check equality */
            if (!BMD.areEqual(E,D)) {
                ret = false ;
                commentator().report() << " *** BMD ERROR (" << (LeftSide?"left":"right") << ',' << (UnitDiag?" L":" U") << " is unit) *** " << std::endl;
            }
            else {
                commentator().report() << " direct triangular multiplication ok." << std::endl;
            }
        }
        /****  Transpose ****/
        {
            /*  L */
            TriangularMatrix L1 (A1, Tag::Shape::Lower,(UnitDiag?Tag::Diag::Unit:Tag::Diag::NonUnit));

            /*  U */
            TriangularMatrix U1 (A1, Tag::Shape::Upper,(!UnitDiag?Tag::Diag::Unit:Tag::Diag::NonUnit));

            TransposedTriangular L(L1);
            TransposedTriangular U(U1);
            /*  make product */
            E = B ;
            // Matrix G(m,n);
            // G = E ;
            Matrix G((const Matrix&)E); //!@warning on n'oublie pas l'esperluette !!!
            if(LeftSide) {
                BMD.mulin_right(L,E) ; // B <- AB
                BMD.mulin_right(U,G) ;
            }
            else {
                BMD.mulin_left(G,L) ;  // B <- BA
                BMD.mulin_left(E,U) ;
            }
            BMD.addin(E,G);

            /*  check equality */
            if (!BMD.areEqual(E,D)) {
                ret = false ;
                commentator().report() << " *** BMD ERROR Transpose (" << (LeftSide?"left":"right") << ',' << (UnitDiag?" L":" U") << " is unit) *** " << std::endl;
            }
            else {
                commentator().report() << " transposed triangular multiplication ok." << std::endl;
            }
        }
    }
	commentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAddShapeTrans");
	return ret ;
}


template <class Field>
static bool testMulPermutation (const Field& F, size_t n, uint64_t b, size_t iter, size_t seed)
{

	//typedef typename Field::Element     Element;
	typedef typename Field::RandIter   RandIter;
	typedef BlasMatrix<Field>          Matrix;
    typedef BlasPermutation<size_t>    Permutation;
    typedef TransposedBlasMatrix<BlasPermutation<size_t>>    TransposedPermutation;
    // typedef typename Matrix::subMatrixType subMatrix;
    // typedef typename Matrix::constSubMatrixType constSubMatrix;
    
    ostream &report = commentator().report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);

    std::ostringstream msg_field;
    F.write(msg_field);
    std::string msg="Testing square mul on blasPermutation with "+msg_field.str();
    commentator().start (msg.c_str(),"testMulPerm");

    Givaro::Integer samplesize(1); samplesize<<=b;
    RandIter G(F,seed,samplesize);

	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);


    Matrix A(F,n,n), B(F,n,n),C(F, n,n),D(F,n,n),T(F,n,n),R(F,n,n);

    for (size_t k=0;k<iter; ++k) {

		commentator().progress(k);

        // Create 3 random n*n matrices
        A.random(G); B.random(G); 
    
        // Create 2 random permutation
        std::vector<size_t> x(n),y(n);
        RandomPermutation(x.data(),n);
        RandomPermutation(y.data(),n);    
        Permutation X(x),Y(y);
        TransposedPermutation Xt(X);
	    
        // compute C=A*B
        BMD.mul(C,A,B);

        // compute R = (AX)(X^TB)
        BMD.mul(D,A,X);
        BMD.mul(T,Xt,B);
        BMD.mul(R,D,T);

        /*  check equality C=R*/
        if (!BMD.areEqual(C,R)) {
            ret = false ;
            report << " *** BMD ERROR mul with permutation failed (AX)(X^TB)" << std::endl;
            report<< "A="<<A<<endl;
            report<< "X=";X.write(report,false)<<endl;
            report<< "AX="<<D<<endl;
            report<< "B="<<B<<endl;
            report<< "XtB="<<T<<endl;
            report << "C="<<C<<endl;
            report << "R="<<R<<endl;
        }
        else {
            report << "  multiplication with permutation ok." << std::endl;
        }

        // compute R = (AX^T)(XB)
        BMD.mul(D,A,Xt);
        BMD.mul(T,X,B);
        BMD.mul(R,D,T);
        /*  check equality C=R*/
        if (!BMD.areEqual(C,R)) {
            ret = false ;
            report << " *** BMD ERROR mul with permutation failed (AX^T)(XB)" << std::endl;
            report<< "A="<<A<<endl;
            report<< "X=";X.write(report,false)<<endl;
            report<< "AXt="<<D<<endl;
            report<< "B="<<B<<endl;
            report<< "XB="<<T<<endl;
            report << "C="<<C<<endl;
            report << "R="<<R<<endl;
        }
        else {
            report << "  multiplication with permutation ok." << std::endl;
        }

        
    }
    
    commentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulPerm");
    return ret;
}
    
// returns true if ok, false if not.
template<class Field>
int launch_tests(const Field & F, size_t n, uint64_t b, size_t iter, size_t seed)
{
    ostream &report = commentator().report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);
	bool pass = true ;
	if (!testMulAdd (F,n,b,iter,seed))                     pass=false;

    size_t m = n+n/2 ; size_t k = 2*n+1 ;
    if (!testMulAddShapeTrans (F,n,m,k,b,iter,seed))       pass=false;
    if (!testMulAddShapeTrans (F,n,k,m,b,iter,seed))       pass=false;
    if (!testMulAddShapeTrans (F,m,n,k,b,iter,seed))       pass=false;
    if (!testMulAddShapeTrans (F,m,k,n,b,iter,seed))       pass=false;
    if (!testMulAddShapeTrans (F,k,n,m,b,iter,seed))       pass=false;
    if (!testMulAddShapeTrans (F,k,m,n,b,iter,seed))       pass=false;
    if (!testTriangMulShapeTrans<Field,true,true>   (F,m,n,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,true>   (F,n,m,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,true>  (F,m,n,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,true>  (F,n,m,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,false>  (F,m,n,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,false>  (F,n,m,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,false> (F,m,n,b,iter,seed))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,false> (F,n,m,b,iter,seed))     pass=false;

    if (!testMulPermutation (F,n,b,iter,seed))                     pass=false;
	if (not pass) F.write(report) << endl;
	return pass ;
}


int main(int argc, char **argv)
{

	static size_t  n = 33;    // matrix dimension
	static integer q = 65537; // field size
	static int iterations = 1; //nbr iterartion
    static long seed = time(NULL);
    static Argument args[] = {
                              { 'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT,     &n },
                              { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]",  TYPE_INTEGER, &q },
                              { 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
                              { 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
                              END_OF_ARGUMENTS
    };
    
	parseArguments (argc, argv, args);

    size_t bits=200;
    PrimeIterator<IteratorCategories::HeuristicTag> Rd(bits,seed);
    integer p= *Rd;    
    
	bool pass = true;

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	commentator().start("BlasMatrixDomain Multiplication test suite", "BlasMatrixDomain");

	Givaro::Modular<double>  F1(q);
    Givaro::Modular<integer> F2(p);			
    Givaro::ZRing<integer>   Z;

	pass &= launch_tests(F1,n,q.bitsize(),iterations,seed);
    pass &= launch_tests(F2,n,p.bitsize(),iterations,seed);
    pass &= launch_tests(Z,n,bits,iterations,seed);
    
	commentator().stop(MSG_STATUS (pass), (const char *) 0,"BlasMatrixDomain");
	return pass ? 0 : -1;
}




// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
