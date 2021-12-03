/* tests/test-local-smith.C
 * Copyright (C) LinBox
 *
 * Written by J-G. Dumas
 *
 * --------------------------------------------------------
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
 *
 */


/*! @file  tests/test-smith-form-local.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */ 

#ifdef LINBOX_DEBUG
# ifndef LINBOX_LOCAL_SMITH_OUTPUT_
# define LINBOX_LOCAL_SMITH_OUTPUT_
# endif
#endif

#include <linbox/linbox-config.h>
#include <linbox/util/contracts.h>
#include <linbox/matrix/sparse-matrix.h>
#include <givaro/modular.h>
#include <givaro/modular-ruint.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>
#include <linbox/algorithms/smith-form-sparseelim-poweroftwo.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/smith-form.h>
#include <recint/rint.h>
#include <recint/ruint.h>
#include <linbox/util/commentator.h>
#include <iostream>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/util/timer.h>

using namespace LinBox;

/** @brief Test 1: Invariant factors of random sparse matrices.
 *
 * Construct a random matrix with random RPM 
 *           and random local Smith form.
 * Then uses sparse elimination modulo power of prime.
 *
 * Return true on success and false on failure
 */

template<typename Base>
bool check_ranks(const std::vector<std::pair<Base,size_t> >& local,
                 const std::map<int, size_t>& map_values,
                 const Base& p) {
    bool pass(true);

    auto li(local.begin());
    auto mi(map_values.begin());
    for (; li != local.end(); ++li, ++mi) {
        if ( (li->second != mi->second) ||
             (li->first != Givaro::power(p,mi->first)) ) {
            std::cerr << "*** ERROR *** (" <<
                li->first << ',' << li->second << ')' 
                      << " ---->---- (" <<
                mi->first << ',' << mi->second << ')'
                      << std::endl;
            return pass = false;
        }
    }
	commentator().stop (MSG_DONE, nullptr, "SELS");
    return pass;
}


template<typename Base, typename SparseMat>
bool sparse_local_smith(SparseMat& B,
                        size_t R, size_t M, size_t N,
                        const Base& p, int exp,
                        const std::map<int, size_t>& map_values) {
    typedef typename std::remove_reference<decltype(B.field())>::type ModRing;
    PowerGaussDomain< ModRing > PGD( B.field() );
    std::vector<std::pair<Base,size_t> > local;
    Permutation<ModRing> Q(B.field(),B.coldim());
    PGD(local, B, Q, Givaro::power(p,exp), p, PRESERVE_UPPER_MATRIX);

	std::ostream &report = commentator().report ();
    report << "Computed " << p << "-local smith form: (" ;
    for (auto ip = local.begin(); ip != local.end(); ++ip)
        report << '[' << ip->first << ',' << ip->second << "] ";
    report << ")" << std::endl;

    bool pass = check_ranks(local,map_values,p);

	commentator().start ("Check local smith rank", "SELSR");
    
        // Map the resulting PRESERVED upper matrix, mod p
        // check that the rank is preserved, mod p
    ModRing F(p); 
    SparseMatrix<ModRing> BF(F,M,N);
    MatrixHom::map(BF,B);
    size_t rr;
    LinBox::rank(rr,BF, Method::SparseElimination() );

    report << "Residue rank: " << rr << std::endl;
	commentator().stop (MSG_DONE, nullptr, "SELSR");

    return pass && (rr == R);
}


template<typename Base, typename SparseMat>
bool sparse_local_smith_poweroftwo(SparseMat& B,
                        size_t R, size_t M, size_t N,
                        const Base& p, int exp,
                        const std::map<int, size_t>& map_values) {
    LinBox::PowerGaussDomainPowerOfTwo< Base > PGD;
    LinBox::GF2 F2;
    Permutation<GF2> Q(F2,B.coldim());
    std::vector<std::pair<Base,size_t> > local;
    PGD(local, B, Q, exp, PRESERVE_UPPER_MATRIX);

	std::ostream &report = commentator().report();
    report << "Computed 2-local smith form : (" ;
    for (auto ip = local.begin(); ip != local.end(); ++ip)
        report << '[' << ip->first << ',' << ip->second << "] ";
    report << ")" << std::endl;

	bool pass = check_ranks(local,map_values,p);

    commentator().start ("Check binary local smith rank", "SEBLSR");
    
        // Map the resulting PRESERVED upper matrix, mod p
        // check that the rank is preserved, mod p
    GaussDomain<GF2>::Matrix BF(F2,M,N);
    for(auto iter=B.IndexedBegin(); iter != B.IndexedEnd(); ++iter) {
            bool val; Givaro::Caster(val, iter.value());
            if (val) BF.setEntry(iter.rowIndex(), iter.colIndex(), val);
    }
    size_t rr;
    LinBox::rank(rr,BF, Method::SparseElimination() );

    report << "Rank mod 2: " << rr << std::endl;
	commentator().stop (MSG_DONE, nullptr, "SEBLSR");

    return pass && (rr == R);
}

template<typename Base, typename Compute = Base>
bool test_sparse_local_smith(size_t seed, size_t R, size_t M, size_t N,
                             const Base& p, int exp,
                             double density) {
	std::ostream &report = commentator().report ();
	commentator().start ("Testing sparse elimination local smith", "SELS");
	/*
    commentator().indent(std::cerr);
    std::cerr << "test_sparse_local_smith(" << density << ", r" << R << '[' << M << 'x' << N << "]) mod " << p << '^' << exp << std::endl;
	*/
    commentator().indent(report);
    report << "test_sparse_local_smith(" << density << ", r" << R << '[' << M << 'x' << N << "]) mod " << p << '^' << exp << std::endl;

    ASSERT(R<=M && R<=N);

	commentator().start ("Random sparse matrix with random RPM", "RSMRRPM");
    typedef Givaro::Modular<Base,Compute> ModRing;
    ModRing F(Givaro::power(p,exp));

    std::vector<int> v(R);
    std::srand(seed);
    std::generate(v.begin(), v.end(),
                  [&exp]()->int { return std::rand()% exp; });

    std::map<int, size_t> map_values;
    std::for_each(v.begin(), v.end(),
                  [&](int value) { map_values[value]++; } );

    report << "Generated Smith form exponents:\n";
    for(auto & it : map_values)
        report << it.first << ": " << it.second << " times\n";
    report << std::endl;
    
    std::vector<int> indices(R);
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 gen; gen.seed(seed);
    std::shuffle(indices.begin(), indices.end(), gen);

    std::uniform_real_distribution<> nonzero(0., 1.);

    typename ModRing::RandIter G(F,seed);
    typename ModRing::Element pp;
    DenseMatrix<ModRing> A(F,M,N);
    DenseMatrix<ModRing> L(F,M,N);
    for (size_t i=0; i<R; ++i){
        const size_t j = indices[i];
        F.init(pp,Givaro::power(p,v[i]));
        L.setEntry(i,j,pp);
        for (size_t l=i+1; l < M; ++l) {
            if (nonzero(gen)<density) {
                G.random (L.refEntry(l,j));
                F.mulin(L.refEntry(l,j),pp);
            }
        }
    }
#ifdef LINBOX_DEBUG
        L.write(std::cerr<<"L:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif
    DenseMatrix<ModRing> U(F,N,N);
    for (size_t i=0 ; i<std::min(M,N); ++i){
        U.setEntry(i,i,F.one);
        for (size_t j= i+1; j<N ;++j)
            if (nonzero(gen)<density) G.random ( U.refEntry(i,j) );
    }
#ifdef LINBOX_DEBUG
        U.write(std::cerr<<"U:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif

    FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, F.one, L.getPointer(), N, U.getPointer(), N, F.zero, A.getPointer(), N);

#ifdef LINBOX_DEBUG
        A.write(std::cerr<<"A:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif
	commentator().stop (MSG_DONE, nullptr, "RSMRRPM");

    if (p == Base(2)) {
        typedef Givaro::ZRing<Compute> CompRing;
        CompRing RR;
        SparseMatrix<CompRing, SparseMatrixFormat::SparseSeq > B (RR,M,N);
        MatrixHom::map(B, A);
        return sparse_local_smith_poweroftwo(B,R,M,N,Compute(p),exp,map_values);
    } else {
        typedef Givaro::Modular<Compute> CompRing;
        CompRing FF(Givaro::power(p,exp));
        SparseMatrix<CompRing, SparseMatrixFormat::SparseSeq > B (FF,M,N);
        MatrixHom::map(B, A);
        return sparse_local_smith(B,R,M,N,Compute(p),exp,map_values);
    }
}

template<size_t K>
bool ruint_test(size_t seed, size_t R, size_t M, size_t N,
                const Integer& q, int exp, double density) {
        // This is just to generate the random test via fgemm_mp
        // as rns needs minimal size
    int ne = ((1<<(K-1))/(q.bitsize()-1))+1; ne = (ne>exp?ne:exp);

    return test_sparse_local_smith<RecInt::ruint<K>,RecInt::ruint<K+1>>(
        seed,R,M,N,RecInt::ruint<K>(q),ne,density);
}


int main (int argc, char **argv) {
    bool pass0(true);
    static int64_t m = 35;
    static int64_t n = 38;
    static int64_t r = 23;
    static Integer  q = 5;
    static int32_t e = 12;
    static double d = 0.3;
    static size_t iter = 2;
    static int rseed = (int)time(NULL);

	static Argument args[] = {
		{ 'm', "-m M", "Set dimension of test matrices to MxN.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set dimension of test matrices to MxN.", TYPE_INT,     &n },
		{ 'r', "-r R", "Set rank of test matrices to R.", TYPE_INT,     &r },
		{ 'i', "-i i", "Set # repetitions.", TYPE_INT,     &iter },
		{ 'q', "-q Q", "Operate over the ring Z/q^eZ.", TYPE_INTEGER, &q },
		{ 'e', "-e e", "Operate over the ring Z/q^eZ.", TYPE_INT, &e },
		{ 'd', "-d d", "Density of random sparse matrices.", TYPE_DOUBLE, &d },
        { 's', "-s S", "Random generator seed.", TYPE_INT,     &rseed }	,
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	r = std::min(r,m);
	r = std::min(r,n);

    commentator().setMaxDetailLevel (-1);
    commentator().setMaxDepth (-1);
    //commentator().setReportStream (std::cerr);
	std::ostream &report = commentator().report ();

	FFLAS::writeCommandString(report << argv[0] << ' ', args) << std::endl;
    { // sparseelim
        pass0 &= test_sparse_local_smith(rseed,r,m,n,Givaro::Integer(2),e,d);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,Givaro::Integer(2),e,0.01); // to check code against issue #286 (first row is zero)
        pass0 &= test_sparse_local_smith(rseed,r,m,n,Givaro::Integer(q),e,d/2.);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,uint64_t(2),e,d);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,uint64_t(2),e,0.01); // to check code against issue #286 (first row is zero)
        pass0 &= test_sparse_local_smith(rseed,r,m,n,uint64_t(q),e,d/2.);
        pass0 &= ruint_test<6>(rseed,r,m,n,2,e,d);
        pass0 &= ruint_test<6>(rseed,r,m,n,q,e,d/2.);
        pass0 &= ruint_test<7>(rseed,r,m,n,2,e,d);
        pass0 &= ruint_test<7>(rseed,r,m,n,q,e,d/2.);
        pass0 &= ruint_test<8>(rseed,r,m,n,2,e,d);
        pass0 &= ruint_test<8>(rseed,r,m,n,q,e,d/2.);

        Givaro::IntPrimeDom IPD; IPD.seeding(rseed);
        std::mt19937 gen; gen.seed(rseed);
        for(size_t i(0); i<iter; ++i) {
            int64_t im((gen()%m)+1), in((gen()%n)+1);
            int64_t ir(std::min(int64_t((gen()%r)+1),im));
            ir = std::min(ir,in);
            int32_t ie((gen()%e)+1);
            double id(d>0.5?0.5:d);
            Integer iq;
            Integer::random(iq,std::max(q,Integer(131071)));
            IPD.nextprimein(iq);
            pass0 &= test_sparse_local_smith(rseed,ir,im,in,iq,ie,id);
        }
    }
    
	return pass0 ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
