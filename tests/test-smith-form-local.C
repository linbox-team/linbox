/* tests/test-local-smith.C
 * Copyright (C) LinBox
 *
 * Written by David Saunders
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

#ifdef DEBUG
# ifndef LINBOX_LOCAL_SMITH_OUTPUT_
# define LINBOX_LOCAL_SMITH_OUTPUT_
# endif
#endif

#include "linbox/linbox-config.h"
#include <linbox/util/contracts.h>
#include <linbox/matrix/sparse-matrix.h>
#include <givaro/modular.h>
#include <givaro/modular-ruint.h>
#include <linbox/algorithms/smith-form-sparseelim-local.h>
#include <linbox/algorithms/smith-form-sparseelim-poweroftwo.h>
#include <linbox/solutions/rank.h>
#include <recint/rint.h>
#include <recint/ruint.h>



#include <functional>

#include "test-common.h"
#include "test-field.h"

#include "linbox/util/commentator.h"
#include "linbox/ring/local-pir-modular.h"
#include "linbox/ring/pir-modular-int32.h"
#include "linbox/ring/local2_32.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/util/timer.h"

using namespace LinBox;

/** @brief Test 1: Invariant factors of random dense matrices.
 *
 * Construct a random matrix which is equivalent to a random diagonal matrix,
 * and check its Smith form.
 *
 * R - PIR over which to perform computations
 * stream - Stream that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class LocalPIR>
typename LocalPIR::Element &
normal(typename LocalPIR::Element &d, LocalPIR & R)
{	typename LocalPIR::Element x; R.init(x, 2);  R.mulin(x, d);
	return R.gcdin(d, x);
}

template <class LocalPIR>
class pplt { // prime power less than
public:
	pplt(LocalPIR R) : _R_(R){}
	bool operator() (typename LocalPIR::Element a, typename LocalPIR::Element b)
	{
		if ( b == 0 ) return true;
		else if ( a == 0 ) return false;
		else return a <= b;
	}
	//protected:
	LocalPIR _R_;
};

template <class LocalPIR>
static bool testLocalSmith (const LocalPIR &R, vector<typename LocalPIR::Element>& d, string s)
{
	typedef typename LocalPIR::Element Elt;
	typedef BlasMatrix<LocalPIR> Blackbox;

	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << s << endl;

	MatrixDomain<LocalPIR> MR (R);
	VectorDomain<LocalPIR> VD (R);

	bool ret = true;
	size_t i,j;
	size_t n = d.size();

	report << "Input vector:  ";
	VD.write (report, d);
	report << endl;

	// set up A equiv diag d.
	Blackbox L (R, n, n), D (R, n, n), U (R, n, n), A (R, n, n);
	for( i = 0; i < n; ++i )
		{ R.assign(D[i][i], d[i]); L[i][i]=U[i][i]=R.one; }
	for (i = 0; i < n; ++ i)
		for ( j = 0; j < i; ++ j) {
			R.assign(D[i][j], R.assign(D[j][i], 0));
			R.assign(L[i][j], rand() % 10);
			R.assign(L[j][i], 0);
			R.assign(U[j][i], rand() % 10);
			R.assign(U[i][j], 0);
		}
		
	MR.mul(A,L,D);
	MR.mulin(A,U);

	list< Elt > Inv;
	SmithFormLocal< LocalPIR > SmithForm;
	//timer.start();
	SmithForm( Inv, A, R );
	//timer.stop();
	//report << "Time " << timer <<"\n"; report.flush();

	report << "Computed invariants: ";
	report << "[";
	typedef typename list<Elt>::iterator listptr;
	for (listptr p = Inv.begin(); p != Inv.end(); ++p)
		report << *p << ", ";
	//report << "\b\b]" << endl;
	report << "normalize done" << endl; report.flush();

	// figure true invariants
	pplt<LocalPIR> lt(R);
	for (size_t i = 0; i < d.size(); ++i) normal(d[i], R);
	stable_sort(d.begin(), d.end(), lt);
	report << "True invariants: ";
	VD.write (report, d) << endl; report.flush();

	typename vector<Elt>::iterator q;
	listptr p;
	for (p = Inv.begin(), q = d.begin(); q != d.end(); ++p, ++q)
	{
		if ( !R.areEqual (*p, *q ) ) {
			report << "ERROR: Computed invariants incorrect" << endl;
			ret = false;
		}
		//commentator().progress();
	}
	return ret;
}

template<typename Base>
bool check_ranks(const std::vector<std::pair<size_t,Base> >& local,
                 const std::map<int, size_t>& map_values,
                 const Base& p) {
    bool pass(true);

    auto li(local.begin());
    auto mi(map_values.begin());
    for (; li != local.end(); ++li, ++mi) {
        if ( (li->first != mi->second) ||
             (li->second != Givaro::power(p,mi->first)) ) {
            std::cerr << "*** ERROR *** (" <<
                li->first << ',' << li->second << ')' 
                      << " ---->---- (" <<
                mi->second << ',' << mi->first << ')' 
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
    std::vector<std::pair<size_t,Base> > local;
    Permutation<ModRing> Q(B.field(),B.coldim());
    PGD(local, B, Q, Givaro::power(p,exp), p, PRESERVE_UPPER_MATRIX);

#ifdef LINBOX_LOCAL_SMITH_OUTPUT_
    for (auto ip = local.begin(); ip != local.end(); ++ip)
        std::cout << '[' << ip->first << ',' << ip->second << "] ";
    cout << endl;
#endif

    ModRing F(p);
    SparseMatrix<ModRing> BF(F,M,N);
    MatrixHom::map(BF,B);
    size_t rr;
    LinBox::rank(rr,BF, Method::SparseElimination() );

#ifdef LINBOX_LOCAL_SMITH_OUTPUT_
    cout << "Residue rank: " << rr << endl;
#endif

    return check_ranks(local,map_values,p) && (rr == R);
}


template<typename Base, typename SparseMat>
bool sparse_local_smith_poweroftwo(SparseMat& B,
                        size_t R, size_t M, size_t N,
                        const Base& p, int exp,
                        const std::map<int, size_t>& map_values) {
    LinBox::PowerGaussDomainPowerOfTwo< Base > PGD;
    LinBox::GF2 F2;
    Permutation<GF2> Q(F2,B.coldim());
    std::vector<std::pair<size_t,Base> > local;
    PGD(local, B, Q, exp, PRESERVE_UPPER_MATRIX);

#ifdef LINBOX_LOCAL_SMITH_OUTPUT_
    for (auto ip = local.begin(); ip != local.end(); ++ip)
        std::cout << '[' << ip->first << ',' << ip->second << "] ";
    cout << ")" << endl;
#endif

    return check_ranks(local,map_values,p);
}

template<typename Base, typename Compute = Base>
bool test_sparse_local_smith(size_t seed, size_t R, size_t M, size_t N,
                             const Base& p, int exp,
                             double density) {
	commentator().start ("Testing sparse elimination local smith", "SELS");
    commentator().indent(std::cerr);
    std::cerr << "test_sparse_local_smith(" << density << ", r" << R << '[' << M << 'x' << N << "]) mod " << p << '^' << exp << std::endl;

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

#ifdef LINBOX_LOCAL_SMITH_OUTPUT_
    for(auto & it : map_values)
        std::cout << it.first << ": " << it.second << " times\n";
#endif

    std::vector<int> indices(R);
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937 gen; gen.seed(seed);
    std::shuffle(indices.begin(), indices.end(), gen);

    std::uniform_real_distribution<> nonzero(0., 1.);

    typename ModRing::RandIter G(F,0,seed);
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
#ifdef DEBUG
        L.write(cerr<<"L:=",Tag::FileFormat::Maple) << ';' << endl;
#endif
    DenseMatrix<ModRing> U(F,N,N);
    for (size_t i=0 ; i<std::min(M,N); ++i){
        U.setEntry(i,i,F.one);
        for (size_t j= i+1; j<N ;++j)
            if (nonzero(gen)<density) G.random ( U.refEntry(i,j) );
    }
#ifdef DEBUG
        U.write(cerr<<"U:=",Tag::FileFormat::Maple) << ';' << endl;
#endif

    FFLAS::fgemm (F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans, M,N,N, F.one, L.getPointer(), N, U.getPointer(), N, F.zero, A.getPointer(), N);

#ifdef DEBUG
        A.write(cerr<<"A:=",Tag::FileFormat::Maple) << ';' << endl;
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
    bool pass0(true), pass1(true), pass2(true);
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
    FFLAS::writeCommandString(std::cout, args) << std::endl;

    commentator().setMaxDetailLevel (-1);
    commentator().setMaxDepth (-1);
    commentator().setReportStream (std::cerr);
	commentator().start("Local Smith Form test suite", "LocalSmith");
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

    { // sparseelim
        pass0 &= test_sparse_local_smith(rseed,r,m,n,Givaro::Integer(2),e,d);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,Givaro::Integer(q),e,d/2.);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,uint64_t(2),e,d);
        pass0 &= test_sparse_local_smith(rseed,r,m,n,uint64_t(q),e,d/2.);
        pass0 &= ruint_test<6>(rseed,r,m,n,q,e,d/2.);
        pass0 &= ruint_test<6>(rseed,r,m,n,2,e,d);
        pass0 &= ruint_test<7>(rseed,r,m,n,q,e,d/2.);
        pass0 &= ruint_test<7>(rseed,r,m,n,2,e,d);
        pass0 &= ruint_test<8>(rseed,r,m,n,q,e,d/2.);
        pass0 &= ruint_test<8>(rseed,r,m,n,2,e,d);

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


#if 1
  { // zero-th local ring type: modulus p^e as integer.
	typedef LocalPIRModular<integer> Ring;
	Ring R (q, e);
	vector<Ring::Element> d(n);

	commentator().start ("Testing local smith over LocalPIRModular<integer>", "local smith");
	for( int32_t i = 0; i < n; ++i ) R.init(d[i],i);
	pass1 = testLocalSmith<Ring> (R, d, "LocalPIRModular<integer>");
	commentator().stop ("local smith");
	if (not pass1) report << "LocalPIRModular smith FAIL" << std::endl;

  }
#endif
  { // first local ring type: modulus p^e as int32_t.
	typedef PIRModular<int32_t> LocalPIR;
	//typedef PIRModular<dense> LocalPIR;
	//LocalPIR R (81); 
	LocalPIR R (q, e);
	vector<LocalPIR::Element> d(n);

	commentator().start ("Testing local smith on singular dense mat over PIRModular", "testSingular");
	for( int32_t i = 0; i < n; ++i ) R.init(d[i],i);
	pass1 = testLocalSmith<LocalPIR> (R, d, "PIRModular<int32_t>");
	commentator().stop ("testSingular");
	if (not pass1) report << "PIRModular sing FAIL" << std::endl;

	commentator().start ("Testing local smith on nonsingular dense mat over PIRModular", "testNonsingular");
	LocalPIR::RandIter r(R);
	LocalPIR::Element e; R.init(e);
	for( int32_t i = 0; i < n; ++i ) {	
		r.random(e);
		do { r.random(e); } while (R.isZero(e));
		R.assign(d[i], e);
	}
	bool p = testLocalSmith<LocalPIR> (R, d, "PIRModular<int32_t>");
	if (not p) report << "PIRModular nonsing FAIL" << std::endl;
	commentator().stop ("testNonsingular");
	pass1 = pass1 and p;
	if (not pass1) report << "PIRModular FAIL" << std::endl;
  }

  { // second local ring type: m = 2^32
	typedef Local2_32 LocalPIR;
	LocalPIR R;
	vector<LocalPIR::Element> d(n);

	commentator().start ("Testing local smith on singular dense mat over Local2_32", "testSingular");
	for( int64_t i = 0; i < n; ++i )
		d[i] = (LocalPIR::Element) i;
	if (!testLocalSmith<LocalPIR> (R, d, "Local2_32")) pass2 = false;
	commentator().stop ("testSingular");

	commentator().start ("Testing local smith on nonsingular dense mat over Local2_32", "testNonsingular");
	for( int64_t i = 0; i < n; ++i )
		d[i] = (LocalPIR::Element) i+1;
	if (!testLocalSmith<LocalPIR> (R, d, "Local2_32")) pass2 = false;
	commentator().stop ("testNonsingular");
	if (not pass2) report << "Local2_32 FAIL" << std::endl;
  }

	commentator().stop("Local Smith Form test suite");
	return pass0 and pass1 and pass2 ? 0 : -1;
}


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

