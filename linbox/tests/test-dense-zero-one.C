/* Copyright (C) LinBox
 *
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */



/*! @file  tests/test-dense-zero-one.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */




#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/field/modular.h"

#include "linbox/blackbox/dense-zero-one.h"


#include "linbox/util/timer.h"
#include "test-common.h"

namespace LinBox
{
template <class _Field>
struct BlackboxDomain : public _Field
{
	// types Scalar, Block, Blackbox
	typedef _Field Field; // transitional
	typedef typename _Field::Element Scalar;
	typedef typename _Field::Element Element; // transitional
	typedef typename _Field::RandIter RandIter;
	typedef BlasMatrix<Field> Block;
	/*struct Block: public BlasMatrix<Field> {
		Block( int m, int n ): BlasMatrix<Field> ( m, n ) {}
		Block& subBlock( Block & B, size_t i, size_t j, size_t m, size_t n ) {
			return B = BlasMatrix<Field>(*this, i, j, m, n);
		}
		const Block& subBlock( const Block & B, size_t i, size_t j, size_t m, size_t n ) {
			return B = BlasMatrix<Field>(static_cast<BlasMatrix<Field> >(*this), i, j, m, n);
		}
		Block & operator= (BlasMatrix<Field> & rhs){
			*this = rhs;
			return *this;
		}
	}; */
	BlackboxDomain(){}
	BlackboxDomain(long q): Field(q) {}

	// A <-- A + B.  A and B must have same shape.
	Block & addin( Block & A, const Block & B ) const { return MatrixDomain<Field>(*this).addin(A,B); }
	// for the Scalar addin.
	using Field::addin;

	// A <-- B * C.  A, B, C must have compatible shapes.
	Block & mul( Block & C, const Block & A, const Block & B ) const {
		return BlasMatrixDomain<Field>(*this).mul(C,A,B);
	}

	// for the Scalar mul.
	using Field::mul;

	// C <-- B * C.  A, B, C must have compatible shapes.
	Block & axpyin( Block & C, const Block & A, const Block & B ) const {
		return BlasMatrixDomain<Field>(*this).axpyin(C,A,B);
	}

	// for the Scalar axpyin
	using Field::axpyin;

	// A <-- B * C.  A, B, C must have compatible shapes.
	bool areEqual( const Block & A, const Block & B ) const {
		return MatrixDomain<Field>(*this).areEqual(A,B);
	}
	// for the Scalar areEqual.
	using Field::areEqual;

	// Set the entries in a block to zero.
	void zero( Block& B ) const {
		for ( typename Block::Iterator raw = B.Begin(); raw != B.End(); ++raw )
			init(*raw, 0);
	}

	// Set the entries in a block to random field elements.
	void random( Block& B) {
		RandIter r(*this);
		for ( typename Block::RowIterator row = B.rowBegin(); row != B.rowEnd(); ++row )
			for ( typename Block::Row::iterator place = row->begin(); place != row->end(); ++place )
				r.random(*place);
	}

}; //BlackboxDomain

} // LinBox

using namespace std;
using namespace LinBox;

template <class Blackbox>
bool testAssociativity(Blackbox& A)
{
	typedef typename Blackbox::MatrixDomain Dom;
	Dom MD = A.domain();
	size_t m = A.rowdim(), n = A.coldim() - 100;
	size_t k = (m + n)/2;
	typename Dom::Block B(A.field(),k,m), C(A.field(),m,n);
	MD.random(B); MD.random(C);

	typename Dom::Block D(A.field(),m,n), E(A.field(),k,n);

	A.apply(D, C); // D = AC
	MD.mul(E,B,D); // E = B(AC)

	typename Dom::Block F(A.field(),k,m), G(A.field(),k,n);

	A.unpackingApplyTranspose(F,B); // F = BA
	MD.mul(G,F,C); // G = (BA)C
	return MD.areEqual(E,G);


} // testAssociativity

template <class Blackbox>
void testTiming(Blackbox & A)
{
	typedef typename Blackbox::MatrixDomain Dom;
	typedef typename Dom::Block Block;

	Dom MD = A.domain();
	size_t m = A.rowdim(), n = A.coldim();
	size_t k = (m + n)/2;

	LinBox::UserTimer timer;

	Block B(n,k), C(m,k), D(k,m), E(k,n), F(k,k);
	MD.random(B); MD.random(D);

	vector<typename Dom::Element> v1, v2(m);
	typename Dom::RandIter r(MD);
	typename Dom::Element x;
	for(size_t i = 0; i != n; ++i){
		r.random(x);
		v1.push_back(x);
	}


	//Tests:
	cout << "Timing tests:" << endl << endl;

	timer.clear(); timer.start();
	for(size_t j = 0; j != m; ++j) A.apply(v2,v1);
	timer.stop();
	cout << "apply using vectors time: " << timer << endl;

	timer.clear(); timer.start();
	A.applyTranspose(C,B);
	timer.stop();
	cout << "apply using row addin time: " << timer << endl;

	timer.clear(); timer.start();
	A.unpackingApplyTranspose(C,B);
	timer.stop();
	cout << "apply using block axpy time: " << timer << endl;

	timer.clear(); timer.start();
	MD.mul(F, D, C);
	timer.stop();
	cout << "Matrix Domain mul time: " << timer << endl;

	cout << "End of timing tests" << endl << endl;

} // testTiming


template <class Blackbox>
void blockSizeTimingTest(Blackbox & A, size_t size)
{
	typedef typename Blackbox::MatrixDomain Dom;
	typedef typename Dom::Block Block;

	Dom MD = A.domain();
	size_t m = A.rowdim();

	LinBox::UserTimer timer;

	Block B(m,m), C(m,m), D(m,m);
	MD.random(B); MD.random(D);

	cout << size << "         " << m << "       ";

	timer.clear(); timer.start();
	A.unpackingApply(C,B,size);
	timer.stop();
	cout << timer << "              ";

	timer.clear(); timer.start();
	A.unpackingApplyTranspose(C,B,size);
	timer.stop();
	cout << timer << "               ";

	timer.clear(); timer.start();
	MD.mul(C,D,B);
	timer.stop();
	cout << timer << "     ";

	cout << endl;

} //blockSizeTimingTest()


template <class Blackbox>
void stressTest (Blackbox & A)
{
	//Note that the rowdim/coldim of A must be 30000
	typedef typename Blackbox::MatrixDomain Dom;
	typedef typename Dom::Block Block;

	Dom MD = A.domain();
	size_t m = 30000;
	size_t n = 2000;

	LinBox::UserTimer timer;

	Block B(m,n), C(m,n);
	MD.random(B);

	cout << "Test: 30000x30000 matrix multiplied by 30000x2000 matrix\nblock size = 2048\n\n";

	timer.clear(); timer.start();
	A.unpackingApply(C,B,2048);
	timer.stop();
	cout << "unpacking apply time: " << timer << endl;

	Block D(m,m); MD.random(D);

	timer.clear(); timer.start();
	MD.mul(C,D,B);
	timer.stop();
	cout << "domain mul time: " << timer << endl;

	cout << endl;
}  //end stressTest()

template <class Blackbox>
void largeTest (Blackbox & A)
{
	//Use for large blackboxes
	typedef typename Blackbox::MatrixDomain Dom;
	typedef typename Dom::Block Block;

	Dom MD = A.domain();
	size_t m = A.coldim();
	size_t n = 2000;

	LinBox::UserTimer timer;

	Block B(m,n), C(m,n);
	MD.random(B);

	cout << "Test: " << A.rowdim() << "x" << m << "blackbox multiplied by " << m << "x" << n << "block\nblock size: 2048\n\n";

	timer.clear(); timer.start();
	A.unpackingApply(C,B,2048);
	timer.stop();
	cout << "unpacking apply time: " << timer << endl;
}  //end largeTest


int main (int argc, char* argv[])
{

	static size_t n = 100;
        static integer f = 4093;
        static integer d = 1000003;
	static bool t = false;

        static Argument args[] = {
                { 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
                { 'f', "-f Q", "Operate over the \"field\" GF(Q).", TYPE_INTEGER, &f },
                { 'd', "-d Q", "Operate over the \"field\" GF(Q).", TYPE_INTEGER, &d },
                { 't', "-t T", "If T is flagged, only run timing tests.", TYPE_BOOL, &t },
		END_OF_ARGUMENTS
        };

        parseArguments (argc, argv, args);

	// create a dense-zero-one
	typedef Modular<float> FieldF;
	typedef Modular<double> FieldD;

	BlackboxDomain<FieldF> F(f);
	FieldF::RandIter r1(F);

	BlackboxDomain<FieldD> D(d);
	FieldD::RandIter r2(D);

	DenseZeroOne<BlackboxDomain<FieldF> > A(F, n, n);
	DenseZeroOne<BlackboxDomain<FieldD> > B(D, n, n);

	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j) {
			if (rand()%2){
				A.setEntry(i, j, F.one);
				B.setEntry(i, j, D.one);
			}

			else {
				A.setEntry(i, j, F.zero);
				B.setEntry(i, j, D.zero);
			}
		}

	//A.write(cout) << endl << endl;
	//B.write(cout) << endl << endl;

#if 0
	/* basic everyday test */
	cout << endl;
	cout << "Domain: Modular<float>, GF(" << f << ")" << endl;
	testTiming(A);

	cout << "Domain: Modular<double>, GF(" << d << ")" << endl;
	testTiming(B);

	/* block size tests */
	cout << endl;
	cout << "Domain: Modular<double>, GF(" << d << ")" << endl << endl;
	cout << "block Size     n     unpackingApply   unpackingApplyTranspose    Domain mul" << endl << endl;
	for (size_t sizeU = 256; sizeU != 4096; sizeU *= 2)
		for (int count = 500; count != 4500; count += 500){
			DenseZeroOne<BlackboxDomain<FieldD> > C(D, count, count);

			for (size_t i = 0; i < count; ++i)
				for (size_t j = 0; j < count; ++j) {
					if (rand()%2)
						C.setEntry(i, j, D.one);
					else
						C.setEntry(i, j, D.zero);
				}

				blockSizeTimingTest(C, sizeU);
		}

#endif

#if 0
	/* stress test */
	cout << "Domain: Modular<double>, GF(" << d << ")" << endl;
	DenseZeroOne<BlackboxDomain<FieldD> > C(D, 30000, 30000);
	for (size_t i = 0; i != 30000; ++i)
		for (size_t j = 0; j != 30000; ++j) {
			if (rand()%2)
				C.setEntry(i, j, D.one);
			else
				C.setEntry(i, j, D.zero);
		}
	stressTest(C);

	/* large tests */
	cout << "Domain: Modular<double>, GF(" << d << ")" << endl;
	for (size_t dim = 120000; dim < 150000; dim *= 2){
		DenseZeroOne<BlackboxDomain<FieldD> > C(D, dim, dim);
		for (size_t i = 0; i != dim; ++i)
			for (size_t j = 0; j != dim; ++j) {
				if (rand()%2)
					C.setEntry(i, j, D.one);
				else
					C.setEntry(i, j, D.zero);
			}
		cout << "blackbox created\n";
		largeTest(C);
	}
#endif

	std::cout << "most test is commented out" << std::endl;
	if (!t){
		bool pass = testAssociativity(B);

		return pass ? 0 : -1;
	}

}


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

