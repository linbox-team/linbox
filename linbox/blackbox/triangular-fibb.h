/* linbox/blackbox/triangular-fibb.h
 * bds
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

#ifndef __LINBOX_triangular_fibb_H
#define __LINBOX_triangular_fibb_H

#include <utility>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/linbox-tags.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/blackbox/fibb.h"

namespace LinBox
{

	/** \brief

	  \ingroup blackbox
	 */
	template<class Field_>
	struct TriangularFIBB : public  FIBB<Field_>
	{
		typedef Field_			Field;
	protected:
		typedef TriangularBlasMatrix<BlasMatrix<Field>> Rep_t;
		Rep_t * rep_;
		const Field* field_;
	public:
		typedef TriangularFIBB<Field>	Self_t;
		typedef typename Field::Element	Element;
		typedef DenseMatrix<Field> Matrix;

		Self_t& init(Rep_t& T) const
		{ rep_ = &T; field_ = &(T.field()); return *this; }

		TriangularFIBB(const Field& F): field_(&F) {}
		/// Constructor from a TriangularBlasMatrix
		TriangularFIBB(Rep_t& T): rep_(&T), field_(&T.field()) {}

		void random() { rep_->random(); }
		/*
		void random(size_t n)
		{
			size_t n = rowdim();
			rep_->identity((int)n);
			MersenneTwister r((unsigned int)std::time(nullptr));
			// Knuth construction
			for (size_t i = 0; i < n-1; ++i) {
				size_t j = i + r.randomInt()%(n-i);
				std::swap(rep_->_indices[i], rep_->_indices[j]);
			}
		}
		*/


		/* Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		TriangularFIBB (const TriangularFIBB &Mat) :
			rep_(Mat.rep_)
		{}

		// Destructor
		~TriangularFIBB (void) {}

		Matrix& applyRight(Matrix& B, const Matrix& A) const
		{	return applyMatrix(B, A, FFLAS::FflasLeft);   }

		Matrix& applyLeft(Matrix& B, const Matrix& A) const
		{	return applyMatrix(B, A, FFLAS::FflasRight);   }

		Matrix& applyMatrix(Matrix& B, const Matrix& A, const FFLAS::FFLAS_SIDE side) const {
			B=A;
			FFLAS::ftrmm<Field> (field(), 
				side,
				FFLAS::FFLAS_UPLO(rep_->getUpLo()), //const FFLAS_UPLO Uplo,
		    	FFLAS::FflasNoTrans, //const FFLAS_TRANSPOSE TransA,
				FFLAS::FFLAS_DIAG(rep_->getDiag()), //const FFLAS_DIAG Diag,
				B.rowdim(), //const size_t M,
				B.coldim(), //const size_t N,
				field().one, //const typename Field::Element alpha,
				rep_->getPointer(), //typename Field::Element_ptr A,
				rep_->getStride(), //const size_t lda,
				B.getPointer(), //typename Field::Element_ptr B,
				B.getStride() //const size_t ldb)
			);
			return B;
		}

#if 0
		// will want better base case
		Matrix& applyRight(Matrix& Y, const Matrix& X) const
		{
			Element x; field().init(x);
			if (rowdim() == 1) ...//TODO == quit here
			for (size_t i = 0; i < Y.rowdim(); ++i){
				size_t k = rep_->_indices[i];
				for (size_t j = 0; j < Y.coldim(); ++j)
					Y.setEntry(i,j, X.getEntry(x, k, j));
			}
		/* desired form
			for (size_t i = 0; i < rowdim(); ++i)
			{
				Matrix Yrow(Y, i, 0, 1, Y.coldim());
				Matrix Xrow(X, rep_->_indices[i], 0, 1, X.coldim());
				Yrow.copy(Xrow); // right kind of copy?
			}
		*/
			return Y;
		}
		Matrix& applyLeft(Matrix& Y, const Matrix& X) const
		{
			Element x; field().init(x);
			for (size_t i = 0; i < Y.coldim(); ++i){
				size_t k = rep_->_indices[i];
				for (size_t j = 0; j < Y.rowdim(); ++j)
					Y.setEntry(j,k, X.getEntry(x, j, i));
			}
		/* desired form
			for (size_t i = 0; i < coldim(); ++i)
			{
				Matrix Ycol(Y, 0, rep_->_indices[i], Y.rowdim(), 1);
				Matrix Xcol(X, 0, i, X.rowdim(), 1);
				Ycol.copy(Xcol);
		*/
			return Y;
		}
#endif
		/* FIBB functions */

		BBType bbTag() const { return triangular; }

		size_t& rank(size_t& r) const
		{	size_t m = rowdim(), n = coldim();
			size_t k = m < n ? m : n;
			if (rep_->getDiag() == Tag::Diag::Unit) {
				r = k;
			} else {
				// assume trapezoid for now, later fix to allow echelon
				Element x;
				r = 0;
				while (r < k and not field().isZero(rep_->getEntry(x,r,r)))
					++r;
			}
			return r;
		}

		Element& det(Element& d) const
		{	size_t m = rowdim(), n = coldim();
			if (m != n) {
				field().assign(d,field().zero);
			} else if (rep_->getDiag() == Tag::Diag::Unit) {
				field().assign(d,field().one);
			} else {
				Element x;
				field().assign(d, field().one);
				for (size_t i = 0; i < m; ++i)
					field().mulin(d,rep_->getEntry(x,i,i));
			}
			return d;
		}

		Matrix& solveRight(Matrix& B, const Matrix& A) const
		{	return solveMatrix(B, A, Tag::Side::Left);	}

		Matrix& solveLeft(Matrix& B, const Matrix& A) const
		{	return solveMatrix(B, A, Tag::Side::Right);	}

		Matrix& solveMatrix(Matrix& B, const Matrix& A, const Tag::Side side) const
		{
			//B.copy(A);
            B=A;
			FFLAS::ftrsm<Field> (field(), FFLAS::FFLAS_SIDE(side),
				FFLAS::FFLAS_UPLO(rep_->getUpLo()), //const FFLAS_UPLO Uplo,
		    	FFLAS::FflasNoTrans, //const FFLAS_TRANSPOSE TransA,
				FFLAS::FFLAS_DIAG(rep_->getDiag()), //const FFLAS_DIAG Diag,
				B.rowdim(), //const size_t M,
				B.coldim(), //const size_t N,
				field().one, //const typename Field::Element alpha,
				rep_->getPointer(), //typename Field::Element_ptr A,
				rep_->getStride(), //const size_t lda,
				B.getPointer(), //typename Field::Element_ptr B,
				B.getStride() //const size_t ldb)
			);
			return B;
		}

#if 0
		Matrix& solveRight(Matrix& Y, const Matrix& X) const
		{	Element x; field().init(x);
			for (size_t i = 0; i < Y.rowdim(); ++i){
				size_t k = rep_->_indices[i];
				for (size_t j = 0; j < Y.coldim(); ++j)
					Y.setEntry(k,j, X.getEntry(x, i, j));
			}
		/* desired form
		 	for (size_t i = 0; i < rowdim(); ++i)
			{
				Matrix Yrow(Y, rep_->_indices[i], 0, 1, Y.coldim());
				Matrix Xrow(X, i, 0, 1, X.coldim());
				Yrow.copy(Xrow);
			}
		*/
			return Y;
		}
		Matrix& solveLeft(Matrix& Y, const Matrix& X) const
		{	Element x; field().init(x);
			for (size_t i = 0; i < Y.coldim(); ++i){
				size_t k = rep_->_indices[i];
				for (size_t j = 0; j < Y.rowdim(); ++j)
					Y.setEntry(j,i, X.getEntry(x, j, k));
			}
		/* desired form
			for (size_t i = 0; i < coldim(); ++i)
			{
				Matrix Ycol(Y, 0, i, Y.rowdim(), 1);
				Matrix Xcol(X, 0, rep_->_indices[i], X.rowdim(), 1);
				Ycol.copy(Xcol);
			}
		*/
			return Y;
		}
#endif
		Matrix& nullspaceRandomRight(Matrix& N) const
		{	N.zero(); return N; }
		Matrix& nullspaceRandomLeft(Matrix& N) const
		{	N.zero(); return N; }
		BlasMatrix<Field>& nullspaceBasisRight(BlasMatrix<Field>& N) const
		{	N.resize(rowdim(), 0); return N; }
		BlasMatrix<Field>& nullspaceBasisLeft(BlasMatrix<Field>& N) const
		{	N.resize(0, coldim()); return N; }
		/* end FIBB section */

		template<typename _Tp1>
		struct rebind {
			typedef TriangularFIBB<_Tp1> other;
			void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
				Ap->setStorage( A.getStorage() );
			}
		};



		/* Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		/// rowdim
		size_t rowdim (void) const
		{
			return rep_->rowdim();
		}

		/* Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		/// coldim
		size_t coldim (void) const
		{
			return rep_->coldim();
		}

		const Field& field() const { return rep_->field(); }
		const Field& ring() const { return rep_->field(); }

		//!@bug needs a MM version
		std::ostream &write(std::ostream &os) const
		{ return rep_->write(os, Tag::FileFormat::Plain); }

		std::ostream &write(std::ostream &os, Tag::FileFormat format) const
		{ return rep_->write(os, format); }

		//!@bug there is no read here. (needed by test-blackbox.h)
		std::istream &read(std::istream &is)
		{ return is; }
		std::istream &read(std::istream &is) const
		{ return is; }
		std::istream &read(std::istream &is, Tag::FileFormat format) const
		{
			throw NotImplementedYet();
			return is ;
		}

}; // template <Vector> class TriangularFIBB


} // namespace LinBox

#endif // __LINBOX_triangular_fibb_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
