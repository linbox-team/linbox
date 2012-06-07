/* linbox/matrix/factorized-matrix.inl
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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

#ifndef __LINBOX_factorized_matrix_INL
#define __LINBOX_factorized_matrix_INL

namespace LinBox
{
	namespace Protected
	{
		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixLeftSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixLeftSolve

		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixRightSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixRightSolve

		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixLeftLSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixLeftLSolve

		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixRightLSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixRightLsolve

		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixLeftUSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixLeftUSolve

		/// @internal
		template <class Field, class Operand>
		class FactorizedMatrixRightUSolve {
		public:
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& X, const Operand& B ) const;
			Operand& operator() ( const Field& F,
					      const LQUPMatrix<Field>& A,
					      Operand& B ) const;
		}; // end of class FactorizedMatrixRightUSolve

		/*
		 * Solvers with Matrices: Operand=BlasMatrix<Field>
		 */

		template <class Field>
		class FactorizedMatrixLeftSolve<Field, BlasMatrix<Field> > {
		public:

			BlasMatrix<Field>& operator() (const Field& F,
						       const LQUPMatrix<Field>& A,
						       BlasMatrix<Field>& X,
						       const BlasMatrix<Field>& B) const
			{
				linbox_check (A.coldim() == X.rowdim());
				linbox_check (A.rowdim() == B.rowdim());
				linbox_check (B.coldim() == X.coldim());
				int info;

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasLeft, A.rowdim(), A.coldim(), B.coldim(), A.getRank(),
						A.getPointer(), A.getStride(), A.getP().getPointer(), A.getQ().getPointer(),
						X.getPointer(), X.getStride(),
						B.getPointer(), B.getStride(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return X;
			}

			BlasMatrix<Field>& operator() (const Field& F,
						       const LQUPMatrix<Field>& A,
						       BlasMatrix<Field>& B) const
			{

				int info;
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.coldim() == B.rowdim());

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasLeft, B.rowdim(), B.coldim(), A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						B.getPointer(), B.getStride(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return B;
			}
		}; // end of class FactorizedMatrixLeftSolve

		template <class Field>
		class FactorizedMatrixRightSolve<Field, BlasMatrix<Field> > {
		public:

			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& X,
							const BlasMatrix<Field>& B ) const
			{
				linbox_check (A.rowdim() == X.coldim());
				linbox_check (A.coldim() == B.coldim());
				linbox_check (B.rowdim() == X.rowdim());
				int info;

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasRight, A.rowdim(), A.coldim(), B.rowdim(), A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						X.getPointer(), X.getStride(),
						B.getPointer(), B.getStride(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return X;
			}

			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& B ) const
			{

				int info;
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.rowdim() == B.coldim());

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasRight, B.rowdim(), B.coldim(), A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						B.getPointer(), B.getStride(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return B;
			}

		}; // end of class FactorizedMatrixRightSolve

		template <class Field>
		class FactorizedMatrixLeftLSolve<Field, BlasMatrix<Field> > {
		public:
			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& X,
							const BlasMatrix<Field>& B ) const
			{
				linbox_check (A.rowdim() == B.rowdim());
				X = B;
				return  (*this)(F, A, X);
			}

			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& B ) const
			{

				linbox_check (A.rowdim() == B.rowdim());

				FFPACK::solveLB2 ((typename Field::Father_t)F, FFLAS::FflasLeft, B.rowdim(), B.coldim(), A.getRank(),
						  A.getPointer(), A.getStride(),
						  A.getQ().getPointer(),
						  B.getPointer(), B.getStride());

				return B;
			}
		}; // end of class FactorizedMatrixLeftLSolve

		template <class Field>
		class FactorizedMatrixRightLSolve<Field, BlasMatrix<Field> > {
		public:
			BlasMatrix<Field>& operator() (const Field& F,
						       const LQUPMatrix<Field>& A,
						       BlasMatrix<Field>& X,
						       const BlasMatrix<Field>& B) const
			{
				linbox_check (A.rowdim() == B.coldim());
				X = B;
				return  (*this)( F, A, X );
			}

			BlasMatrix<Field>& operator() (const Field& F,
						       const BlasMatrix<Field>& A,
						       BlasMatrix<Field>& B) const
			{

				linbox_check( A.rowdim() == B.coldim() );

				FFPACK::solveLB2 ((typename Field::Father_t)F, FFLAS::FflasRight, B.rowdim(), B.coldim(), A.getRank(),
						  A.getPointer(), A.getStride(),
						  A.getQ().getPointer(), B.getPointer(), B.getStride());
				return B;
			}
		}; // end of class FactorizedMatrixRightLsolve

		template <class Field>
		class FactorizedMatrixLeftUSolve<Field, BlasMatrix<Field> > {
		public:

			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& X,
							const BlasMatrix<Field>& B ) const
			{

				linbox_check (A.coldim() == X.rowdim());
				linbox_check (B.coldim() == X.coldim());
				linbox_check (A.rowdim() == B.rowdim());

				bool consistent = true;
				size_t ldb = B.getStride();
				size_t ldx = X.getStride();
				typename Field::Element * Bp = B.getPointer();
				typename Field::Element * Xp = X.getPointer();

				for (size_t i = A.getRank(); i < B.rowdim(); ++i)
					for (size_t j = 0; j < B.coldim(); ++j)
						if (!F.isZero (*(Bp + i*ldb + j)))
							consistent = false;
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				// The last rows of B are now supposed to be 0

				for (size_t i=0; i < A.getRank(); ++i)
					FFLAS::fcopy ((typename Field::Father_t)F, B.coldim(), Xp + i*ldx, 1, Bp + i*ldx,1);

				FFLAS::ftrsm ((typename Field::Father_t)F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      A.getRank(), X.coldim(), F.one, A.getPointer(), A.getStride(), X.getPointer(), X.getStride());

				for (size_t i=A.getRank(); i < X.rowdim(); ++i)
					for (size_t j = 0; j < X.coldim(); ++j)
						F.assign (*(Xp + i*ldx + j), F.zero);

				return X;

			}
			BlasMatrix<Field>& operator() ( const Field& F,
							const BlasMatrix<Field>& A,
							BlasMatrix<Field>& B ) const
			{

				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.coldim() == B.rowdim());
				typename Field::Element * Bp = B.getPointer();
				size_t ldb = B.getStride();
				bool consistent = true;

				for (size_t i = A.getRank(); i < B.rowdim(); ++i)
					for (size_t j = 0; j < B.coldim(); ++j)
						if (!F.isZero (*(Bp + i*ldb + j)))
							consistent = false;
				if (!consistent)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				FFLAS::ftrsm ((typename Field::Father_t)F, FFLAS::FflasLeft, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      A.getRank(), B.coldim(), F.one, A.getPointer(), A.getStride(), Bp, ldb);

				return B;
			}

		}; // end of class FactorizedMatrixLeftUSolve

		template <class Field>
		class FactorizedMatrixRightUSolve<Field, BlasMatrix<Field> > {
		public:
			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& X,
							const BlasMatrix<Field>& B ) const
			{
				linbox_check (X.coldim() == A.rowdim());
				linbox_check (X.rowdim() == B.rowdim());
				linbox_check (A.coldim() == B.coldim());
				typename Field::Element * Bp = B.getPointer();
				typename Field::Element * Xp = X.getPointer();
				size_t R = A.getRank();
				size_t ldb = B.getStride();
				size_t ldx = X.getStride();

				for (size_t i = 0; i < X.getrowdim(); ++i)
					FFLAS::fcopy ((typename Field::Father_t)F, R, Xp + i*ldx, 1, Bp + i*ldb,1);

				FFLAS::ftrsm ((typename Field::Father_t)F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      X.rowdim(), R, F.one, A.getPointer(), A.getStride(), X.getPointer(), X.getStride());

				bool consistent = true;
				if (B.coldim() > X.coldim()) {
					typename Field::Element* W = new typename Field::Element [B.rowdim() * (B.coldim() - R)];
					size_t ldw = B.rowdim();

					FFLAS::fgemm ((typename Field::Father_t) F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						      X.rowdim(), B.coldim() - R, R,
						      F.one, Xp, X.getStride(), A.getPointer() + R, A.getStride,
						      F.zero, W, ldw);

					for (size_t i = 0; i < B.rowdim(); ++i)
						for (size_t j = 0; j < B.coldim()-R; ++j)
							if (!F.areEqual (*(W + i*ldw + j), *(Bp + R + i*ldb +j)))
								consistent = false;
					delete[] W;
				}
				else {
					FFLAS::fgemm ((typename Field::Father_t)F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						      X.rowdim(), B.coldim() - R, R,
						      F.one, Xp, X.getStride(), A.getPointer() + R, A.getStride,
						      F.zero, Xp + R, X.getStride());

					for (size_t i = 0; i < B.rowdim(); ++i)
						for (size_t j = R; j < B.coldim(); ++j)
							if (!F.areEqual (*(X + i*ldx + j), *(Bp + i*ldb +j)))
								consistent = false;
				}
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				for (size_t i = 0; i < X.rowdim(); ++i)
					for (size_t j = R; j < X.coldim(); ++j)
						F.assign (*(Xp + i*ldx + j), F.zero);
				return X;
			}

			BlasMatrix<Field>& operator() ( const Field& F,
							const LQUPMatrix<Field>& A,
							BlasMatrix<Field>& B ) const
			{
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (B.coldim() == A.rowdim());

				typename Field::Element * Bp = B.getPointer();
				size_t ldb = B.getStride();
				size_t R = A.getRank();

				FFLAS::ftrsm ((typename Field::Father_t)F, FFLAS::FflasRight, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      B.rowdim(), R, F.one, A.getPointer(), A.getStride(), B.getPointer(), B.getStride());

				FFLAS::fgemm ((typename Field::Father_t)F, FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					      B.rowdim(), B.coldim() - R, R,
					      F.one, Bp, B.getStride(), A.getPointer() + R, A.getStride,
					      F.mOne, Bp + R, B.getStride());

				bool consistent = true;
				for (size_t i = 0; i < B.rowdim(); ++i)
					for (size_t j = R; j < B.coldim(); ++j)
						if (!F.isZero (*(Bp + i*ldb + j)))
							consistent = false;
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return B;

			}
		}; // end of class FactorizedMatrixRightUSolve


		/*
		 * Solvers with vectors: Operand=std::vector<Element>
		 */

		template <class Field>
		class FactorizedMatrixLeftSolve<Field, std::vector<typename Field::Element> > {
		public:
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.coldim() == x.size());
				linbox_check (A.rowdim() == b.size());
				int info;

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasLeft, A.rowdim(), A.coldim(), 1, A.getRank(),
						A.getPointer(), A.getStride(), A.getP().getPointer(), A.getQ().getPointer(),
						&x[0], 1, &b[0], 1, &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return x;
			}

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{

				int info;
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.coldim() == b.size());

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasLeft, b.size(), 1, A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						&b[0], 1, &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return b;
			}
		}; // end of class FactorizedMatrixLeftSolve

		template <class Field>
		class FactorizedMatrixRightSolve<Field, std::vector<typename Field::Element> > {
		public:

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.rowdim() == x.size());
				linbox_check (A.coldim() == b.size());
				int info;

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasRight, A.rowdim(), A.coldim(), 1, A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						&x[0], x.size(), &b[0], b.size(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return x;
			}

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{

				int info;
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.rowdim() == b.size());

				FFPACK::fgetrs ((typename Field::Father_t)F, FFLAS::FflasRight, 1, b.size(), A.getRank(),
						A.getPointer(), A.getStride(),
						A.getP().getPointer(), A.getQ().getPointer(),
						&b[0], b.size(), &info);
				if (info > 0)
					throw LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return b;
			}
		}; // end of class FactorizedMatrixRightSolve

		template <class Field>
		class FactorizedMatrixLeftLSolve<Field, std::vector<typename Field::Element> > {
		public:
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check( A.rowdim() == b.size() );
				x = b;
				return (*this)( F, A, x );
			}

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.rowdim() == b.size());

				FFPACK::solveLB2 ((typename Field::Father_t)F, FFLAS::FflasLeft, b.size(), 1, A.getRank(),
						  A.getPointer(), A.getStride(),
						  A.getQ().getPointer(), &b[0], 1);

				return b;

#if 0 /* BB: unreachable  ! */
				   size_t n = b.size(); // bds: b not B
				   linbox_check( A.rowdim() == n );
				   size_t r = A.getRank();

				// To be changed: solveLB is designed for matrices, not for vectors
				if ( A.rowdim() <= A.coldim() ) {
				FFPACK::solveLB((typename Field::Father_t) F, FFLAS::FflasLeft, n, 1, r, A.getPointer(), A.getStride(),
				A.getQ().getPointer(), &b[0], b.size() );
				}
				else
				FFPACK::solveLB2((typename Field::Father_t) F, FFLAS::FflasLeft, n, 1, r, A.getPointer(), A.getStride(),
				A.getQ().getPointer(), b.getPointer(), b.getStride() );
				return b;

#endif
			}
		}; // end of class FactorizedMatrixLeftLSolve

		template <class Field>
		class FactorizedMatrixRightLSolve<Field, std::vector<typename Field::Element> > {
		public:
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check( A.rowdim() == b.size() );
				x = b;
				return (*this)( F, A, x );
			}

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.rowdim() == b.size());

				FFPACK::solveLB2 ((typename Field::Father_t)F, FFLAS::FflasRight, 1, b.size(),  A.getRank(),
						  A.getPointer(), A.getStride(),
						  A.getQ().getPointer(), &b[0], b.size());
				return b;
			}
		}; // end of class FactorizedMatrixRightLsolve

		template <class Field>
		class FactorizedMatrixLeftUSolve<Field, std::vector<typename Field::Element> > {
		public:
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.coldim() == x.size());
				linbox_check (A.rowdim() == b.size());

				bool consistent = true;
				typename Field::Element * bp = &b[0];           ;
				typename Field::Element * xp = &x[0];

				for (size_t i = A.getRank(); i < b.size(); ++i)
					if (!F.isZero (b[i]))
						consistent = false;
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				// The last rows of B are now supposed to be 0

				for (size_t i=0; i < A.getRank(); ++i)
					F.assign (x[i], b[i]);

				FFLAS::ftrsv ((typename Field::Father_t)F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      A.getRank(), A.getPointer(), A.getStride(), xp, 1);

				for (size_t i=A.getRank(); i < x.size(); ++i)
					F.assign (x[i], F.zero);
				return x;

			}

			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{

				linbox_check (A.coldim() == A.rowdim());
				linbox_check (A.coldim() == b.size());
				bool consistent = true;

				for (size_t i = A.getRank(); i < b.size(); ++i)
					if (!F.isZero (b[i]))
						consistent = false;
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				FFLAS::ftrsv ((typename Field::Father_t)F, FFLAS::FflasUpper, FFLAS::FflasNoTrans, FFLAS::FflasNonUnit,
					      A.getRank(), A.getPointer(), A.getStride(), &b[0], 1);

				return b;
			}

		}; // end of class FactorizedMatrixLeftUSolve

		template <class Field>
		class FactorizedMatrixRightUSolve<Field, std::vector<typename Field::Element> > {
		public:
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& x,
									   const std::vector<typename Field::Element>& b ) const
			{
				linbox_check (x.size() == A.rowdim());
				linbox_check (A.coldim() == b.size());
				typename Field::Element * bp = b.getPointer();
				typename Field::Element * xp = x.getPointer();
				size_t R = A.getRank();

				for (size_t i = 0; i < R; ++i)
					F.assign (x[i], b[i]);

				FFLAS::ftrsv ((typename Field::Father_t)F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit,
					      R, A.getPointer(), A.getStride(), xp, 1);

				bool consistent = true;
				if (b.size() > x.size()) {
					typename Field::Element* W = new typename Field::Element [b.size() - R];

					FFLAS::fgemv ((typename Field::Father_t)F, FFLAS::FflasTrans,
						      R, b.size() - R,
						      F.one, A.getPointer() + R, A.getStride, xp, 1,
						      F.zero, W, 1);

					for (size_t i = 0; i < b.size() - R; ++i)
						if (!F.areEqual (W[i], b[i + R]))
							consistent = false;
					delete[] W;
				}
				else {
					FFLAS::fgemv ((typename Field::Father_t)F, FFLAS::FflasTrans,
						      R, b.size() - R,
						      F.one, A.getPointer() + R, A.getStride, xp, 1,
						      F.zero, xp + R, 1);

					for (size_t i = R; i < b.size(); ++i)
						if (!F.areEqual (x[i], b[i]))
							consistent = false;
				}

				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				for (size_t j = R; j < x.size(); ++j)
					F.assign (x[j], F.zero);
				return x;
			}
			std::vector<typename Field::Element>& operator() ( const Field& F,
									   const LQUPMatrix<Field>& A,
									   std::vector<typename Field::Element>& b ) const
			{
				linbox_check (A.coldim() == A.rowdim());
				linbox_check (b.size() == A.rowdim());

				typename Field::Element * bp = &b[0];
				size_t R = A.getRank();

				FFLAS::ftrsv ((typename Field::Father_t)F, FFLAS::FflasUpper, FFLAS::FflasTrans, FFLAS::FflasNonUnit,
					      R, F.one, A.getPointer(), A.getStride(), bp, 1);

				FFLAS::fgemv ((typename Field::Father_t)F, FFLAS::FflasTrans,
					      R, b.size() - R,
					      F.one, A.getPointer() + R, A.getStride, bp, 1,
					      F.mOne, bp + R, 1);

				bool consistent = true;
				for (size_t j = R; j < b.size(); ++j)
					if (!F.isZero (b[j]))
						consistent = false;
				if (!consistent)
					throw  LinboxMathInconsistentSystem ("Linear system is inconsistent");

				return b;
			}
		}; // end of class FactorizedMatrixRightUSolve


	}
}

namespace LinBox
{
	template <class Field>
	LQUPMatrix<Field>::LQUPMatrix (const BlasMatrix<Field>& A) :
		_field(A.field()), _factLU(*(new BlasMatrix<Field> (A))) ,
		_permP(*(new BlasPermutation<size_t>(A.coldim()))),
		_permQ(*(new BlasPermutation<size_t>(A.rowdim()))),
		_m(A.rowdim()), _n(A.coldim()),
		_alloc(true),_plloc(true)
	{
		if (!A.coldim() || !A.rowdim()) {
			// throw LinBoxError("LQUP does not accept empty matrices");
			_rank = 0 ;
		}
		else {
			_rank= FFPACK::LUdivine((typename Field::Father_t) _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans,
						 _m, _n,
						 _factLU.getPointer(),_factLU.getStride(),
						 _permP.getWritePointer(), _permQ.getWritePointer(),
						 FFPACK::FfpackLQUP );
		}
		_permP.setOrder(_rank);
		_permQ.setOrder(_rank);

	}

	template <class Field>
	LQUPMatrix<Field>::LQUPMatrix (BlasMatrix<Field>& A) :
		_field(A.field()), _factLU(static_cast<BlasMatrix<Field>&> (A)) ,
		_permP(*(new BlasPermutation<size_t>(A.coldim()))),
		_permQ(*(new BlasPermutation<size_t>(A.rowdim()))),
		// _permP(A.coldim()), _permQ(A.rowdim()),
		_m(A.rowdim()), _n(A.coldim()),
		_alloc(false),_plloc(true)
	{
		if (!A.coldim() || !A.rowdim()) {
			// throw LinBoxError("LQUP does not accept empty matrices");
			_rank = 0 ;
		}
		else {
			_rank= FFPACK::LUdivine((typename Field::Father_t) _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans,
						 _m, _n,
						 _factLU.getPointer(),_factLU.getStride(),
						 _permP.getWritePointer(), _permQ.getWritePointer(),
						 FFPACK::FfpackLQUP );
		}
		_permP.setOrder(_rank);
		_permQ.setOrder(_rank);

	}

	template <class Field>
	LQUPMatrix<Field>::LQUPMatrix (const BlasMatrix<Field>& A,
				       BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) :
		_field(A.field()), _factLU(*(new BlasMatrix<Field> (A))) ,
		_permP(P), _permQ(Q),
		_m(A.rowdim()), _n(A.coldim()),
		_alloc(true),_plloc(false)
	{

		linbox_check(_permQ.getOrder()==A.rowdim());
		linbox_check(_permP.getOrder()==A.coldim());


		_rank= FFPACK::LUdivine((typename Field::Father_t) _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans,
					 _m, _n,
					 _factLU.getPointer(),_factLU.getStride(),
					 _permP.getWritePointer(), _permQ.getWritePointer(),
					 FFPACK::FfpackLQUP );

		_permP.setOrder(_rank);
		_permQ.setOrder(_rank);

	}

	template <class Field>
	LQUPMatrix<Field>::LQUPMatrix (BlasMatrix<Field>& A,
				       BlasPermutation<size_t> & P, BlasPermutation<size_t> & Q) :
		_field(A.field()), _factLU(static_cast<BlasMatrix<Field>&> (A)) ,
		_permP(P), _permQ(Q),
		_m(A.rowdim()), _n(A.coldim()),
		_alloc(false),_plloc(false)
	{


		linbox_check(_permQ.getOrder()<=A.rowdim());
		linbox_check(_permP.getOrder()<=A.coldim());
		if (_permQ.getOrder() == 0)
			_permQ.resize(A.rowdim());
		if (_permP.getOrder() == 0)
			_permP.resize(A.coldim());


		_rank= FFPACK::LUdivine((typename Field::Father_t) _field,FFLAS::FflasNonUnit,  FFLAS::FflasNoTrans,
					 _m, _n,
					 _factLU.getPointer(),_factLU.getStride(),
					 _permP.getWritePointer(), _permQ.getWritePointer(),
					 FFPACK::FfpackLQUP );
		_permP.setOrder(_rank);
		_permQ.setOrder(_rank);

	}

	template <class Field>
	LQUPMatrix<Field>::~LQUPMatrix ()
	{
		if (_alloc)
			delete &_factLU;
		if (_plloc) {
			delete &_permP;
			delete &_permQ;
		}
	}

	template <class Field>
	Field& LQUPMatrix<Field>::field()
	{
		return _field;
	}

	template <class Field>
	size_t LQUPMatrix<Field>::rowdim() const
	{
		return _m;
	}

	template <class Field>
	size_t LQUPMatrix<Field>::coldim() const
	{
		return _n;
	}

	template <class Field>
	size_t LQUPMatrix<Field>::getRank() const
	{
		return _rank;
	}

	template <class Field>
	const BlasPermutation<size_t>& LQUPMatrix<Field>::getP() const
	{
		return _permP;
	}

	template <class Field>
	BlasPermutation<size_t> & LQUPMatrix<Field>::getP( BlasPermutation<size_t> & P ) const
	{
		P = BlasPermutation<size_t>(_permP.getPointer(),_rank);
		return P;
	}

	template <class Field>
	const BlasPermutation<size_t>& LQUPMatrix<Field>::getQ() const
	{
		return _permQ;
	}

	template <class Field>
	BlasPermutation<size_t> & LQUPMatrix<Field>::getQ( BlasPermutation<size_t> & Qt ) const
	{
		Qt = BlasPermutation<size_t>(_permQ.getPointer(),_rank);
		return Qt ;
	}


	// get the Matrix L
	template <class Field>
	inline TriangularBlasMatrix<Field>&
	LQUPMatrix<Field>::getL( TriangularBlasMatrix<Field>& L,
				 bool _QLUP ) const
	{

		linbox_check( L.coldim() == _m);
		linbox_check( L.rowdim() == _m);
		linbox_check( L.getUpLo() == LinBoxTag::Lower);
		linbox_check( L.getDiag() == LinBoxTag::Unit);

#if 0
		if (_m > _n) {
			size_t i = 0 ;
			for ( ; i<_n; ++i ){
				size_t j=0;
				for (; j<i ; ++j )
					L.setEntry( i, j, _factLU.getEntry(i,j) );
				for (; j<_m; ++j )
					L.setEntry( i, j, zero );
			}
			for ( ; i<_m; ++i ){
				size_t j=0;
				for (; j<_n ; ++j )
					L.setEntry( i, j, _factLU.getEntry(i,j) );
				for (; j<_m; ++j )
					L.setEntry( i, j, zero );
			}


		}
		else {
			for ( size_t i=0; i<_m; ++i ){
				size_t j=0;
				for (; j<i ; ++j )
					L.setEntry( i, j, _factLU.getEntry(i,j) );
				for (; j<_m; ++j )
					L.setEntry( i, j, zero );
			}
		}
#endif
#if 1 /*  slower */
		for ( size_t i=0; i<_m; ++i ){
			size_t j=0;
			for (; j< ((i<_n)?i:_n) ; ++j )
				L.setEntry( i, j, _factLU.getEntry(i,j) );
			for (; j<_m; ++j )
				L.setEntry( i, j, _field.zero );
		}
#endif

		if (!_permQ.isIdentity())
			FFPACK::applyP((typename Field::Father_t) _field, FFLAS::FflasRight, FFLAS::FflasNoTrans,
					_m,0,(int)_permQ.getOrder(),
					L.getWritePointer(), _m, _permQ.getPointer() );
		for ( size_t i=0; i<_m; ++i )
			L.setEntry( i, i, _field.one );
		if (_QLUP) {
			if (!_permQ.isIdentity()) {
				FFPACK::applyP((typename Field::Father_t) _field, FFLAS::FflasLeft, FFLAS::FflasNoTrans,
						_m,0,(int)_permQ.getOrder(),
						L.getWritePointer(), _m, _permQ.getPointer() );
				FFPACK::applyP((typename Field::Father_t) _field, FFLAS::FflasRight, FFLAS::FflasTrans,
						_m,0,(int)_permQ.getOrder(),
						L.getWritePointer(), _m, _permQ.getPointer() );

			}
		}

		return L;

	}

	// get the matrix U
	template <class Field>
	inline TriangularBlasMatrix<Field>&
	LQUPMatrix<Field>::getU( TriangularBlasMatrix<Field>& U ) const
	{

		linbox_check( U.rowdim() == _m);
		linbox_check( U.coldim() == _n);
		linbox_check( U.getUpLo() == LinBoxTag::Upper);
		linbox_check( U.getDiag() == LinBoxTag::NonUnit);
		for ( size_t i=0; i<_m; ++i )
			for ( size_t j=i; j<_n; ++j )
				U.setEntry( i, j, _factLU.getEntry(i,j) );
		return U;
	}

	// get the Matrix S (from the LSP factorization of A deduced from LQUP)
	template <class Field>
	inline BlasMatrix<Field>&
	LQUPMatrix<Field>::getS( BlasMatrix<Field>& S) const
	{

		linbox_check( S.rowdim() == _m);
		linbox_check( S.coldim() == _n);
		for ( size_t i=0; i<_m; ++i ){
			for ( size_t j=0; j<i; ++j )
				S.setEntry( i, j, _field.zero );
			for ( size_t j=i; j<_n; ++j )
				S.setEntry( i, j, _factLU.getEntry(i,j) );
		}

		if (!_permQ.isIdentity())
			FFPACK::applyP((typename Field::Father_t) _field, FFLAS::FflasLeft, FFLAS::FflasTrans,
					_n, 0,(int) _permQ.getOrder(),
					S.getWritePointer(), _n, _permQ.getPointer() );
		return S;
	}

	template <class Field>
	typename Field::Element* LQUPMatrix<Field>::getPointer() const
	{
		return _factLU.getPointer();
	}

	/*! @internal get  the stride in \c _factLU
	*/
	template <class Field>
	size_t LQUPMatrix<Field>::getStride() const
	{
		return _factLU.getStride();
	}

	// solve AX=B
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::left_solve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixLeftSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve AX=B (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::left_solve(Operand& B) const
	{
		return Protected::FactorizedMatrixLeftSolve<Field,Operand>()( _field, *this, B );
	}

	// solve XA=B
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_solve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixRightSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve XA=B (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_solve(Operand& B) const
	{
		return Protected::FactorizedMatrixRightSolve<Field,Operand>()( _field, *this, B );
	}

	// solve LX=B (L from LQUP)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::left_Lsolve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixLeftLSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve LX=B (L from LQUP) (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::left_Lsolve(Operand& B) const
	{
		return Protected::FactorizedMatrixLeftLSolve<Field,Operand>()( _field, *this, B );
	}

	// solve XL=B (L from LQUP)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_Lsolve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixRightLSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve XL=B (L from LQUP) (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_Lsolve(Operand& B) const
	{
		return Protected::FactorizedMatrixRightLSolve<Field,Operand>()( _field, *this, B );
	}

	// solve UX=B (U from LQUP is r by r)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::left_Usolve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixLeftUSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve UX=B (U from LQUP) (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::rleft_Usolve(Operand& B) const
	{
		return Protected::FactorizedMatrixLeftUSolve<Field,Operand>()( _field, *this, B );
	}

	// solve XU=B (U from LQUP)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_Usolve(Operand& X, const Operand& B) const
	{
		return Protected::FactorizedMatrixRightUSolve<Field,Operand>()( _field, *this, X, B );
	}

	// solve XU=B (U from LQUP) (X is stored in B)
	template <class Field>
	template <class Operand>
	Operand& LQUPMatrix<Field>::right_Usolve(Operand& B) const
	{
		return Protected::FactorizedMatrixRightUSolve<Field,Operand>()( _field, *this, B );
	}


} //end of namespace LinBox


#endif // __LINBOX_factorized_matrix_INL


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

