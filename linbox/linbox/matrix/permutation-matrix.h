
/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
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

/** @file matrix/permutation-matrix.h
 * @ingroup matrix
 * A permutation class for operations on permutations, their representations
 * and matrix row/column permuting.
 *
 * We provide a \ref LinBox::BlasPermutation class that stores the
 * permutation packed in a Lapack style and a \ref LinBox::MatrixPermutation
 * class that represents a permutation naturally.  Converstions are provided.
 */

#ifndef __LINBOX_matrix_permutation_H
#define __LINBOX_matrix_permutation_H

#include <vector>
#include <ostream>


// BlasPermutation
namespace LinBox
{


	// forward declaration
	template<class _UnsignedInt>
	class MatrixPermutation ;

	template<class _UnsignedInt>
	class BlasPermutation ;

	/** Lapack-style permutation.
	 * @ingroup permutation
	 *
	 * A Lapack permutation is represented with a vector \f$[p_1,p_2,\cdots, p_r]\f$
	 * such that \f$p_i > i\f$. Converting it to a classic representation
	 * of a permutation corresponds to taking an identity permutation and
	 * then successively permuting \f$(i,p_i)\f$.
	 * Example : if <code>P=[1,4,4]</code> and <code>V=[1,2,3,4,5]</code>, then <code>P.V=[1,4,2,3,5]</code>.
	 * @internal
	 * @pre if \c Q_ is built, then \c P_=Q_
	*/
	template<class _UnsignedInt> // unsigned * ou Integer
	class BlasPermutation /*  : PermutationInterface<_UnsignedInt> */ {
		typedef BlasPermutation<_UnsignedInt> BlasPerm ;
	public :
		BlasPermutation() ;
		~BlasPermutation() ;

		BlasPermutation(size_t n) ;

		BlasPermutation(const _UnsignedInt * V, const _UnsignedInt & n) ;
		BlasPermutation(const std::vector<_UnsignedInt> & V);
		BlasPermutation(const MatrixPermutation<_UnsignedInt> &M);


		//! copy operator (with copy)
		BlasPermutation<_UnsignedInt>& operator= (const BlasPermutation<_UnsignedInt> & P)
		{
			r_       = P.r_;
			n_       = P.n_;
			P_       = P.P_;
			Q_       = P.Q_;
			inv_     = P.inv_ ;

			return (*this) ;
		}

		/*! Returns the size of the permuation.
		 * If given, we return \p n as we see \p P in \f$S_n\f$.
		 * We default to the order of the permutation (minimal such \p n)
		 */
		_UnsignedInt getSize() const ;
		// _UnsignedInt getOrder() ;

		/*! Returns the order of the permuation */
		_UnsignedInt getOrder() const ;
		void setOrder( size_t r)  ;

		//! returns a copy of the raw storage.
		std::vector<_UnsignedInt>  getStorage() const
		{
			return P_;
		};

		// resize a blas permutation.
		void resize(_UnsignedInt s, bool with_zeros=true)
		{
			if (s < r_) {
				r_ = s ;
#ifndef NDEBUG
				std::cout << "*** Warning *** you are resizing a Blas Permutation (possibly corrupting it)" << std::endl;
#endif
			}
				n_ = s ;
				P_.resize(s);
				if (Q_.size())
					Q_.resize(s);
		}

		// template<class OutVector, class InVector>
		// OutVector &apply (OutVector &y, const InVector &x)  ;
		// template<class OutVector, class InVector>
		// OutVector &applyTranspose (OutVector &y, const InVector &x) ;

		/*  properties */
		bool isIdentity() const
		{
			return (!r_);
		}

		/*  convert */
		/*! Converts a \c BlasPermutation to a \c MatrixPermutation.
		 * @param[out] P MatrixPermutation to be created. Need not be initialized.
		 */
		MatrixPermutation<_UnsignedInt> & Convert(MatrixPermutation<_UnsignedInt> &P);

		/*  apply */
		// /*! \f$ M \gets P M\f$   */
		// template<class Matrix>
		// Matrix & applyRows(Matrix &M);
		// /*! \f$ M \gets M P\f$   */
		// template<class Matrix>
		// Matrix & applyCols(Matrix &M);

		// /*! \f$ M \gets M P^t\f$   */
		// template<class Matrix>
		// Matrix & applyTransposeRows(Matrix &M);
		// /*! \f$ M \gets P^t M\f$   */
		// template<class Matrix>
		// Matrix & applyTransposeCols(Matrix &M);

		//_UnsignedInt & operator[] (const _UnsignedInt &i) ;
		_UnsignedInt  operator[] (const _UnsignedInt i) const ;

		// /*! col \p i and col \p j are swapped
		 // */
		// void TransposeCols(_UnsignedInt i, _UnsignedInt j);

		// /*! row \p i and row \p j are swapped
		 // */
		// void TransposeRows(_UnsignedInt i, _UnsignedInt j);

		const _UnsignedInt* getPointer() const
		{
			linbox_check(P_.size());
			return &P_[0];
		}

		_UnsignedInt* getWritePointer()
		{
			linbox_check(P_.size());
			return &P_[0];
		}

		/*  invert */
		void Transpose();
		void Invert();
		BlasPerm & Transpose(BlasPerm &Mt);
		BlasPerm & Invert(BlasPerm &Mt);

		/* clean */
		void Compress() ;

		/*  print */
		/*! writes on output stream \p o */
		std::ostream & write(std::ostream & o, bool Lapack=true) const ;

		/*! writes \p P on output stream \p o */
		template<class _Uint>
		friend std::ostream & operator<<(std::ostream &o, BlasPerm & P) ;


	protected :
		_UnsignedInt			        r_ ;	// size of compressed permutation
		mutable _UnsignedInt			n_ ;	// dim of permutation
		std::vector<_UnsignedInt>	        P_ ;	// internal blas permutation
		mutable std::vector<_UnsignedInt>       Q_ ;    // corresponding matrix permutation
		bool                                    inv_ ;  // matrix is inverted ?

		void BuildQ_() const ;
		void InvertQ_();
		std::vector<_UnsignedInt> &InvertQ_(std::vector<_UnsignedInt> & Qinv);
		void BuildP_(std::vector<_UnsignedInt>&Q, std::vector<_UnsignedInt>&Qinv);
		bool CheckP_();
		void InitQ_() const ;


	};
} // LinBox

// MatrixPermutation
namespace LinBox
{

	/*! Permutation classique.
	 * @ingroup permutation
	 */
	template<class _UnsignedInt>
	class MatrixPermutation /*  : PermutationInterface<_UnsignedInt> */ {
		typedef MatrixPermutation<_UnsignedInt> Self_t ;
	private :
		_UnsignedInt			n_ ; // order of permutation
		std::vector<_UnsignedInt>	P_ ; // _M_[i] = j ssi P(i) = j

	public :
		MatrixPermutation();
		~MatrixPermutation() {};
		MatrixPermutation(const _UnsignedInt * V, const _UnsignedInt & n) ;
		MatrixPermutation(const std::vector<_UnsignedInt> & V) ;

		_UnsignedInt  operator[] (const _UnsignedInt i) const ;
		_UnsignedInt getSize() const ;
		// _UnsignedInt getSize() ;

		void resize( _UnsignedInt n ) ;

		void Transpose();
		void Invert();
		Self_t & Transpose(Self_t &Mt);
		Self_t & Invert(Self_t &Mt);



		/*  print */
		/*! writes on output stream \p o */
		std::ostream & write(std::ostream & o) const ;

		/*! writes \p P on output stream \p o */
		template<class _Uint>
		friend std::ostream & operator<<(std::ostream &o, Self_t & P) ;

		template<class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x) const ;
		template<class OutVector, class InVector>
		OutVector &applyTranspose (OutVector &y, const InVector &x) const ;


		void TransposeRows(_UnsignedInt i, _UnsignedInt j);
		void TransposeCols(_UnsignedInt i, _UnsignedInt j);

	};

} // LinBox


#include "permutation-matrix.inl"

#endif //__LINBOX_matrix_permutation_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

