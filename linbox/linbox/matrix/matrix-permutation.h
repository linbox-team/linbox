/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/** @file matrix/matrix-permutation.h
 * @ingroup matrix
 * A permutation class for operations on permutations, their representations
 * and matrix row/column permuting.
 *
 * We provide a \ref LinBox::PackedPermutation class that stores the
 * permutation packed in a Lapack style and a \ref LinBox::MatrixPermutation
 * class that represents a permutation naturally.  Converstions are provided.
 */

#ifndef __LINBOX_matrix_permutation_H
#define __LINBOX_matrix_permutation_H

#include <vector>
#include <ostream>


// PackedPermutation
namespace LinBox
{

	template<class _Uint>
	class PermutationInterface {
	public :
		virtual _Uint getSize() = 0;
	};

	// forward declaration
	template<class _Uint>
	class MatrixPermutation ;

	template<class _Uint>
	class PackedPermutation ;

	/** Lapack-style permutation.
	 * @ingroup permutation
	 *
	 * The permutation \f$P=[a,b,c]\f$ tells that \f$1\f$ and \f$a\f$ are transposed, then
	 * \f$2\f$ and \f$b\f$, and finally \f$3\f$ and \f$c\f$.
	 * It is always true that \f$P[i] \geq i\f$.
	 * Example : if <code>P=[1,4,4]</code> and <code>V=[1,2,3,4,5]</code>, then <code>P.V=[1,4,2,3,5]</code>.
	 * @pre if Q_ is built, then P_=Q_
	*/
	template<class _UnsignedInt> // unsigned * ou Integer
	class PackedPermutation : PermutationInterface<_UnsignedInt> {
		typedef PackedPermutation<_UnsignedInt> BlasPerm ;
	public :
		PackedPermutation() ;
		~PackedPermutation() ;

		PackedPermutation(const _UnsignedInt * V, const _UnsignedInt & n, bool clean = false) ;
		PackedPermutation(const std::vector<_UnsignedInt> & V, bool clean = false);
		PackedPermutation(const MatrixPermutation<_UnsignedInt> &M, bool clean = false);
#if 0
		PackedPermutation(MatrixPermutation & P);
		PackedPermutation(TranspositionPermutation & P);

		//        void Invert() ;
		//        BlasPerm & getInverse(BlasPerm & P) ;


		template<class _Perm>
		_Perm & convert(_Perm & NewP);

		template<class _Perm>
		_Perm & convertTranspose(_Perm & NewP);

		void getOrder() ;
		void resize() ;
		void size() ;


		void Tranpsose() ;
		BlasPerm & getTranspose(BlasPerm & P) ;
		void Tranpsose(Index &i, Index & j);

		void applyP( BlasPerm & P, enum Side s, enum Trans = NoTranspose) ;
		void applyPT( BlasPerm & P, enum Side s) ;
#endif

		// operator =
		PackedPermutation<_UnsignedInt>& operator= (const PackedPermutation<_UnsignedInt> & P)
		{
			r_       = P.r_;
			n_       = P.n_;
			P_       = P.P_;
			Q_       = P.Q_;
			cleaned_ = P.cleaned_;

			return (*this) ;
		}

		/*  size */
		_UnsignedInt getSize() ;
		_UnsignedInt getOrder() ;

		_UnsignedInt getCompression() ;
		std::vector<_UnsignedInt> & getStorage() ;
		void resize(_UnsignedInt & s, bool with_zeros=true) ;

		template<class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x)  ;
		template<class OutVector, class InVector>
		OutVector &applyTranspose (OutVector &y, const InVector &x) ;



		/*  convert */
		/** converts a \c PackedPermutation to a \c MatrixPermutation */
		MatrixPermutation<_UnsignedInt> & Convert(const MatrixPermutation<_UnsignedInt> &P);

		/*  apply */
		/*! \f$ M \gets P M\f$   */
		template<class Matrix>
		Matrix & applyRows(Matrix &M);
		/*! \f$ M \gets M P\f$   */
		template<class Matrix>
		Matrix & applyCols(Matrix &M);

		/*! \f$ M \gets M P^t\f$   */
		template<class Matrix>
		Matrix & applyTransposeRows(Matrix &M);
		/*! \f$ M \gets P^t M\f$   */
		template<class Matrix>
		Matrix & applyTransposeCols(Matrix &M);

		//_UnsignedInt & operator[] (const _UnsignedInt &i) ;
		_UnsignedInt & operator[] (const _UnsignedInt i) const ;

		/*! col \p i and col \p j are swapped
		 */
		void TransposeCols(_UnsignedInt i, _UnsignedInt j);

		/*! row \p i and row \p j are swapped
		 */
		void TransposeRows(_UnsignedInt i, _UnsignedInt j);

		/*  invert */
		void Transpose();
		void Invert();
		BlasPerm & Transpose(BlasPerm &Mt);
		BlasPerm & Invert(BlasPerm &Mt);

		/*  print */
		/*! writes on output stream \p o */
		std::ostream & write(std::ostream & o, bool Lapack=true) ;

		/*! writes \p P on output stream \p o */
		template<class _Uint>
		friend std::ostream & operator<<(std::ostream &o, BlasPerm & P) ;


	protected :
		_UnsignedInt			r_ ;	// size of compressed permutation
		_UnsignedInt			n_ ;	// dim of permutation
		std::vector<_UnsignedInt>	P_ ;	// blas permutation
		mutable std::vector<_UnsignedInt>       Q_ ;    // corresponding matrix permutation
		mutable bool                            cleaned_ ;

		// hmmmm...
		// using stl vectors instead of pointers for the sake of simplicity...
		// this allows permutation up to MAX_INT size. Not so restricting for now...

		void Compress_() ;
		void BuildQ_() const ;
		void InvertQ_();
		std::vector<_UnsignedInt> &InvertQ_(std::vector<_UnsignedInt> & Qinv);
		void BuildP_(std::vector<_UnsignedInt>&Q, std::vector<_UnsignedInt>&Qinv);
		bool CheckP_();
		void InitQ_() ;


	};
} // LinBox

// MatrixPermutation
namespace LinBox
{

	/*! Permutation classique.
	 * @ingroup permutation
	 */
	template<class _UnsignedInt>
	class MatrixPermutation  : PermutationInterface<_UnsignedInt> {
		typedef MatrixPermutation<_UnsignedInt> MatPerm ;
	private :
		_UnsignedInt			n_ ; // order of permutation
		std::vector<_UnsignedInt>	P_ ; // _M_[i] = j ssi P(i) = j

	public :
		MatrixPermutation();
		~MatrixPermutation() {};
		MatrixPermutation(const _UnsignedInt * V, const _UnsignedInt & n) ;
		MatrixPermutation(const std::vector<_UnsignedInt> & V) ;

		_UnsignedInt & operator[] (const _UnsignedInt i) const ;
		_UnsignedInt getSize() const ;
		_UnsignedInt getSize() ;

		void resize( _UnsignedInt n ) ;

		void Transpose();
		void Invert();
		MatPerm & Transpose(MatPerm &Mt);
		MatPerm & Invert(MatPerm &Mt);



		/*  print */
		/*! writes on output stream \p o */
		std::ostream & write(std::ostream & o) ;

		/*! writes \p P on output stream \p o */
		template<class _Uint>
		friend std::ostream & operator<<(std::ostream &o, MatPerm & P) ;

		template<class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x) const ;
		template<class OutVector, class InVector>
		OutVector &applyTranspose (OutVector &y, const InVector &x) const ;


		void TransposeRows(_UnsignedInt i, _UnsignedInt j);
		void TransposeCols(_UnsignedInt i, _UnsignedInt j);

	};

} // LinBox

#if 0 /* stuff to be removed */
namespace LinBox
{

	//! produit de permutations
	template<class _UnsignedInt>
	class TranspositionPermutation {
		typedef std::pair<_UnsignedInt,_UnsignedInt> Transposition  ;
	private :
		_UnsignedInt _n_ ;			// order of permutation
		_UnsignedInt _r_ ;			//  number of transpositions
		std::vector<Transposition> _T_ ;	// if _T_[k] = (i,j) then P(i) = j

	public :
	};

	//!@todo produit de cycles à supports disjoints


	//! on regroupe dedans les parties "BB" des perms précédentes.

	template<class _Field>
	class BlackboxInterface {
	public :
		typedef BlackboxInterface<_Field> Self_t ;
		virtual size_t rowdim() const ;
		virtual size_t coldim() const ;
		template<class OutVec,InVec>
		virtual OutVec &apply (OutVec & , const InVec & ) const ;
		template<class OutVec,InVec>
		virtual OutVec &applyTranspose (OutVec & , const InVec & ) const ;
		virtual _Field & field() const ;
		template<class _Tp1>
		virtual struct rebind {
			typedef BlackboxInterface<_Tp1> other ;
			virtual void operator() (other &, cosnt Self_t &, const _Tp1 &) ;
		} ;
	};

	template<class _Perm, class _Field>
	class BlackBoxPermutation : public BlackboxInterface { // ?????????????
	private :
		_Perm  _P_ ;
		_Field _F_ ; //????????????????????????
	public :
		size_t rowdim() ;
		size_t coldim() ;
		void permute( Index i, Index j);
		const Field & field() { return _F_ ;}
		apply();
		applyTranspose() ;
		template<typename _Tp1>
		struct rebind {} ;

		setStorage();
		getStorage();

	};

#if 0
	template<class _Perm>
	class PermutationDomain {
	private :
		_Perm _P_ ;
	public :
		// P.permute(A) ??? :)
		template<class _Matrix>
		_Matrix & permute (_Matrix & M,  enum Side s, enum Trans t = NoTranspose) ;
		template<class _Matrix>
		_Matrix & invmute (_Matrix & M, enum Side s) ;



	};
#endif
}
#endif

#include "matrix-permutation.inl"

#endif //__LINBOX_matrix_permutation_H
