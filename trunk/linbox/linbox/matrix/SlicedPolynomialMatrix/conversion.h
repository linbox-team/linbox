#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_conversion_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_conversion_H

#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/gfq.h>

namespace LinBox
{
	//FieldF == SlicedPolynomialMatrix::IntField
	template <class _FieldGF, class _Storage, class _MatrixElement = double, class _FieldF>
	class SlicedPolynomialMatrixtoPolynomialMatrix
	{
	public:
		PolynomialMatrix<PMType::polfirst,PMStorage::plain,_FieldF>
			&operator() (PolynomialMatrix<PMType::polfirst,PMStorage::plain,_FieldF> &PM,
			SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement> &SPM);
	};

	//FieldF == SlicedPolynomialMatrix::IntField
	template <class _FieldF, class _FieldGF, class _Storage, class _MatrixElement = double>
	class PolynomialMatrixtoSlicedPolynomialMatrix
	{
	public:
		SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement>
			&operator() (SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement> &SPM),
			PolynomialMatrix<PMType::polfirst,PMStorage::plain,_FieldF> &PM);
	};

	template <class _FieldGF, class _Storage1, class _MatrixElement = double, class TT, class _Storage2>
	class SlicedPolynomialMatrixtoBlasMatrix
	{
	private:
		typedef SlicedPolynomialMatrixtoBlasMatrix<_FieldGF, _Storage1,
		_MatrixElement, TT, _Storage2> spm;
	public:
		BlasMatrix<Givaro::GFqDom<TT>, _Storage2>
			&operator() (BlasMatrix<Givaro::GFqDom<TT>, _Storage2> &BM,
			SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement> &SPM);
	};

	template <class TT, class _Storage2, class _FieldGF, class _Storage1, class _MatrixElement = double>
	class BlasMatrixtoSlicedPolynomialMatrix
	{
	public:
		SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement>
			&operator() (SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement> &SPM,
			BlasMatrix<Givaro::GFqDom<TT>, _Storage2> &BM);
	};
} /* end of namespace LinBox */

#include "conversion.inl"

#endif