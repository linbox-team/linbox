#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_conversion_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_conversion_INL

namespace LinBox
{
	template <class _FieldGF, class _Storage, class _MatrixElement, class _FieldF>
	PolynomialMatrix<PMType::polfirst,PMStorage::plain,_FieldF>&
		SlicedPolynomialMatrixtoPolynomialMatrix<_FieldGF, _Storage,
		_MatrixElement, _FieldF>::operator()
			(PolynomialMatrix<PMType::polfirst, PMStorage::plain,_FieldF>& PM,
			SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement> &SPM)
	{
		size_t rd = SPM.rowdim();
		size_t cd = SPM.coldim();
		size_t l = SPM.length();
		for (size_t i = 0; i < rd; i++)
		{
			for (size_t j = 0; j < cd; j++)
			{
				for (size_t k = 0; k < l; k++)
				{
					PM.ref(i, j, k) = SPM.getEntry(k, i, j);
				}
			}
		}
		return PM;
	}

	template <class _FieldF, class _FieldGF, class _Storage, class _MatrixElement>
	SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement>&
		PolynomialMatrixtoSlicedPolynomialMatrix<_FieldF, _FieldGF,
		_Storage, _MatrixElement>::operator()
			(SlicedPolynomialMatrix<_FieldGF, _Storage, _MatrixElement> &SPM,
			PolynomialMatrix<PMType::polfirst,PMStorage::plain,_FieldF> &PM)
	{
		size_t rd = SPM.rowdim();
		size_t cd = SPM.coldim();
		size_t l = SPM.length();
		for (size_t k = 0; k < l; k++)
		{
			SPM.setMatrixCoefficient(k, PM[k]);
		}
		return SPM;
	}

	template <class _FieldGF, class _Storage1, class _MatrixElement, class TT, class _Storage2>
	BlasMatrix<Givaro::GFqDom<TT>, _Storage2>&
		SlicedPolynomialMatrixtoBlasMatrix<_FieldGF, _Storage1,
		_MatrixElement, TT, _Storage2>::operator()
			(BlasMatrix<Givaro::GFqDom<TT>, _Storage2> &BM,
			SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement> &SPM)
	{
		size_t rd = SPM.rowdim();
		size_t cd = SPM.coldim();
		size_t l = SPM.length();
		GFqDom<TT> GF(SPM.fieldGF().characteristic(), SPM.fieldGF().exponent());
		GFqDom<TT>::Element el;
		PolynomialMatrix<PMType::polfirst,PMStorage::plain, spm::IntField> PM(
			SPM.fieldF(), rd, cd, l);
		SlicedPolynomialMatrixtoPolynomialMatrix<_FieldGF, _Storage1,
		_MatrixElement, spm::IntField>()(PM, SPM);
		for (size_t i = 0; i < rd; i++)
		{
			for (size_t j = 0; j < cd; j++)
			{
				GF.init(el, PM(i, j));
				BM.setEntry(i, j, el);
			}
		}
		return BM;
	}

	template <class TT, class _Storage2, class _FieldGF, class _Storage1, class _MatrixElement>
	SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement>&
		BlasMatrixtoSlicedPolynomialMatrix<TT, _Storage2,
		_FieldGF, _Storage1, _MatrixElement>::operator()
			(SlicedPolynomialMatrix<_FieldGF, _Storage1, _MatrixElement> &SPM,
			BlasMatrix<Givaro::GFqDom<TT>, _Storage2> &BM)
	{
		size_t rd = SPM.rowdim();
		size_t cd = SPM.coldim();
		size_t l = SPM.length();
		if (l == 1)
		{
			for (size_t i = 0; i < rd; i++)
			{
				for (size_t j = 0; j < cd; j++)
				{
					SPM.setEntry(0, i, j, BM.getEntry(i, j));
				}
			}
			return SPM;
		}
		else
		{
			size_t p = SPM.fieldGF().characteristic();
			for (size_t i = 0; i < rd; i++)
			{
				for (size_t j = 0; j < cd; j++)
				{
					size_t el = log2pol(BM.getEntry(i, j));
					for (size_t k = 0; k < l; k++)
					{
						SPM.setEntry(k, i, j, el % p);
						el /= p;
					}
				}
			}
			return SPM;
		}
	}
} // LinBox

#endif