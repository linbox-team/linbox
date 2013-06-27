
#ifndef __LINBOX_algorithms_iml_row_echelon_transform_INL
#define __LINBOX_algorithms_iml_row_echelon_transform_INL

namespace LinBox {
namespace iml {


template<class FField>
RowEchelonTransform<FField>::RowEchelonTransform(const BlasMatrix<FField> & A)  :
	_field(A.field())
	, _echlA(A)
	, _permQ(A.rowdim()+1)
	, _rowPf(A.rowdim()+1)
	,_det(0)
	,_red(false)
{
	_echlA.copy(A);
}

template<class FField>
RowEchelonTransform<FField>::RowEchelonTransform(BlasMatrix<FField> & A)  :
	_field(A.field())
	, _echlA(A)
	, _permQ(A.rowdim()+1)
	, _rowPf(A.rowdim()+1)
	,_det(0)
	,_red(false)
{
}


template<class FField>
size_t RowEchelonTransform<FField>::reduce_rec( BlasMatrix<FField> & A, size_t m1, size_t m2, size_t k,
		const size_t ks, size_t frows, size_t lrows, size_t redflag,
		size_t eterm,
		std::vector<size_t> & P, std::vector<size_t>&rp,
		Element & d)
		{
	typedef BlasMatrix<FField> Matrix ;
			size_t m = A.coldim();
			size_t n = A.rowdim();
			size_t i, j, r1, r2, r, ri, mm, inv;
			Element a;
			double t, b;
			BlasMatrixDomain<FField> BMD(_field);

			if (m1 == m2)
			{
#if 0
				for (i = k+1; i <= n; i++)
					if (*(A+(i-1)*m+m1-1) != 0)
						break;
#endif
				for (i = k+1; i <= n; i++)
					if (!_field.isZero(A.refEntry(i-1,m1-1)))
					    break;

				if ((i > n) && (eterm == 0)) {
					return 0;
				}
				else if ((i > n) && (eterm == 1)) {
					d = 0; //! @bug inits an integer ?
					return 0;
				}
				if (i > k+1) {
					FFLAS::fswap(_field,m-m1+1,&A.refEntry(k,m1-1),1,&A.refEntry(i-1,m1-1),1);
				}
				if (k-ks > 0){
					FFLAS::fswap(_field,k-ks,&A.refEntry(k,0),1,&A.refEntry(i-1,0),1);
				}
				P[k+1] = i;
				//!@todo what if no inversion ?
				_field.inv(a,A.getEntry(k,m1-1));

				/* mp_a is not relatively prime with modulus */
				// if (!inv) {
					// throw("in RowEchelonTransform: modulus is composite");
				// }
				_field.init(b,A.getEntry(k,m1-1));
				if (b < 0) { b = (double)_field.characteristic()+b; } //! @bug assert non balanced ?
				if ((frows == 1) && (redflag == 1)) {
					for (j = 1; j <= n; j++){
						_field.init(A.refEntry(j-1,k-ks),A.getEntry(j-1,m1-1)*(-a));
					}
					A.setEntry(k,k-ks,a);
				}
				else {
					if (k+2 <= n) {
						for (j = k+2; j <= n; j++) {
							_field.init(A.refEntry(j-1,k-ks),A.getEntry(j-1,m1-1)*(-a));
						}
					}
					if (k > 0) {
						for (j = 1; j <= k; j++){
							A.setEntry(j-1,k-ks,_field.zero);
						}
					}
					A.setEntry(k,k-ks,_field.one);
				}
				++rp[0];
				_field.init(d,d*b);
				ri =  rp[0];
				rp[ri] = m1;
				return 1;
			}
			/* Recursively solve the first subproblem */
			mm = m1+((m2-m1)/2);
			r1 = reduce_rec(A,m1, mm, k, ks, frows, 1, redflag,\
						     eterm, P, rp, d);
			if ((eterm == 1) && (r1 < mm-m1+1)) {
				d= 0;
				return 0;
			}
			/* If r1=0 then don't need to construct second subproblem */
			if (r1 > 0) {
				/* Compute U1.A2 by submatrix multiply */
				if (k+r1 < n) {
					BlasSubmatrix<Matrix> U1(A,k+r1,k-ks,n-k-r1,r1);
					BlasSubmatrix<Matrix> A1(A,k,mm,r1,m2-mm);
					BlasSubmatrix<Matrix> A2(A,k+r1,mm,n-k-r1,m2-mm);
					BMD.axpyin(A2,U1,A1);

				}
				if ((frows == 1) && (redflag == 1)) {
					if (k > 0) {
						BlasSubmatrix<Matrix> U1(A,0,k-ks,k,r1);
						BlasSubmatrix<Matrix> A1(A,k,mm,r1,m2-mm);
						BlasSubmatrix<Matrix> A2(A,0,mm,k,m2-mm);
						BMD.axpyin(A2,U1,A1);

					}
					BlasSubmatrix<Matrix> U1(A,k,k-ks,r1,r1);
					BlasMatrix<FField> A1(A,k,mm,r1,m2-mm);
					BlasSubmatrix<Matrix> A2(A,k,mm,r1,m2-mm);
					BMD.mul(A2,U1,A1);

				}
			}
			/* Recursively solve the second subproblem */
			r2 = reduce_rec(A,mm+1, m2, k+r1, ks, frows, lrows, \
						     redflag, eterm, P, rp, d);
			r = r1+r2;
			if ((eterm == 1) && (r < m2-m1+1)) {
				d = 0;
				return 0;
			}
			/* Only need to combine when both subproblems nontrivial */
			if ((r2 > 0) && (r1 > 0)) {
				if ((k+r+1 <= n) && (lrows == 1)) {
					/* Bottom block of U */
					BlasSubmatrix<Matrix> U1(A,k+r,k-ks+r1,n-k-r,r-r1);
					BlasSubmatrix<Matrix> A1(A,k+r1,k-ks,r-r1,r1);
					BlasSubmatrix<Matrix> A2(A,k+r,k-ks,n-k-r,r1);
					BMD.axpyin(A2,U1,A1);

				}
				if (frows == 1) {
					if (redflag == 1)
						i = 1;
					else
						i = k+1;

					/* Rows i..k of top block of U plus first r1 rows of middle block */
					{
						BlasSubmatrix<Matrix> U1(A,i-1,k-ks+r1,k+r1-i+1,r-r1);
						BlasSubmatrix<Matrix> A1(A,k+r1,k-ks,r-r1,r1);
						BlasSubmatrix<Matrix> A2(A,i-1,k-ks,k+r1-i+1,r1);
						BMD.axpyin(A2,U1,A1);
					}



					/* Last r2 rows of middle block */
					{
						BlasSubmatrix<Matrix> U1(A,k+r1,k-ks+r1,r-r1,r-r1);
						BlasMatrix<FField>    A1(A,k+r1,k-ks,r-r1,r1);
						BlasSubmatrix<Matrix> A2(A,k+r1,k-ks,r-r1,r1);
						BMD.mul(A2,U1,A1);
					}


				}
			}
			return r;
		}

} // iml
} // LinBox



#endif // __LINBOX_algorithms_iml_row_echelon_transform_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
