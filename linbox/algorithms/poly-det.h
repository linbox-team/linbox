#ifndef __LINBOX_POLY_DET_H
#define __LINBOX_POLY_DET_H
/* by Alex Stachnik
*/

#include <givaro/extension.h>
#include <linbox/algorithms/poly-interpolation.h>
#include <linbox/solutions/det.h>

namespace LinBox {
/*  
Matrix is a polynomial matrix.  
result is set to its determinant and returned (a polynomial).
d is the number of evaluation points.

The method is to compute dets at each evaluation point and interpolate.
 (note by bds)
 */
template <class Field>
typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element&
computePolyDet(typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element& result,
				DenseMatrix<Givaro::Poly1Dom<Field,Givaro::Dense> >& A, 
               int d)
{
	typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
	typedef typename PolyDom::Element PolyElt;
	typedef typename Field::Element FieldElt;
	//typedef MatrixDomain<Field> FieldMatDom;
	typedef DenseMatrix<Field> FieldMat;

	int n=A.coldim(),m=A.rowdim();

	PolyDom BR=A.field();
	Field F(BR.subDomain()); // coeff field

	std::vector<FieldMat> mats;
	std::vector<FieldElt> pts(d);
	FieldElt fieldElt;
	F.assign(fieldElt,F.zero);
	for (int i=0;i<d;++i) {
		mats.push_back(FieldMat(F,m,n));
		F.init(pts[i],int64_t(i));
	}

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Initialized mats" << std::endl;


	PolyInterpolation<Field,PolyDom> PI(pts,F,BR);
#pragma omp parallel for shared(mats,A,PI,BR,F)
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			PolyElt p;
			A.getEntry(p,i,j);
			std::vector<FieldElt> vals;
			PI.evaluate(vals,p,BR,F);
			for (int k=0;k<d;++k) {
				mats[k].setEntry(i,j,vals[k]);
			}
		}
	}

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Finished evaluations" << std::endl;

	std::vector<FieldElt> dets(d);
#pragma omp parallel for shared(dets,mats)
	for (int k=0;k<d;++k) {
		det(dets[k],mats[k],Method::Elimination());
	}

	PI.interpolate(result,pts,dets,BR,F);
	return result;
}

int roundUpPowerOfTwo(unsigned int n)
{
	if (n==0) {
		return 0;
	} else if (n==1) {
		return 1;
	}
	
	int bits=0;
	int loopN=n;
	while (loopN>0) {
		loopN = loopN >> 1;
		++bits;
	}
	
	unsigned int mask=(1<<(bits-1))-1;
	if ((mask&n)!=0) {
		return 1<<bits;
	} else {
		return 1<<(bits-1);
	}
}

template <class Field,class Matrix>
typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element&
computePolyDetExtension(typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element& result,
                        Field& F,
                        Matrix& A)
{
	typedef Givaro::Poly1Dom<Field,Givaro::Dense> BasePolyDom;
	typedef typename BasePolyDom::Element BasePolyElt;

	BasePolyDom BR=F;

	int n=A.coldim(),m=A.rowdim();

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Computing d" << std::endl;

	int d=0;
	BasePolyElt p;
	for (int i=0;i<m;++i) {
		int rowMaxD=0;
		for (int j=0;j<n;++j) {
			A.getEntry(p,i,j);
			int t=BR.degree(p).value();
			rowMaxD=(rowMaxD<t)?t:rowMaxD;
		}
		d += rowMaxD;
	}
	d=roundUpPowerOfTwo(d+1);

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Found d" << std::endl;

	int a,e=1;
	a=F.cardinality();
	int newCard=a;
	while (newCard<d) {
		newCard *= a;
		++e;
	}
	
	typedef Givaro::Extension<Field> ExtField;
	typedef Givaro::Poly1Dom<ExtField,Givaro::Dense> ExtPolyDom;
	typedef typename ExtPolyDom::Element ExtPoly;
	//typedef typename Matrix::template rebind<ExtPolyDom>::other EPolyMatrix;
	typedef DenseMatrix<ExtPolyDom> EPolyMatrix;

	ExtField EF(F,e);
	ExtPolyDom EPD(EF,"x");
	Hom<Field,Givaro::Extension<Field> > hom(F,EF);

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Constructed new matrix" << std::endl;

	EPolyMatrix Ap(EPD,m,n);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			BasePolyElt oldElt;
			A.getEntry(oldElt,i,j);
			ExtPoly newElt;
			if (BR.isZero(oldElt)) {
				EPD.assign(newElt,EPD.zero);
			} else {
				int eltDegree=BR.degree(oldElt).value();
				EPD.init(newElt,BR.degree(oldElt));
				for (int k=0;k<eltDegree+1;++k) {
					hom.image(newElt[k],oldElt[k]);
				}
			}
			Ap.setEntry(i,j,newElt);
		}
	}

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Converted matrix" << std::endl;

	ExtPoly ep;
	computePolyDet<ExtField>(ep,Ap,d);
	//computePolyDet(ep,EF,Ap,d);

	commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
		<< "Computed over base field" << std::endl;

	if (EPD.isZero(ep)) {
		BR.assign(result,BR.zero);
	} else {
		int eltDegree=EPD.degree(ep).value();
		BR.init(result,EPD.degree(ep));
		for (int k=0;k<eltDegree+1;++k) {
			hom.preimage(result[k],ep[k]);
		}
	}
	return result;
}

}

#endif //__LINBOX_POLY_DET_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
