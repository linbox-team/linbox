
#ifndef __LINBOX_POLY_DET_H
#define __LINBOX_POLY_DET_H

#include "linbox/field/Givaro/givaro-extension.h"
#include "linbox/solutions/det.h"

namespace LinBox {

template <class Field,class Matrix>
typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element&
computePolyDet(typename Givaro::Poly1Dom<Field,Givaro::Dense>::Element& result,
               Field& F,
               Matrix& A, 
               int d)
{
	typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
	typedef typename PolyDom::Element PolyElt;
	typedef typename Field::Element FieldElt;
	typedef MatrixDomain<Field> FieldMatDom;
	typedef typename FieldMatDom::OwnMatrix FieldMat;

	int n=A.coldim(),m=A.rowdim();

	PolyDom BR=A.field();

	std::vector<FieldMat> mats;
	std::vector<FieldElt> pts(d);
	FieldElt fieldElt;
	F.assign(fieldElt,F.zero);
	for (int i=0;i<d;++i) {
		mats.push_back(FieldMat(F,m,n));
		F.assign(pts[i],fieldElt);
		F.next(fieldElt);
	}

	std::vector<FieldElt> vals;
	PolyInterpolation<Field,PolyDom> PI;
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			PolyElt p;
			A.getEntry(p,i,j);
			PI.evaluate(vals,pts,p,BR,F);
			for (int k=0;k<d;++k) {
				mats[k].setEntry(i,j,vals[k]);
			}
		}
	}
	std::vector<FieldElt> dets(d);
	for (int k=0;k<d;++k) {
		det(dets[k],mats[k],Method::Wiedemann());
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

	BasePolyDom BR=A.field();

	int n=A.coldim(),m=A.rowdim();

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

	int a,e=1;
	a=F.cardinality();
	int newCard=a;
	while (newCard<d) {
		newCard *= a;
		++e;
	}
	
	typedef GivaroExtension<Field> ExtField;
	typedef Givaro::Poly1Dom<ExtField,Givaro::Dense> ExtPolyDom;
	typedef typename ExtPolyDom::Element ExtPoly;
	typedef typename Matrix::template rebind<ExtPolyDom>::other EPolyMatrix;

	ExtField EF(F,e);
	ExtPolyDom EPD(EF,"x");
	Hom<Field,GivaroExtension<Field> > hom(F,EF);

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

	ExtPoly ep;
	computePolyDet(ep,EF,Ap,d);

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
