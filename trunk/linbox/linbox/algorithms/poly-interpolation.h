
#ifndef __LINBOX_POLY_INTERPOLATION_H
#define __LINBOX_POLY_INTERPOLATION_H

#include <linbox/ring/givaro-poly.h>

namespace LinBox {

template<class Field, class PolyDom>
class PolyInterpolation {
public:
	typedef GivaroPoly<PolyDom> Ring;
	typedef typename Ring::Element RingElt;
	typedef typename Field::Element FieldElt;
	typedef std::vector<std::vector<RingElt> > ProductTree;



	RingElt& interpolate(RingElt& poly,
	                     const std::vector<FieldElt>& pts,
	                     const std::vector<FieldElt>& vals,
	                     PolyDom& PD,
	                     Field& F)
	{
		int n=pts.size();
		ProductTree Mtree;
		int k=productTree(Mtree,pts,PD);
		std::vector<FieldElt> si;
		RingElt mprime;
		PD.diff(mprime,Mtree[k][0]);
		evaluate(si,mprime,Mtree,PD,F);
		for (int i=0;i<n;++i) {
			F.invin(si[i]);
			F.mulin(si[i],vals[i]);
		}
		scaledSum(poly,si,Mtree,PD);
		return poly;
	}

	// Expects cs.size()==2^k, will segfault otherwise
	// FIXME: Fail gracefully
	RingElt& scaledSum(RingElt& comb,
	                   const std::vector<FieldElt>& cs,
	                   const ProductTree& Mtree,
	                   PolyDom& PD)
	{
		int k=Mtree.size()-1,numPts=cs.size();
		std::vector<RingElt> fRow,tempRow;
		fRow.resize(numPts);
		for (int i=0;i<numPts;++i) {
			PD.assign(fRow[i],cs[i]);
		}

		for (int i=0;i<k;++i) {
			int rowLen=fRow.size();
			for (int j=0;j<rowLen/2;++j) {
				RingElt p,q;
				PD.mul(p,fRow[2*j],Mtree[i][2*j+1]);
				PD.mul(q,fRow[2*j+1],Mtree[i][2*j]);
				PD.addin(p,q);
				tempRow.push_back(p);
			}
			fRow.swap(tempRow);
			tempRow.clear();
		}
		PD.assign(comb,fRow[0]);
		return comb;
	}

	int productTree(ProductTree& Mtree,
	                const std::vector<FieldElt>& pts,
	                PolyDom& PD)
	{
		int n=pts.size(),k;

		std::vector<RingElt> mRow;
		for (int i=0;i<n;++i) {
			RingElt p;
			PD.assign(p,PD.one);
			PD.shiftin(p,1);
			PD.subin(p,pts[i]);
			mRow.push_back(p);
		}
		Mtree.push_back(mRow);
		k=0;
		while (Mtree[k].size()>1) {
			std::vector<RingElt> row;
			for (int i=0;i<Mtree[k].size()/2;++i) {
				RingElt p;
				PD.mul(p,Mtree[k][2*i],Mtree[k][2*i+1]);
				row.push_back(p);
			}
			Mtree.push_back(row);
			++k;
		}
		return k;
	}

	void evaluate(std::vector<FieldElt>& vals,
	              const std::vector<FieldElt>& pts,
	              const RingElt& poly,
	              PolyDom& PD,
	              Field& F)
	{
		ProductTree Mtree;
		productTree(Mtree,pts,PD);
		evaluate(vals,poly,Mtree,PD,F);
	}

	void evaluate(std::vector<FieldElt>& vals,
	              const RingElt& poly,
	              ProductTree& Mtree,
	              PolyDom& PD,
	              Field& F)
	{
		int k=Mtree.size()-1;
		std::vector<RingElt> fRow,tempRow;
		fRow.push_back(poly);

		for (int i=k-1;i>=0;--i) {
			int rowLen=fRow.size();
			for (int j=0;j<rowLen;++j) {
				if (PD.degree(fRow[j])==0) {
					tempRow.push_back(fRow[j]);
					tempRow.push_back(fRow[j]);
				} else {
					RingElt p;
					PD.mod(p,fRow[j],Mtree[i][2*j]);
					tempRow.push_back(p);
					PD.mod(p,fRow[j],Mtree[i][2*j+1]);
					tempRow.push_back(p);
				}
			}
			fRow.swap(tempRow);
			tempRow.clear();
		}

		int numPts=fRow.size();
		vals.resize(numPts);
		for (int i=0;i<numPts;++i) {
			FieldElt d;
			PD.leadcoef(d,fRow[i]);
			vals[i]=d;
		}
	}

	void naiveInterpolate(RingElt& poly,
	                      const std::vector<FieldElt>& vals,
	                      const std::vector<FieldElt>& pts,
	                      GivaroPoly<PolyDom>& R,
	                      PolyDom& PD,
	                      Field& F)
	{
		typedef GivaroPoly<PolyDom> Ring;
		typedef typename Ring::Element RingElt;
		typedef typename Ring::Scalar_t FieldElt;

		int n=vals.size();
		linbox_check(vals.size()==pts.size());
		PD.assign(poly,PD.zero);
		for (int i=0;i<n;++i) {
			RingElt mulElt;
			PD.assign(mulElt,PD.one);
			for (int j=0;j<n;++j) {
				if (i==j) {
					continue;
				}
				RingElt p;
				PD.assign(p,PD.one);
				PD.shiftin(p,1);
				PD.subin(p,pts[j]);
				FieldElt d;
				F.sub(d,pts[i],pts[j]);
				PD.divin(p,d);
				PD.mulin(mulElt,p);
			}
			PD.mulin(mulElt,vals[i]);
			PD.addin(poly,mulElt);
		}
	}
};

}

#endif //__LINBOX_POLY_INTERPOLATION_H
