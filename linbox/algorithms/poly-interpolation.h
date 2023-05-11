#ifndef __LINBOX_POLY_INTERPOLATION_H
#define __LINBOX_POLY_INTERPOLATION_H


#ifdef __LINBOX_USE_OPENMP
#include <omp.h>
#endif

namespace LinBox {


	// Expects pts to have size 2^n, undefined behavior otherwise
template<class Field, class PolyDom>
class PolyInterpolation {
public:
	typedef GivaroPoly<PolyDom> Ring;
	typedef typename Ring::Element RingElt;
	typedef typename Field::Element FieldElt;
	typedef std::vector<std::vector<RingElt> > ProductTree;

	PolyInterpolation(const std::vector<FieldElt>& pts,
	                  Field& F,
	                  PolyDom& PD)
	{
		int n=pts.size();

		std::vector<RingElt> mRow;
		for (int i=0;i<n;++i) {
			RingElt p;
			PD.assign(p,PD.one);
			PD.shiftin(p,1);
			PD.subin(p,pts[i]);
			mRow.push_back(p);
		}
		Mtree_.push_back(mRow);
		k_=0;
		while (Mtree_[k_].size()>1) {
			std::vector<RingElt> row;
			for (uint32_t i = 0; i < Mtree_[k_].size()/2; ++i) {
				RingElt p;
				PD.mul(p,Mtree_[k_][2*i],Mtree_[k_][2*i+1]);
				row.push_back(p);
			}
			Mtree_.push_back(row);
			++k_;
		}
	}

	RingElt& interpolate(RingElt& poly,
	                     const std::vector<FieldElt>& pts,
	                     const std::vector<FieldElt>& vals,
	                     PolyDom& PD,
	                     Field& F)
	{
		int n=pts.size();
		std::vector<FieldElt> si;
		RingElt mprime;
		PD.diff(mprime,Mtree_[k_][0]);
		evaluate(si,mprime,PD,F);
		for (int i=0;i<n;++i) {
			F.invin(si[i]);
			F.mulin(si[i],vals[i]);
		}
		scaledSum(poly,si,PD);
		return poly;
	}

	RingElt& scaledSum(RingElt& comb,
	                   const std::vector<FieldElt>& cs,
	                   PolyDom& PD)
	{
		int numPts=cs.size();
		std::vector<RingElt> fRow,tempRow;
		fRow.resize(numPts);
		for (int i=0;i<numPts;++i) {
			PD.assign(fRow[i],cs[i]);
		}

		for (int i=0;i<k_;++i) {
			int rowLen=fRow.size();
			tempRow.resize(rowLen/2);
			for (int j=0;j<rowLen/2;++j) {
				RingElt p,q;
				PD.mul(p,fRow[2*j],Mtree_[i][2*j+1]);
				PD.mul(q,fRow[2*j+1],Mtree_[i][2*j]);
				PD.addin(p,q);
				PD.assign(tempRow[j],p);
			}
			fRow.swap(tempRow);
			tempRow.clear();
		}
		PD.assign(comb,fRow[0]);
		return comb;
	}

	void evaluate(std::vector<FieldElt>& vals,
	              const RingElt& poly,
	              PolyDom& PD,
	              Field& F)
	{
		std::vector<RingElt> fRow,tempRow;
		fRow.push_back(poly);

		for (int i=k_-1;i>=0;--i) {
			int rowLen=fRow.size();
			tempRow.resize(rowLen*2);
			for (int j=0;j<rowLen;++j) {
				RingElt p;
				PD.mod(p,fRow[j],Mtree_[i][2*j]);
				PD.assign(tempRow[2*j],p);
				PD.mod(p,fRow[j],Mtree_[i][2*j+1]);
				PD.assign(tempRow[2*j+1],p);
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


	static void naiveInterpolate(RingElt& poly,
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
protected:

	ProductTree Mtree_;

	int k_;

};


}

#endif //__LINBOX_POLY_INTERPOLATION_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
