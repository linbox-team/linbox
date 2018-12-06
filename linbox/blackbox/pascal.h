#ifndef __LINBOX_pascal_H
#define __LINBOX_pascal_H

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#ifdef __LINBOX_USE_OPENMP
#include <omp.h>
#endif

#include "linbox/util/debug.h"
#include "linbox/matrix/sliced3.h"
#include "linbox/blackbox/blockbb.h"

#define PASCAL_BASECASE_THRESH 81

namespace LinBox {

template <class Field>
class ModularNChooseK {
public:
	typedef typename Field::Element Element;

	int reduceZeros(int x)
	{
		while ((x>0)&&(x%q_==0)) {
			x/=q_;
		}
		return x;
	}

	void initFactorials()
	{
		factList_.resize(maxN_+1);
		Element d;
		F_.assign(factList_[0],F_.one);
		for (int i=1;i<maxN_+1;++i) {
			F_.init(d,reduceZeros(i));
			if (d==0) {
				F_.assign(factList_[i],factList_[i-1]);
			} else {
				F_.mul(factList_[i],factList_[i-1],d);
			}
		}
	}

	void initPowerList()
	{
		powerList_.resize(maxN_+1);
		for (int i=0;i<maxN_+1;++i) {
			powerList_[i]=0;
		}
		for (int k=q_;k<maxN_;k*=q_) {
			for (int i=k-1;i<maxN_;i+=k) {
				++powerList_[i+1];
			}
		}
		
		int sum=0;
		for (int i=0;i<maxN_+1;++i) {
			sum+=powerList_[i];
			powerList_[i]=sum;
		}
	}

	ModularNChooseK(Field& F, int q, int maxN) :
		F_(F), q_(q), maxN_(maxN)
	{
		
		initPowerList();
		initFactorials();
	}

	Element& compute(Element& d,int n,int k)
	{
		int denomPower=powerList_[k]+powerList_[n-k];
		int numerPower=powerList_[n];
		if (numerPower>denomPower) {
			F_.assign(d,F_.zero);
		} else {
			linbox_check(numerPower==denomPower);
			F_.mul(d,factList_[n-k],factList_[k]);
			F_.div(d,factList_[n],d);
		}
		return d;
		
	}
protected:
	Field F_;
	int q_;
	std::vector<int> powerList_;
	std::vector<Element> factList_;
	int maxN_;
};

template <class Field>
class PascalBlackbox {
public:
	typedef typename Field::Element Element;

	template <class Vector>
	PascalBlackbox(int r,
	               int c,
	               const Vector& polyCoeffs,
	               Field& F) :
		polyCoeffs_(polyCoeffs),
		rowdim_(r),coldim_(c),
		F_(F)
	{
		initPowersOfThree();
		initBaseCase();
	}

	void initBaseCase()
	{
		ModularNChooseK<Field> Choose(F_,3,std::max(rowdim_+coldim_,2*PASCAL_BASECASE_THRESH)+1);
		nonZeroRows_.clear();
		nonZeroCols_.clear();
		nonZerosUnflipped_.clear();
		nonZerosFlipped_.clear();
		int r=PASCAL_BASECASE_THRESH;
		Element two;
		F_.add(two,F_.one,F_.one);
		for (int k=0;k<(2*r)-1;++k) {
			for (int j=(k<r)?0:(1+k-4);j<std::min(k+1,r);++j) {
				Element d;
				int i=k-j;
				Choose.compute(d,k,j);
				if (!F_.isZero(d)) {
					nonZeroRows_.push_back(i);
					nonZeroCols_.push_back(j);
					nonZerosUnflipped_.push_back(d);
					F_.mulin(d,two);
					nonZerosFlipped_.push_back(d);
				}
			}
		}
	}

	template <class Mat1, class Mat2>
	Mat1& applyLeft(Mat1& lhs,const Mat2& rhs) const
	{
		int r=std::max(PASCAL_BASECASE_THRESH,
		               nextPowerOfThree(std::max(rowdim_,coldim_)));
		lhs.zero();
		applyLeft(0,0,r,lhs,const_cast<Mat2&>(rhs),false,1);
		return lhs;
	}

	template<class Mat1, class Mat2>
	void applyLeft(int i0, int j0, int r, Mat1& lhs,Mat2& rhs, bool flip, int numThreads) const
	{
		if (i0 > rowdim_ || j0 > coldim_) return;
		if (r<=PASCAL_BASECASE_THRESH) {
			if (flip) {
				applyHelper(i0,j0,r,lhs,rhs,nonZerosFlipped_);
			} else {
				applyHelper(i0,j0,r,lhs,rhs,nonZerosUnflipped_);
			}
		} else {
			int r3=r/3;

//CP: removing all explicit OMP calls
//    to be replaced by calls to FFLAS-FFPACK's DSL

// #ifdef __LINBOX_USE_OPENMP
// 			int threadLimit=omp_get_max_threads();
// 			//int threadLimit=4;
// #endif

// #pragma omp parallel sections if (numThreads+2<=threadLimit) shared(lhs,rhs)
// 			{
// #pragma omp section
				applyLeft(i0,j0,r3,lhs,rhs,flip,numThreads+2);
// #pragma omp section
				applyLeft(i0+r3,j0,r3,lhs,rhs,flip,numThreads+2);
// #pragma omp section
				applyLeft(i0+2*r3,j0,r3,lhs,rhs,flip,numThreads+2);
// 			}

// #pragma omp parallel sections if (numThreads+1<=threadLimit) shared(lhs,rhs)
// 			{
// #pragma omp section
				applyLeft(i0,j0+r3,r3,lhs,rhs,flip,numThreads+1);
// #pragma omp section
				applyLeft(i0+r3,j0+r3,r3,lhs,rhs,!flip,numThreads);
			// }


			applyLeft(i0,j0+2*r3,r3,lhs,rhs,flip,numThreads+1);

		}
	}

	template<class SlicedDom>
	void applyHelper(int i0,int j0,int r,Sliced<SlicedDom>& lhs,Sliced<SlicedDom>& rhs,
	                 const std::vector<uint16_t>& nonZeros) const
	{
		int nnz=nonZeroRows_.size();
		for (int l=0;l<nnz;++l) {
			int i=nonZeroRows_[l]+i0;
			int j=nonZeroCols_[l]+j0;
			int k=i+j;
			if (i>=rowdim_ || j>=coldim_) {
				continue;
			}
			if (F_.isZero(polyCoeffs_[k])) {
				continue;
			}
			Element d;
			d=nonZeros[l];
			F_.mul(d,d,polyCoeffs_[k]);
			typename Sliced<SlicedDom>::RawIterator Ab(lhs.rowBegin(i)),
				Ae(lhs.rowEnd(i)), Bb(rhs.rowBegin(j));
			lhs.axpyin(Ab,Ae,d,Bb);
		}
	}

	template<class Mat1, class Mat2>
	void applyHelper(int i0,int j0,int r,Mat1& lhs,Mat2& rhs,const std::vector<uint16_t>& nonZeros) const
	{
		int nnz=nonZeroRows_.size();
		int w=rhs.coldim();
		for (int l=0;l<nnz;++l) {
			int i=nonZeroRows_[l]+i0;
			int j=nonZeroCols_[l]+j0;
			int k=i+j;
			if (i>=rowdim_ || j>=coldim_) {
				continue;
			}
			if (F_.isZero(polyCoeffs_[k])) {
				continue;
			}
			Element d;
			d=nonZeros[l];
			F_.mul(d,d,polyCoeffs_[k]);
			typename Mat2::constSubMatrixType Xr(rhs,j,0,1,w);
			typename Mat1::subMatrixType Yr(lhs,i,0,1,w);
			lhs._MD.saxpyin(Yr,d,Xr);
		}
	}

	size_t rowdim() const {return rowdim_;}
	size_t coldim() const {return coldim_;}

	const Field& field() const {
		return F_;
	}
	
protected:
	std::vector<Element> polyCoeffs_;
	int rowdim_,coldim_;
	Field F_;
	std::vector<int> powersOfThree_;
	std::vector<uint16_t> nonZeroRows_,nonZeroCols_,nonZerosFlipped_,nonZerosUnflipped_;

	void initPowersOfThree()
	{
		int x=1;
		while (x<std::numeric_limits<int>::max()/3) {
			powersOfThree_.push_back(x);
			x*=3;
		}
		powersOfThree_.push_back(x);
	}
	int nextPowerOfThree(int x) const
	{
		int i=0;
		while (powersOfThree_[i] < x) {++i;}
		return powersOfThree_[i];
	}
};

template<class Field>
struct is_blockbb<PascalBlackbox<Field>> {
	static const bool value = true;
};

}

#endif //__LINBOX_pascal_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
