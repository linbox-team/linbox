/* Copyright (C) 2013 LinBox
 * Written by BB <bbboyer@ncsu.edu>
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

#include "linbox-config.h"

#include "benchmark.h"


// Least Squares
#ifdef __LINBOX_HAVE_LAPACK
extern "C" {
#if 1 // from lapack (not clapack)

	void dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
			double *b, int *ldb, double *work, int *lwork, int *info);

	void dgelsy_(int *m, int *n, int *nrhs, double *a, int *lda,
			double *b, int *ldb, int *JPVT, double *RCOND, int *RANK,
			double *work, int *lwork, int *info);

	void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,
			double *b, int *ldb, double *s, double *RCOND, int *RANK,
			double *work, int *lwork, int *info);
#endif

#if 0
	int clapack_dgels (const enum CBLAS_ORDER 	Order,
			   const enum CBLAS_TRANSPOSE 	TA,
			   int M,
			   int N,
			   int NRHS,
			   double * 	A,
			   int lda,
			   double * 	B,
			   const int 	ldb
			  );
#endif
}
#endif // __LINBOX_HAVE_LAPACK

// Curve fitting
namespace LinBox {

	// fit X[nn-1,nn],Y[nn-1,nn] and return evaluation at x.
	double fit2(const dvector_t & X, const dvector_t & Y, int nn, double x)
	{
		index_t n = (index_t) nn;
		assert(n>0);
		if ( n==1 ) {
			if ( X[0]==X[1] ) {
				// std::cerr << "two of your evaluation points are at the same X" << std::endl;
				// this is NOT supposed to happen.
				return (Y[0]+Y[1])/2 ;
			}
		}
		if (X[n]==X[n-1]) // discard the last one.
			return fit2(X,Y,(int)n-1,x);

		double a = (Y[n-1]-Y[n])/(X[n-1]-X[n]) ;
		double b = (X[n-1]*Y[n]-X[n]*Y[n-1])/(X[n-1]-X[n]) ;
		return a*x+b ;
	}

#ifdef __LINBOX_HAVE_LAPACK

	double fit_lapack3(const dvector_t &X, const dvector_t &Z, double x)
	{

		dvector_t Y = Z ;

		int n = (int) Z.size();
		linbox_check((size_t)n == X.size());

		int deg = (int) std::min((int)4,n);
		dvector_t V((size_t)(deg*n));

		int ldv = deg ;

#if 0 // Clapack (not working)
		for(index_t i = 0 ; i < (index_t)n; ++i) {
			for (index_t j = 0 ; j < (index_t)ldv; ++j) {
				V[i*ldv+j] = std::pow(X[i],j);
			}
		}

		clapack_dgels(CblasRowMajor, CblasNoTrans, n, deg, 1, &V[0],
			      deg, &Y[0], 1);

#endif

#if 1 /* basic least squares */
		// std::cout << V.size() << std::endl;
		for(index_t i = 0 ; i < (index_t)n; ++i) {
			for (index_t j = 0 ; j < (index_t)ldv; ++j) {
				V[i+j*(index_t)n] = std::pow(X[i],j);
			}
		}

		// std::cout << V << std::endl;

		int info;
		int ldun = 1 ;
		{

			int lwork = 2*n*deg*4 ;
			dvector_t work((size_t)lwork);
			char N[] = "N";
			dgels_(N, &n, &ldv, &ldun, &(V[0]) , &n, &(Y[0]), &n, &work[0], &lwork, &info);
		}

#endif

#if 0 /* least squares using SVN and V not nec. full rank */
		{
			int lwork = 2*deg+std::max(2*deg,n)*4 ;
			dvector_t work(lwork);
			dvector_t s(deg);
			double rcond = 1e-8 ;
			int rank ;
			dgelss_( &n, &ldv, &ldun, &(V[0]) , &n, &(Y[0]), &n, &s[0], &rcond, &rank, &work[0], &lwork, &info);
		}
#endif

#if 0 /* weighted least squares */
		//DGGGLM

		// TODO

#endif


		// std::cout << Y << std::endl;
		// horner eval the poly
		double res = 0.0;

		for(int i=deg-1; i >= 0; i--) {
			res = res * x + Y[(index_t)i];
		}
		return res;
	}
#endif // __LINBOX_HAVE_LAPACK


	double fit3(const dvector_t & X, const dvector_t & Y,int n, double x)
	{
#ifndef __LINBOX_HAVE_LAPACK /* à la main */
		linbox_check(n>1);
		linbox_check((size_t)n< X.size());
		linbox_check((size_t)n< Y.size());
		if (n==2) {
			if (X[1]==X[2])
				return fit2(X,Y,1,x) ;
			if (X[0]==X[2]) {
				return fit2(X,Y,1,x) ;
			}
			if (X[0]==X[1]) {
				dvector_t X1(2); X1[0]=X[1]; X1[2]=X[2];
				dvector_t Y1(2); Y1[0]=Y[1]; Y1[2]=Y[2];
				return fit2(X1,Y1,1,x) ;
			}
		}
		if (X[n]==X[n-1]) { // discard last
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));

			return fit3(X1,Y1,n-1,x) ;
		}
		if (X[n]==X[n-2]) { // discard last
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));

			return fit3(X1,Y1,n-1,x) ;
		}
		if (X[n-1]==X[n-2]) { // discard last but one
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));
			X1[n-1]=X[n];
			Y1[n-1]=Y[n];

			return fit3(X1,Y1,n-1,x) ;
		}

		// todo: use Lagrange ?
		// std::cout << X[n-2] << ',' << X[n-1] << ',' << X[n] << std::endl;
		double d  = (-X[n]+X[n-1])*(-X[n]+X[n-2])*(X[n-2]-X[n-1]) ;
		double a1 = -X[n]*Y[n-2]+X[n-2]*Y[n]+X[n-1]*Y[n-2]-X[n-1]*Y[n]+X[n]*Y[n-1]-X[n-2]*Y[n-1];
		double a2 = -X[n-2]*X[n-2]*Y[n]+X[n-2]*X[n-2]*Y[n-1]+X[n-1]*X[n-1]*Y[n]-Y[n-2]*X[n-1]*X[n-1]+Y[n-2]*X[n]*X[n]-Y[n-1]*X[n]*X[n];
		double a3 = X[n-2]*X[n-2]*X[n-1]*Y[n]-X[n-2]*X[n-2]*X[n]*Y[n-1]-X[n-1]*X[n-1]*X[n-2]*Y[n]+Y[n-1]*X[n-2]*X[n]*X[n]+X[n-1]*X[n-1]*X[n]*Y[n-2]-Y[n-2]*X[n-1]*X[n]*X[n];

		// std::cout <<" (("<<a1<<"*x+"<<a2<<")*x+"<<a3<<")/"<<d << std::endl;
		return ((a1*x+a2)*x+a3)/d ;
#else // __LINBOX_HAVE_LAPACK
		int m = min(n,5);
		dvector_t X1((index_t)m) ;
		dvector_t Y1((index_t)m) ;
		for (int i = 0 ; i < m ; ++i) X1[(index_t)i] = X[(index_t)(n-m+i)] ;
		for (int i = 0 ; i < m ; ++i) Y1[(index_t)i] = Y[(index_t)(n-m+i)] ;
		return fit_lapack3(X1,Y1,x);

#endif // __LINBOX_HAVE_LAPACK
	}

} // LinBox

// Terminal progression
namespace LinBox {

	void showAdvanceLinear(index_t curr, index_t min, index_t max)
	{
		std::cout << std::setprecision(4) << "\033[2K" << "\033[30D" << min <<std::flush;
		std::cout << '<' << curr << '<' << max << " (" << std::flush;
		std::cout << double(curr-min)/double(max-min)*100 << "%)" << std::flush;
	}

	void showFinish(index_t curr, index_t all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "finished : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}

	void showSkip(index_t curr, index_t all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "skipped : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}

} // LinBox

// MFLOPS helper
namespace LinBox {

	double computeMFLOPS(const double & tim, const double mflo, const index_t rpt)
	{
		linbox_check(rpt);
		// linbox_check(tim != 0.);
		if (tim == 0.) return NAN ;
		return (double) ((mflo*rpt)/tim);
	}

	dvector_t &
	insertTime(dvector_t & tim3, const double & tps)
	{
		linbox_check(tim3.size() == 3);
		if (tim3[0] > tps) {
			tim3[2] = tim3[1] ;
			tim3[1] = tim3[0] ;
			tim3[0] = tps ;
		}
		else if (tim3[1] > tps) {
			tim3[2] = tim3[1] ;
			tim3[1] = tps ;
		}
		else if (tim3[2] > tps) {
			tim3[2] = tps ;
		}
		return tim3 ;
	}

	double computeMFLOPS(const dvector_t & tim, const double mflo, LINBOX_enum(Tag::TimeSelect) ts )
	{
		linbox_check(tim.size());
		switch (ts) {
		case (Tag::TimeSelect::average) :
			{
				double tps = 0 ;
				for (size_t i = 0 ; i < tim.size() ; ++i)
					tps += tim[i] ;
				return computeMFLOPS(tps,mflo,(index_t)tim.size());
			}
		case (Tag::TimeSelect::bestThree) :
			{
				if (tim.size() <4)
					return computeMFLOPS(tim,mflo,Tag::TimeSelect::average);

				dvector_t tps (3);
				double t1,t2 ;
				if (tim[0]<tim[1]) {
					t1 = tim[0];
					t2 = tim[1];
				}
				else {
					t1 = tim[1];
					t2 = tim[0];
				}
				if (tim[3] < t1) {
					tps[0] = tim[3] ;
					tps[1] = t1 ;
					tps[2] = t2 ;
				}
				else if (tim[2] < t1) {
					tps[0] = t1 ;
					tps[1] = tim[3] ;
					tps[2] = t2 ;
				}
				else {
					tps[0] = t1 ;
					tps[1] = t2;
					tps[2] = tim[3] ;
				}

				for (size_t i = 3 ; i < tim.size() ; ++i)
					insertTime(tps,tim[i]);

				return computeMFLOPS(tim,mflo,Tag::TimeSelect::average);

			}
		case (Tag::TimeSelect::bestOne) :
			{
				double t1 = tim[0] ;
				for (size_t i = 1 ; i < tim.size() ; ++i)
					if (tim[i] < t1)
						t1 = tim[i] ;
				return computeMFLOPS(t1,mflo,1);

			}
		case (Tag::TimeSelect::median) :
			{
				if (tim.size() == 1)
					return computeMFLOPS(tim[0],mflo,1) ;

				dvector_t tps (tim);
				std::sort(tps.begin(),tps.end());
				index_t mid = (index_t)tps.size()/2 ;
				double t1 ;
				if (Givaro::isOdd(tps.size()))
					t1 = tps[mid] ;
				else
					t1 = (tps[mid-1]+tps[mid])/2;
				return computeMFLOPS(t1,mflo,1);
			}
		case (Tag::TimeSelect::medmean) :
			{
				if (tim.size() < 3)
					return computeMFLOPS(tim,mflo,Tag::TimeSelect::median); ;

				index_t q1 = (index_t)std::round((double)tim.size()/(double)4) ;
				index_t q3 = (index_t)tim.size()-q1 ;
				dvector_t tps (tim);
				std::sort(tps.begin(),tps.end());
				dvector_t tps2 (tim.begin()+q1,tim.begin()+q3);
				return computeMFLOPS(tps2,mflo,Tag::TimeSelect::average);
			}

		default :
			{
				throw("not among tags");
			}

		} // switch(ts)
	}

}


	//
	//  TimeWatcher
	//

namespace LinBox {
	TimeWatcher::TimeWatcher (dvector_t & pts, dvector_t & vals) :
		Points_(&pts)
		,Values_(&vals)
		,MaxRepet_(12)
		,MinRepet_(2)
		,MaxTime_(0.5)
		,AbortTime_(10)
		,aborted_(false)
	{
		linbox_check(vals.size() == pts.size());
	}

	TimeWatcher::TimeWatcher () :
		Points_(NULL)
		,Values_(NULL)
		,MaxRepet_(12)
		,MinRepet_(2)
		,MaxTime_(0.5)
		,AbortTime_(10)
		,aborted_(false)
	{
	}


	void TimeWatcher::init(dvector_t & pts, dvector_t & vals)
	{
		linbox_check(vals.size() == pts.size());
		Points_ = &pts ;
		Values_ = &vals ;
	}


	dvector_t & TimeWatcher::refX()
	{
		return *Points_ ;
	}

	dvector_t & TimeWatcher::refY()
	{
		return *Values_ ;
	}

	// we don't assume that t(0sec) = 0 unless nothing has been computed yet.
	double TimeWatcher::predict(double x)
	{
		index_t Current_ = size();
		if (Current_ == 0)
			return 0 ;
		if (Current_ ==1 )
			return refX()[Current_-1] ; // unknown. could be known if we suppose t(0)=0
		if (Current_ == 2 ) {
			return fit2(refX(),refY(),1,x);
		}
		return fit3(refX(),refY(), (int) Current_-1,x);
	}

	bool TimeWatcher::keepon( index_t & repet, double tim, bool usePrediction )
	{

		if (aborted_)
			return false ;

		if (usePrediction) {
			if (predict(tim) < AbortTime_)
				return true ;
			else{
				aborted_ = true ;
				return false;
			}
		}

		if (repet<MinRepet_ || (tim < MaxTime_ && repet < MaxRepet_) ) {
			++repet ;
			return true;
		}
		return false ;
	}

} // LinBox

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
