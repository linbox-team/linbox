/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
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

/*! @file benchmarks/benchmark-fgemm.C
 * @ingroup benchmarks
 * @brief Benchmarking dense matrix multiplication on finite fields.
 * This file benchmarks the FFLAS::fgemm implementation for various fields,
 * shape and parameters. Actually, we use the wrapper member \c mul of LinBox::BlasMatrixDomain.
 * @todo make graphs look better (legends, units,...)
 */

#include "benchmarks/benchmark.h"
#include "linbox/util/error.h"
#include "fflas-ffpack/fflas-ffpack.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"


using LinBox::TimeWatcher;

#if 0
bool keepon(index_t & repet, const double & tim, double maxtime=0.2)
{
	if (repet<2 || tim < maxtime) {
		++repet ;
		return true;
	}
	return false ;
}
#endif

double fgemm_mflops(int m, int n, int k)
{
	return 2*(double)m/100*(double)n/100*(double)k/100 ;
}

double compute_mflops(const Timer & t, const double mflo, const int rpt = 1)
{
	linbox_check(rpt);
	return (double) ((mflo*rpt)/t.usertime());
}

double compute_mflops(const double & t, const double mflo, const int rpt = 1)
{
	linbox_check(rpt);
	return (double) ((mflo*rpt)/t);
}

/*! @internal
 * @brief launches the benchmarks for the square case.
 * @param F field
 * @param min min size to bench
 * @param max max size to bench
 * @param step step between two sizes
 * @param Data where data is stored
 * @param series_nb index of the current series.
 */
template<class Field>
void launch_bench_square(Field & F // const problem
			 , index_t min, index_t max, int step
			 , LinBox::PlotData<index_t> & Data
			 , index_t series_nb)
{
	TimeWatcher TW(10,series_nb);
	index_t l = 0 ;
	Timer fgemm_sq_tim ;
	Timer chrono ;
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	for ( index_t i = min ; i < max ; i += step , ++l ) {
		showAdvanceLinear(i,min,max);
		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		LinBox::BlasMatrix<Field> A (F,ii,ii);
		LinBox::BlasMatrix<Field> B (F,ii,ii);
		LinBox::BlasMatrix<Field> C (F,ii,ii);
		if (!series_nb)
			Data.setAbsciName(l,i); // only write abscissa for serie 0
		index_t j = 0 ; // number of repets.

		RandMat.random(A);
		RandMat.random(B);
		RandMat.random(C);
		fgemm_sq_tim.clear() ;
		while( TW.keepon(j, fgemm_sq_tim) ) {
			chrono.clear() ; chrono.start() ;
			BMD.mul(C,A,B) ; // C = AB
			chrono.stop();
			fgemm_sq_tim += chrono ;
		}
		if (!j){
			std::cout << "multiplication did not happen" << std::endl;
		}
#ifdef _LB_DEBUG
		else {
			std::cout << i << ',' << j << std::endl;
		}
#endif
		mflops = compute_mflops(fgemm_sq_tim,fgemm_mflops(i,i,i),j);
		Data.setEntry(series_nb,l,mflops);

	}
	std::ostringstream nam ;
	nam << '\"' ;
	F.write(nam);
	nam << '\"' ;
	Data.setSerieName(series_nb,nam.str());

}

/*! @internal
 * @brief launches the benchmarks for the square case directly BLAS.
 * @param min min size to bench
 * @param max max size to bench
 * @param step step between two sizes
 * @param Data where data is stored
 * @param series_nb index of the current series.
 */
template<class Field>
void launch_bench_blas(Field & F
		       , index_t min, index_t max, int step
		       , LinBox::PlotData<index_t> & Data
		       , index_t series_nb
		       )
{
	TimeWatcher TW(10,series_nb);
	// typedef LinBox::Modular<T> Field ;
	// Field F((int)charact);
	index_t l = 0 ;
	Timer fgemm_blas_tim ;
	Timer chrono ;
	double mflops ;
	typedef typename Field::Element  Element ;
	typedef typename Field::RandIter Randiter;
	Randiter R(F) ;
	// LinBox::BlasMatrixDomain<Field> BMD(F) ;
	// LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	index_t mm = max * max ;
	Element *  A = new Element[mm] ;
	Element *  B = new Element[mm] ;
	Element *  C = new Element[mm] ;
	// typedef typename LinBox::UnparametricField<T> FloatingDomain ;
	// FloatingDomain G ;

	for ( index_t i = min ; i < max ; i += step , ++l ) {
		showAdvanceLinear(i,min,max);

		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		index_t mimi = (index_t) ii*ii ;
		if (!series_nb)
			Data.setAbsciName(l,i); // only write abscissa for serie 0

		for (index_t j = 0 ; j < mimi ; ++j) R.random(A[j]);
		for (index_t j = 0 ; j < mimi ; ++j) R.random(B[j]);
		for (index_t j = 0 ; j < mimi ; ++j) R.random(C[j]);

		index_t j = 0 ;
		fgemm_blas_tim.clear() ;
		// double fgemm_blas_tim = 0 ;
		while(TW.keepon(j,fgemm_blas_tim)) {
			chrono.clear(); chrono.start() ;
			FFLAS::fgemm((typename Field::Father_t)F,FFLAS::FflasNoTrans,FFLAS::FflasNoTrans,
					     ii,ii,ii,
					     F.one,
					     A,ii,B,ii,
					     F.zero,
					     C,ii) ;
			chrono.stop() ;
			fgemm_blas_tim += chrono ;
		}
		mflops = compute_mflops(fgemm_blas_tim,fgemm_mflops(i,i,i),j);
// #ifndef NDEBUG
		if (i >=950 && i <= 1050 ) {
			std::cerr << std::endl<< typeid(Field).name() << ' ' << i << ':' << mflops << std::endl;
			std::cerr << fgemm_blas_tim << std::endl;
		}
// #endif

		Data.setEntry(series_nb,l,mflops);
	}

	delete[] A ;
	delete[] B ;
	delete[] C ;
	std::ostringstream nam ;
	nam << '\"' ;
	F.write(nam);
	nam << '\"' ;
	Data.setSerieName(series_nb,nam.str());

}

/*! @internal
 * @brief launches the benchmarks for the square case.
 * C= AB
 * @param F field
 * @param m rows in A
 * @param k cols in A
 * @param n cols in C
 * @param Data where data is stored
 * @param point_nb point to be computed
 */
template<class Field>
void launch_bench_rectangular(Field & F // const problem
			      , int m, int k, int n
			      , LinBox::PlotData<std::string> & Data
			      , index_t point_nb)
{
	TimeWatcher TW(10,0);
	Timer fgemm_rect_tim ;
	Timer chrono ; chrono.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	LinBox::BlasMatrix<Field> A (F,m,k);
	LinBox::BlasMatrix<Field> B (F,k,n);
	LinBox::BlasMatrix<Field> C (F,m,n);
	index_t j = 0 ;
	fgemm_rect_tim.clear() ;
	while (TW.keepon(j,fgemm_rect_tim)) {
		RandMat.random(A);
		RandMat.random(B);
		RandMat.random(C);
		chrono.clear() ; chrono.start() ;
		BMD.mul(C,A,B) ; // C = AB
		chrono.stop();
		fgemm_rect_tim += chrono ;
	}
	if (!j) {
		std::cout << "multiplication did not happen" << std::endl;
	}
#ifdef _LB_DEBUG
	else {
		std::cout << point_nb << std::endl;
	}
#endif
	mflops = compute_mflops(fgemm_rect_tim,fgemm_mflops(m,k,n),j);
	Data.setEntry(0,point_nb,mflops);
	std::ostringstream nam ;
	nam << "\"(" << m << ',' << k << ',' << n << ")\"" ;
	Data.setAbsciName(point_nb,nam.str());
}

/*! @internal
 * @brief launches the benchmarks for various parameters of a, b.
 * D = aAB+bC and C = aAB+bC ("in place" versions)
 * Are tested the following couples \c (a,b) (where \c p is invertible in \p F. This has
 * to be true when \c a=p.) :
 * - b=0 and a=1,-1,p ;
 * - a=1 and b=1,-1,p ;
 * - a=-1 and id. ;
 * - a=p and id. ;
 * .
 * We call xA = ^tA if tA is true, A otherwise.
 * @param F field
 * @param m rows in xA
 * @param k cols in xA and rows in xB
 * @param n cols in xB
 * @param alpha alpha (not zero)
 * @param beta beta
 * @tparam tA is A transposed ?
 * @tparam tB is B transposed ?
 * @param Data where data is stored
 * @param point_nb point to be computed
 * @param inplace in place or not (ie C is overwritten or not ? default = \c false.
 */
template<class Field, bool tA, bool tB>
void launch_bench_scalar(Field & F // const problem
			 , int m, int k, int n
			 , const typename Field::Element & alpha, const typename Field::Element & beta
			 , LinBox::PlotData<std::string> & Data
			 , index_t point_nb
			 , bool inplace = false)
{
	TimeWatcher TW(10,0);
	Timer fgemm_scal_tim ;
	Timer chrono ;
	fgemm_scal_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	typedef typename LinBox::BlasMatrix<Field >  Matrix ;
	typedef typename LinBox::TransposedBlasMatrix<Matrix > TransposedMatrix ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	int mm = tA?k:m ;
	int kk = tA?m:k ;
	int nn = tB?k:n ;
	Matrix A (F,mm,kk);
	Matrix B (F,kk,nn);
	Matrix C (F,m,n);
	Matrix D (F,m,n);
	TransposedMatrix At(A);
	TransposedMatrix Bt(B);
	// LinBox::BlasMatrix<Field > A (mm,kk);
	// LinBox::BlasMatrix<Field > B (kk,nn);
	// LinBox::BlasMatrix<Field > C (m,n);
	// LinBox::BlasMatrix<Field > D (m,n);


	// LinBox::TransposedBlasMatrix<LinBox::BlasMatrix<Field > > At(A);
	// LinBox::TransposedBlasMatrix<LinBox::BlasMatrix<Field > > Bt(B);

	index_t j = 0 ;
	while (TW.keepon(j,fgemm_scal_tim)) {
		RandMat.random(A);
		RandMat.random(B);
		RandMat.random(C);

		chrono.clear() ; chrono.start();

		if (inplace) {
			if (tA) {
				if (tB) {
					BMD.muladdin(beta,C,alpha,At,Bt) ; // C = alphaAB+beta C
				}
				else{
					BMD.muladdin(beta,C,alpha,At,B) ;
				}
			}
			else {
				if (tB) {
					BMD.muladdin(beta,C,alpha,A,Bt) ;
				}
				else{
					BMD.muladdin(beta,C,alpha,A,B) ;
				}

			}
		}
		else {
			if (tA) {
				if (tB) {
					BMD.muladd(D,beta,C,alpha,At,Bt) ; // D = alphaAB+beta C
				}
				else{
					BMD.muladd(D,beta,C,alpha,At,B) ;
				}
			}
			else {
				if (tB) {
					BMD.muladd(D,beta,C,alpha,A,Bt) ;
				}
				else{
					BMD.muladd(D,beta,C,alpha,A,B) ;
				}

			}
		}

		chrono.stop() ; fgemm_scal_tim += chrono ;
	}
	if (!j) {
		std::cout << "multiplication did not happen" << std::endl;
	}
#ifdef _LB_DEBUG
	else {
		std::cout << point_nb << std::endl;
	}
#endif
	mflops = compute_mflops(fgemm_scal_tim,fgemm_mflops(m,k,n),j);
	Data.setEntry(0,point_nb,mflops);
	std::ostringstream nam ;
	nam << "\"(" << m << ',' << k << ',' << n << ") ";
	nam << (inplace?"C":"D") << "=" << alpha << " "  ;
	nam << (tA?"t":"")<< "A." <<  (tB?"t":"")<< "B+" << beta << " C\"" ;
	Data.setAbsciName(point_nb,nam.str());
}


/* ********************** */
/*        Tests           */
/* ********************** */

/*! Benchmark fgemm Y=AX for several sizes of sqare matrices.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_blas( index_t min, index_t max, int step, int charac )
{
	int nb = 1 ;// une col de plus (la première)
	typedef LinBox::Modular<double>          Field0 ; ++nb ;
	typedef LinBox::Modular<float>           Field1 ; ++nb ;
	typedef LinBox::Modular<int32_t>         Field2 ; ++nb ;
	typedef LinBox::ModularBalanced<double>  Field3 ; ++nb ;
	typedef LinBox::ModularBalanced<float>   Field4 ; ++nb ;
	typedef LinBox::ModularBalanced<int32_t> Field5 ; ++nb ;
	typedef LinBox::UnparametricField<double>Field6 ; ++nb ;
	typedef LinBox::UnparametricField<float> Field7 ; ++nb ;

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	int it = 0 ;

	Field0 F0(charac) ;
	launch_bench_blas(F0,min,max,step,Data,it++);
	showFinish(it,nb);
	if (charac < 2048) {
		Field1 F1(charac) ;
		launch_bench_blas(F1,min,max,step,Data,it++);
		showFinish(it,nb);
	}
	else {
		showSkip(it,nb);
	}
	Field2 F2(charac) ;
	launch_bench_blas(F2,min,max,step,Data,it++);
	showFinish(it,nb);
	Field3 F3(charac) ;
	launch_bench_blas(F3,min,max,step,Data,it++);
	showFinish(it,nb);
	if (charac < 2048) {
		Field4 F4(charac) ;
		launch_bench_blas(F4,min,max,step,Data,it++);
	showFinish(it,nb);
	}
	else {
		showSkip(it,nb);
	}
	Field5 F5(charac) ;
	launch_bench_blas(F5,min,max,step,Data,it++);
	showFinish(it,nb);
	Field6 F6(charac) ;
	launch_bench_blas(F6,min,max,step,Data,it++);
	showFinish(it,nb);
	if (charac < 2048) {
		Field7 F7((float)charac) ;
		launch_bench_blas(F7,min,max,step,Data,it++);
		showFinish(it,nb);
	}
	else {
		showSkip(it,nb);
	}

	linbox_check(it <= nb);


	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::pdf);
	// Style.setTerm(LinBox::PlotStyle::png);
	// Style.setTerm(LinBox::PlotStyle::svg);
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("FFLAS::fgemm","Mflops","dimensions");
	// Style.setType(LinBox::PlotStyle::histogram);
	// Style.setStyle("set style histogram cluster gap 1");
	// Style.addStyle("set style fill solid border -1");
	// Style.addStyle("set boxwidth 0.9");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,nb));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("fgemm_blas");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark square fgemm Y=AX for several fields.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_square( index_t min, index_t max, int step, int charac )
{

	int nb = 1 ;// une col de plus (la première)
	typedef LinBox::Modular<double>          Field0 ; ++nb ;
	typedef LinBox::Modular<float>           Field1 ; ++nb ;
	typedef LinBox::Modular<int32_t>         Field2 ; ++nb ;
	typedef LinBox::ModularBalanced<double>  Field3 ; ++nb ;
	typedef LinBox::ModularBalanced<float>   Field4 ; ++nb ;
	typedef LinBox::ModularBalanced<int32_t> Field5 ; ++nb ;
	// GivaroZpZ

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	int it = 0 ;

	Field0 F0(charac) ;
	launch_bench_square(F0,min,max,step,Data,it++);
	showFinish(it,nb);
	if (charac < 2048) {
		Field1 F1(charac) ;
		launch_bench_square(F1,min,max,step,Data,it++);
		showFinish(it,nb);
	}
	else {
		showSkip(it,nb);
	}
	Field2 F2(charac) ;
	launch_bench_square(F2,min,max,step,Data,it++);
	showFinish(it,nb);
	Field3 F3(charac) ;
	launch_bench_square(F3,min,max,step,Data,it++);
	showFinish(it,nb);
	if (charac < 2048) {
		Field4 F4(charac) ;
		launch_bench_square(F4,min,max,step,Data,it++);
	showFinish(it,nb);
	}
	else {
		showSkip(it,nb);
	}
	Field5 F5(charac) ;
	launch_bench_square(F5,min,max,step,Data,it++);
	showFinish(it,nb);
	linbox_check(it <= nb);

	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::Term::pdf);
	// Style.setTerm(LinBox::PlotStyle::Term::png);
	// Style.setTerm(LinBox::PlotStyle::Term::svg);
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("BlasMatrixDomain mul","Mflops","dimensions");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,nb));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("bmdmul_square");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark fgemm Y=AX for several shapes.
 * Let n=k^2.
 * we test the following shapes :
 * - (l,nk,nk), (nk,l,nk), (nk,nk,l) : like vector-product
 * - (kl,nk,n), (nk,kl,n),(nk,n,kl)     : one small rectangular matrix
 * - (kl,n,nk), (n,kl,nk),(n,nk,kl)     : same
 * - (nl,n,n),(n,nl,n),(n,n,nl)         : square (or close to)
 * .
 * @param k parameter.
 * @param charac characteristic of the field.
 * @param l small parameter (ie close to 1)
 */
void bench_rectangular( index_t k, int charac, index_t l = 2 )
{
	int n  = k*k ;
	int nk = n*k ;
	int kl = k*l ;
	int nl = n*l ;
	typedef LinBox::Modular<double> Field ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 12 ;
	LinBox::PlotData<std::string>  Data(nb,1);
	Data.setSerieName(0,"mflops");
	launch_bench_rectangular(F,l,nk,nk,Data,it++);
	launch_bench_rectangular(F,nk,l,nk,Data,it++);
	launch_bench_rectangular(F,nk,nk,l,Data,it++);

	launch_bench_rectangular(F,kl,nk,n,Data,it++);
	launch_bench_rectangular(F,nk,kl,n,Data,it++);
	launch_bench_rectangular(F,nk,n,kl,Data,it++);

	launch_bench_rectangular(F,n,nk,kl,Data,it++);
	launch_bench_rectangular(F,nk,n,kl,Data,it++);
	launch_bench_rectangular(F,nk,kl,n,Data,it++);

	launch_bench_rectangular(F,nl,n,n,Data,it++);
	launch_bench_rectangular(F,n,nl,n,Data,it++);
	launch_bench_rectangular(F,n,n,nl,Data,it++);
	//!@todo resize if it>nb !!

	linbox_check(it==nb);

	LinBox::PlotStyle Style;
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("fgemm","x","y");
	Style.setPlotType(LinBox::PlotStyle::Plot::histo);
	Style.setXtics(LinBox::PlotStyle::Options::oblique);// make long legends oblique.
	Style.addPlotType("set style histogram cluster gap 1");
	Style.addPlotType("set style fill solid border -1");
	Style.addPlotType("set boxwidth 0.9");


	// Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(2);

	LinBox::PlotGraph<std::string> Graph(Data,Style);
	Graph.setOutFilename("fgemm_rect");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark fgemm \f$D\gets\alpha A B+\beta C\f$ for general \f$\alpha,\beta\f$.
 * @param k parameter.
 * @param charac characteristic of the field.
 * @param inplace "inplace" matmul (ie. C=D is overwritten)
 */
void bench_scalar( index_t k, int charac, bool inplace )
{
	typedef LinBox::Modular<double> Field ;
	typedef Field::Element Element;
	typedef Field::RandIter Randiter ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 12 ;
	LinBox::PlotData<std::string>  Data(nb,1);
	Data.setSerieName(0,"fgemm");

	Element one, zero, mOne, alpha, beta ;
	F.init(one,1);
	F.init(mOne,-1);
	F.init(zero,0);


	Randiter R(F) ;
	linbox_check(charac >=5) ;

	do { R.random(alpha) ; } // non trivial alpha
	while (F.areEqual(alpha,one) || F.areEqual(alpha,mOne) || F.areEqual(alpha,zero)) ;

	do { R.random(beta) ; }// non trivial beta
	while (F.areEqual(beta,one) || F.areEqual(beta,mOne) || F.areEqual(beta,zero)) ;


	// D = AB + beta C
	launch_bench_scalar<Field,0,0>(F,k,k,k,one,zero,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,one,mOne,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,one,one ,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,one,beta,Data,it++,inplace);

	// D = -AB + beta C
	launch_bench_scalar<Field,0,0>(F,k,k,k,mOne,zero,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,mOne,mOne,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,mOne,one ,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,mOne,beta,Data,it++,inplace);

	// D = alpha AB + beta C
	launch_bench_scalar<Field,0,0>(F,k,k,k,alpha,zero,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,alpha,mOne,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,alpha,one ,Data,it++,inplace);
	launch_bench_scalar<Field,0,0>(F,k,k,k,alpha,beta,Data,it++,inplace);


	linbox_check(it==nb);

	LinBox::PlotStyle Style;
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("fgemm","x","y");
	Style.setPlotType(LinBox::PlotStyle::Plot::histo);
	Style.setXtics(LinBox::PlotStyle::Options::oblique);// make long legends oblique.
	Style.addPlotType("set style histogram cluster gap 1");
	Style.addPlotType("set style fill solid border -1");
	Style.addPlotType("set boxwidth 0.9");

	// Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(2);

	LinBox::PlotGraph<std::string> Graph(Data,Style);
	Graph.setOutFilename("fgemm_scal");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark fgemm \f$D\gets\alpha A^x B^y+\beta C\f$ for \f$x,y=1,\top\f$ (transpose or not).
 * @param k parameter.
 * @param charac characteristic of the field.
 * @param inplace "inplace" matmul (ie. C=D is overwritten)
 */
void bench_transpose( index_t k, int charac, bool inplace )
{
	typedef LinBox::Modular<double> Field ;
	typedef Field::Element Element;
	typedef Field::RandIter Randiter ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 8 ;
	LinBox::PlotData<std::string>  Data(nb,1);
	Data.setSerieName(0,"fgemm");

	Element one, zero, mOne, alpha, beta ;
	F.init(one,1);
	F.init(mOne,-1);
	F.init(zero,0);


	Randiter R(F) ;
	linbox_check(charac >=5) ;

	do { R.random(alpha) ; } // non trivial alpha
	while (F.areEqual(alpha,one) || F.areEqual(alpha,mOne) || F.areEqual(alpha,zero)) ;

	do { R.random(beta) ; }// non trivial beta
	while (F.areEqual(beta,one) || F.areEqual(beta,mOne) || F.areEqual(beta,zero)) ;


	// D = A^xB^y
	launch_bench_scalar<Field,0,0>(F,k,k,k,one,zero,Data,it++,inplace);
	launch_bench_scalar<Field,1,0>(F,k,k,k,one,zero,Data,it++,inplace);
	launch_bench_scalar<Field,0,1>(F,k,k,k,one,zero,Data,it++,inplace);
	launch_bench_scalar<Field,1,1>(F,k,k,k,one,zero,Data,it++,inplace);

	// D = alpha A^xB^y + beta C
	launch_bench_scalar<Field,0,0>(F,k,k,k,alpha,beta,Data,it++,inplace);
	launch_bench_scalar<Field,1,0>(F,k,k,k,alpha,beta,Data,it++,inplace);
	launch_bench_scalar<Field,0,1>(F,k,k,k,alpha,beta,Data,it++,inplace);
	launch_bench_scalar<Field,1,1>(F,k,k,k,alpha,beta,Data,it++,inplace);


	linbox_check(it==nb);

	LinBox::PlotStyle Style;
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("fgemm","x","y");
	Style.setPlotType(LinBox::PlotStyle::Plot::histo);
	Style.setXtics(LinBox::PlotStyle::Options::oblique);// make long legends oblique.
	Style.addPlotType("set style histogram cluster gap 1");
	Style.addPlotType("set style fill solid border -1");
	Style.addPlotType("set boxwidth 0.9");

	// Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(2);

	LinBox::PlotGraph<std::string> Graph(Data,Style);
	Graph.setOutFilename("fgemm_trans");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}


/*  Benchmark fgemm Y = a AX + b Y for various (a,b) couples */

/*  main */

int main( int ac, char ** av)
{
	/*  Argument parsing/setting */

	static index_t       min  = 100;     /*  min size */
	static index_t       max  = 1500;    /*  max size (not included) */
	static index_t       step = 100;    /*  step between 2 sizes */
	static std::list<int> lst  ;       /*  what bench to start ? */
	lst.push_front(1);// ={1,2} vivement le nouveau std...
	lst.push_front(2);

	static Argument as[] = {
		{ 'm', "-m min" , "Set minimal size of matrix to test."    , TYPE_INT , &min },
		{ 'M', "-M Max" , "Set maximal size."                      , TYPE_INT , &max },
		{ 's', "-s step", "Sets the gap between two matrix sizes.", TYPE_INT , &step },
		{ 'l', "-l list", "Only launches a subset of available benchmarks\n - 1: compare to raw blas\n - 2:various square sizes\n - 3:various shapes\n - 4: various parameters (a,b)\n - 5 : various transp. combinations", TYPE_INTLIST, &lst },
		END_OF_ARGUMENTS
	};

	LinBox::parseArguments (ac, av, as);

	if (min >= max) {
		throw LinBox::LinBoxError("min value should be smaller than max...");
	}
	if (min + step >= max) {
		std::cout << "Warning : your x axis has only one point. You should have a smaller step." << std::endl;
	}

	//! @todo use commentator.
	lst.unique();
	lst.sort();
	if (lst.empty()) {
		std::cerr << "Warning, you are not benchmarking anything. Please check the -l value." << std::endl;
	}
	std::list<int>::iterator it = lst.begin();

	/* against pure blas routine */
	if (*it == 1) {
		std::cout << "bench 1 : blas" << std::endl;
		bench_blas(min,max,step,13);
		bench_blas(min,max,step,2011);
		bench_blas(min,max,step,65537);
		if (++it == lst.end()) return EXIT_SUCCESS ;
	}

	/* square for various fields */
	if (*it == 2) {
		std::cout << "bench 2 : square" << std::endl;
		bench_square(min,max,step,13);
		bench_square(min,max,step,2011);
		bench_square(min,max,step,65537);
		if (++it == lst.end()) return EXIT_SUCCESS ;
	}

	/* various shapes. */
	if (*it == 3) {
		std::cout << "bench 3 : shapes" << std::endl;
		int cube = (int) std::pow(max,double(1/3.));
		bench_rectangular(cube,2011);
		if (++it == lst.end()) return EXIT_SUCCESS ;
	}

	/* various parameters */
	if (*it == 4) {
		std::cout << "bench 4 : scalars" << std::endl;
		bench_scalar(max,65537,false);
		bench_scalar(max,65537,true);
		if (++it == lst.end()) return EXIT_SUCCESS ;
	}

	/* various tranpositions */
	if (*it == 5) {
		std::cout << "bench 2 : transp" << std::endl;
		bench_transpose(max,65537,true);
		bench_transpose(max,65537,true);
		if (++it != lst.end())  {
			std::cerr << "*** error *** your list contains (at least) one unknown element : " << *it << '!' << std::endl;
		}
	}
	else {
		std::cerr << "*** error *** your list contains (at least) one unknown element : " << *it << '!' << std::endl;
	}
	return EXIT_SUCCESS ;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

