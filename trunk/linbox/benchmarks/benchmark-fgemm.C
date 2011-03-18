/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* Copyright (C) 2011 LinBox
 * Written by BB <brice.boyer@imag.fr>
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/*! @file benchmarks/benchmark-fgemm.C
 * @ingroup benchmarks
 * @brief Benchmarking dense matrix multiplication on finite fields.
 * This file benchmarks the FFPACK::fgemm implementation for various fields,
 * shape and parameters. Actually, we use the wrapper member \c mul of LinBox::BlasMatrixDomain.
 */

#include "benchmarks/benchmark.h"
#include "linbox/util/error.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/ffpack/ffpack.h"
#include "linbox/fflas/fflas.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/matrix/blas-matrix.h"

/* ********************** */
/*        Outils          */
/* ********************** */

bool keepon(index_t & repet, Timer & tim)
{
	if (tim.userElapsedTime() < 0.5) {
		++repet ;
		return true;
	}
	if (repet<1)
		return true
	return false ;
}

double fgemm_mflops(int m, int n, int k)
{
	return (double)m*(double)n/1e6*(double)k ;
}

double compute_mflops(const Timer & t, const double mflo, const int rpt = 1)
{
	linbox_check(rpt);
	return (double) ((mflo*rpt)/t.usertime());
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
	index_t l = 0 ;
	Timer fgemm_sq_tim ;
	fgemm_sq_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	for ( index_t i = min ; i < max ; i += step , ++l ) {
		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		LinBox::BlasMatrix<Element> A (ii,ii);
		LinBox::BlasMatrix<Element> B (ii,ii);
		LinBox::BlasMatrix<Element> C (ii,ii);
		RandMat.random(A);
		RandMat.random(B);
		RandMat.random(C);
		if (!series_nb)
			Data.setAbsciName(l,i); // only write abscissa for serie 0
		index_t j = 0 ; // number of repets.
		// fgemm_sq_tim.clear() ;
		fgemm_sq_tim.start() ;
		while( keepon(j,fgemm_sq_tim) ) {
			BMD.mul(C,A,B) ; // C = AB
		}
		fgemm_sq_tim.stop();
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
template<class T>
void launch_bench_blas(index_t min, index_t max, int step
		       , LinBox::PlotData<index_t> & Data
		       , index_t series_nb
		       , index_t charact)
{
	typedef LinBox::Modular<T> Field ;
	Field F((int)charact);
	index_t l = 0 ;
	Timer fgemm_blas_tim ;
	fgemm_blas_tim.clear();
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
	for (index_t j = 0 ; j < mm ; ++j) R.random(A[j]);
	for (index_t j = 0 ; j < mm ; ++j) R.random(B[j]);
	for (index_t j = 0 ; j < mm ; ++j) R.random(C[j]);
	typedef typename LinBox::UnparametricField<T> FloatingDomain ;
	FloatingDomain G ;

	for ( index_t i = min ; i < max ; i += step , ++l ) {
		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		index_t j = 0 ;
		fgemm_blas_tim.clear() ;
		fgemm_blas_tim.start() ;
		while(keepon(j,fgemm_blas_tim)) {
			LinBox::FFLAS::fgemm(G,LinBox::FFLAS::FflasNoTrans,LinBox::FFLAS::FflasNoTrans,
					     ii,ii,ii,
					     1.,
					     A,ii,B,ii,
					     0.,
					     C,ii) ; // the last 2 arguments are superfluous and meaningless.
		}
		fgemm_blas_tim.stop();
		mflops = compute_mflops(fgemm_blas_tim,fgemm_mflops(i,i,i),j);
		Data.setEntry(series_nb,l,mflops);
	}

	delete[] A ;
	delete[] B ;
	delete[] C ;
	std::ostringstream nam ;
	nam << '\"' ;
	G.write(nam);
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
	Timer fgemm_rect_tim ;
	fgemm_rect_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	LinBox::BlasMatrix<Element> A (m,k);
	LinBox::BlasMatrix<Element> B (k,n);
	LinBox::BlasMatrix<Element> C (m,n);
	RandMat.random(A);
	RandMat.random(B);
	RandMat.random(C);
	index_t j = 0 ;
	fgemm_rect_tim.clear() ;
	fgemm_rect_tim.start() ;
	while (keepon(j,fgemm_rect_tim)) {
		BMD.mul(C,A,B) ; // C = AB
	}
	fgemm_rect_tim.stop();
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


/* ********************** */
/*        Tests           */
/* ********************** */


/*! Benchmark fgemm Y=AX for several sizes of sqare matrices.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_blas( index_t min, index_t max, int step )
{
	//!@todo compare to cblas_dgemm instead.
	typedef LinBox::Modular<double> DoubleField ;
	typedef LinBox::Modular<float>  FloatField ;

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	int it = 0 ; int nb = 5 ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	FloatField F0(13) ;
	launch_bench_square(F0,min,max,step,Data,it++);
	launch_bench_blas<float>(min,max,step,Data,it++,13);

	DoubleField F1(65537) ;
	launch_bench_square(F1,min,max,step,Data,it++);
	launch_bench_blas<double>(min,max,step,Data,it++,65537);

	linbox_check(it+1==nb);


	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::pdf);
	// Style.setTerm(LinBox::PlotStyle::png);
	// Style.setTerm(LinBox::PlotStyle::svg);
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("fgemm","x","y");
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
	if (charac < 2048) {
		Field1 F1(charac) ;
		launch_bench_square(F1,min,max,step,Data,it++);
	}
	Field2 F2(charac) ;
	launch_bench_square(F2,min,max,step,Data,it++);
	Field3 F3(charac) ;
	launch_bench_square(F3,min,max,step,Data,it++);
	if (charac < 2048) {
		Field4 F4(charac) ;
		launch_bench_square(F4,min,max,step,Data,it++);
	}
	Field5 F5(charac) ;
	launch_bench_square(F5,min,max,step,Data,it++);
	linbox_check(it <= nb);

	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::Term::pdf);
	// Style.setTerm(LinBox::PlotStyle::Term::png);
	// Style.setTerm(LinBox::PlotStyle::Term::svg);
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("fgemm","x","y");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,it));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("fgemm_square");

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
	Style.addPlotType("set style histogram cluster gap 1");
	Style.addPlotType("set style fill solid border -1");
	Style.addPlotType("set boxwidth 0.9");
	//! @todo make long legends oblique.

	// Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(2);

	LinBox::PlotGraph<std::string> Graph(Data,Style);
	Graph.setOutFilename("fgemm_rect");

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

	static size_t       min = 100;     /*  min size */
	static size_t       max = 1500;    /*  max size (not included) */
	static size_t       step = 100;    /*  step between 2 sizes */

        static Argument as[] = {
                { 'm', "-m min", "Set minimal size of matrix to test."    , TYPE_INT , &min },
                { 'M', "-M Max", "Set maximal size."                      , TYPE_INT , &max },
                { 's', "-s step", "Sets the gap between two matrix sizes.", TYPE_INT , &step },
		END_OF_ARGUMENTS
        };

	parseArguments (ac, av, as);

	if (min >= max) {
		throw LinBox::LinBoxError("min value should be smaller than max...");
	}
	if (min + step >= max) {
		std::cout << "Warning : your x axis has only one point. You should have a smaller step." << std::endl;
	}

	//! @todo use commentator.


	// against pure blas routine
	bench_blas(min,max,step);

	// square for various fields
	bench_square(min,max,step,13);
	bench_square(min,max,step,2011);
	bench_square(min,max,step,65537);

	// various shapes.
	int cube = (int) std::pow(max,double(1/3.));
	// std::cout << cube << std::endl;
	bench_rectangular(cube,2011);

	return EXIT_SUCCESS ;
}
