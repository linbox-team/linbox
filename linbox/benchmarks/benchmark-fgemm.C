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
#include "linbox/util/timer.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/ffpack/ffpack.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/matrix/blas-matrix.h"

double fgemm_mflops(int m, int n, int k)
{
	return (double)m*(double)n/1e6*(double)k ;
}

double compute_mflops(Timer & t, double mflo, int rpt = 1)
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
	Timer fgemm_tim ;
	fgemm_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	index_t repet = 3 ;
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
		int temoin = 0 ; // makes sure mul is computed at almost no cost.
		index_t j = 0 ;
		fgemm_tim.clear() ;
		fgemm_tim.start() ;
		for (; j < repet ; ++j ) {
			BMD.mul(C,A,B) ; // C = AB
			temoin += (int)std::abs(C.getEntry(0,0))+1 ;
			if (fgemm_tim.usertime()>0.5) {
				++j;
				break;
			}
		}
		fgemm_tim.stop();
		if (temoin == 0) {
			std::cout << "multiplication did not happen" << std::endl;
		}
#ifdef _LB_DEBUG
		else {
			std::cout << i << ',' << j << std::endl;
		}
#endif
		mflops = compute_mflops(fgemm_tim,fgemm_mflops(i,i,i),j);
		Data.setEntry(series_nb,l,mflops);
	}
	std::ostringstream nam ;
	nam << '\"' ;
	F.write(nam);
	nam << '\"' ;
	Data.setSerieName(series_nb,nam.str());

}

/*! Benchmark fgemm Y=AX for several sizes of sqare matrices.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_square( index_t min, index_t max, int step )
{
	typedef LinBox::Modular<double> Field ;
	typedef LinBox::Modular<float> Gield ;

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	int it = 0 ; int nb = 4 ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	Field F(13) ;
	launch_bench_square(F,min,max,step,Data,it++);
	Field G(2011) ;
	launch_bench_square(G,min,max,step,Data,it++);
	Gield H(13) ;
	launch_bench_square(H,min,max,step,Data,it++);
	Gield I(2011) ;
	launch_bench_square(I,min,max,step,Data,it++);
	linbox_check(it==nb);


	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::pdf);
	// Style.setTerm(LinBox::PlotStyle::png);
	// Style.setTerm(LinBox::PlotStyle::svg);
	Style.setTerm(LinBox::PlotStyle::eps);
	Style.setTitle("fgemm","x","y");
	// Style.setType(LinBox::PlotStyle::histogram);
	// Style.setStyle("set style histogram cluster gap 1");
	// Style.addStyle("set style fill solid border -1");
	// Style.addStyle("set boxwidth 0.9");

	Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,it));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("fgemm_square");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*!  Benchmark fgemm Y=AX for several fields.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_square2( index_t min, index_t max, int step, int charac )
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
	// Style.setTerm(LinBox::PlotStyle::pdf);
	// Style.setTerm(LinBox::PlotStyle::png);
	// Style.setTerm(LinBox::PlotStyle::svg);
	Style.setTerm(LinBox::PlotStyle::eps);
	Style.setTitle("fgemm","x","y");
	// Style.setType(LinBox::PlotStyle::histogram);
	// Style.setStyle("set style histogram cluster gap 1");
	// Style.addStyle("set style fill solid border -1");
	// Style.addStyle("set boxwidth 0.9");

	Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,it));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("fgemm_square2");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

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
	Timer fgemm_tim ;
	fgemm_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	index_t repet = 3 ;
	LinBox::BlasMatrix<Element> A (m,k);
	LinBox::BlasMatrix<Element> B (k,n);
	LinBox::BlasMatrix<Element> C (m,n);
	RandMat.random(A);
	RandMat.random(B);
	RandMat.random(C);
	int temoin = 0 ; // makes sure mul is computed at almost no cost.
	index_t j = 0 ;
	fgemm_tim.clear() ;
	fgemm_tim.start() ;
	for (; j < repet ; ++j ) {
		BMD.mul(C,A,B) ; // C = AB
		temoin += (int)std::abs(C.getEntry(0,0))+1 ;
		if (fgemm_tim.usertime()>0.5) {
			++j;
			break;
		}
	}
	fgemm_tim.stop();
	if (temoin == 0) {
		std::cout << "multiplication did not happen" << std::endl;
	}
#ifdef _LB_DEBUG
	else {
		std::cout << point_nb << std::endl;
	}
#endif
	mflops = compute_mflops(fgemm_tim,fgemm_mflops(m,k,n),j);
	Data.setEntry(0,point_nb,mflops);
	std::ostringstream nam ;
	nam << "\"(" << m << ',' << k << ',' << n << ")\"" ;
	Data.setAbsciName(point_nb,nam.str());
}


/*!  Benchmark fgemm Y=AX for several shapes.
 * Let n=k^2.
 * we test the following shapes :
 * (1,nk,nk), (nk,1,nk), (nk,nk,1),
 * (k,nk,n), (nk,k,n),(nk,n,k)
 * (k,n,nk), (n,k,nk),(n,nk,k)
 * (n,n,n),
 * @param k parameter.
 */
void bench_rectangular( index_t k )
{
	int n  = k*k ;
	int nk = n*k ;
	int charac = 2011 ;
	typedef LinBox::Modular<double> Field ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 10 ;
	LinBox::PlotData<std::string>  Data(nb,1);
	Data.setSerieName(0,"mflops");
	launch_bench_rectangular(F,1,nk,nk,Data,it++);
	launch_bench_rectangular(F,nk,1,nk,Data,it++);
	launch_bench_rectangular(F,nk,nk,1,Data,it++);

	launch_bench_rectangular(F,k,nk,n,Data,it++);
	launch_bench_rectangular(F,nk,k,n,Data,it++);
	launch_bench_rectangular(F,nk,n,k,Data,it++);

	launch_bench_rectangular(F,n,nk,k,Data,it++);
	launch_bench_rectangular(F,nk,n,k,Data,it++);
	launch_bench_rectangular(F,nk,k,n,Data,it++);

	launch_bench_rectangular(F,n,n,n,Data,it++);

	linbox_check(it==nb);

	LinBox::PlotStyle Style;
	// Style.setTerm(LinBox::PlotStyle::pdf);
	// Style.setTerm(LinBox::PlotStyle::png);
	// Style.setTerm(LinBox::PlotStyle::svg);
	Style.setTerm(LinBox::PlotStyle::eps);
	Style.setTitle("fgemm","x","y");
	Style.setType(LinBox::PlotStyle::histogram);
	Style.setStyle("set style histogram cluster gap 1");
	Style.addStyle("set style fill solid border -1");
	Style.addStyle("set boxwidth 0.9");
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

/*  Benchmark fgemm Y=AX against available BLAS cblas_dgemm */

/*  Benchmark fgemm Y = a AX + b Y for various (a,b) couples */

/*  main */

int main( int ac, char ** av)
{

	bench_square(100,1500,100);
	bench_square2(100,1500,100,13);
	// bench_square2(100,1500,100,2011);
	bench_square2(100,1500,100,65537);

	bench_rectangular(20);
}
