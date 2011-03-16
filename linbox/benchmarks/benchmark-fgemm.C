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
		Data.setAbsciName(l,i);
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

/*  Benchmark fgemm Y=AX for several sizes of sqare matrices */
void bench_square( index_t min, index_t max, int step )
{
	typedef LinBox::Modular<double> Field ;
	typedef LinBox::Modular<float> Gield ;

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	LinBox::PlotData<index_t>  Data(nb_pts,4);
	Field F(13) ;
	launch_bench_square(F,min,max,step,Data,0);
	Field G(65537) ;
	launch_bench_square(G,min,max,step,Data,1);
	Gield H(13) ;
	launch_bench_square(H,min,max,step,Data,2);
	Gield I(2011) ;
	launch_bench_square(I,min,max,step,Data,3);


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
	Style.setUsingSeries(std::pair<index_t,index_t>(2,5));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("fgemm_square");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}
/*  Benchmark domain against fgemm  */

/*  Benchmark fgemm Y=AX for several fields */

/*  Benchmark fgemm Y=AX for several shapes */

/*  Benchmark fgemm Y=AX against available BLAS cblas_dgemm */

/*  Benchmark fgemm Y = a AX + b Y for various (a,b) couples */

/*  main */

int main( int ac, char ** av)
{

	bench_square(100,1000,50);

}
