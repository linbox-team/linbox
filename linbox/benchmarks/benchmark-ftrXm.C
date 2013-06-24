
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

/*! @file benchmarks/benchmark-ftrXm.C
 * @ingroup benchmarks
 * @brief Benchmarking triangular matrix multiplication on finite fields.
 * This file benchmarks the FFLAS::ftrmm, FFLAS::ftrsm implementation for various fields,
 * shape and parameters. Actually, we use the wrapper member \c mul of LinBox::BlasMatrixDomain.
 * @todo ftrmm has an 'alpha' but mul/mulin in BMd don't... That could be useful for \f$\alpha=-1\f$...
 * @todo benchmark ftrsm too here.
 */

#include "benchmarks/benchmark.h"
#include "linbox/util/error.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
// parse arguments
#include "fflas-ffpack/utils/args-parser.h"


#define _LB_LEFT   true
#define _LB_RITE  false

#define _LB_TRANS  false
#define _LB_NOTRS  true

#define _LB_UNIT   true
#define _LB_DIAG   false

#define _LB_TSUP   true
#define _LB_TLOW   false

using LinBox::TimeWatcher;




double ftrmm_mflops(index_t m, index_t n, index_t k)
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
template<class Field, bool LeftSide, bool UnitDiag, bool TriSup>
void launch_bench_square(Field & F // const problem
			, index_t min, index_t max, int step
			, LinBox::PlotData<index_t> & Data
			, index_t series_nb)
{
	TimeWatcher TW(10,series_nb);
	index_t l = 0 ;
	Timer fgemm_sq_tim ;
	Timer chrono ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	LinBox::NonzeroRandIter<Field> Rn(F,R);
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);
	// index_t repet = 3 ;
	for ( index_t i = min ; i < max ; i += step , ++l ) {
	double mflops ;
		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		LinBox::TriangularBlasMatrix<Field> A (F,ii,ii,
							 (TriSup?LinBox::LinBoxTag::Upper:LinBox::LinBoxTag::Lower),
							 (UnitDiag?LinBox::LinBoxTag::Unit:LinBox::LinBoxTag::NonUnit));
		LinBox::BlasMatrix<Field> B (F,ii,ii);
		if (!series_nb)
			Data.setAbsciName(l,i); // only write abscissa for serie 0
		index_t j = 0 ; // number of repets.
		fgemm_sq_tim.clear() ;
		while( TW.keepon(j,fgemm_sq_tim) ) {
			RandMat.random(A);
			for (size_t k=0 ; k<(size_t)ii ; ++k) Rn.random(A.refEntry(k,k)) ;
			RandMat.random(B);
			chrono.clear() ; chrono.start() ;
			if(LeftSide)
				BMD.mulin_right(A,B) ; // B <- AB
			else
				BMD.mulin_left(B,A) ;  // B <- BA
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
		mflops = compute_mflops(fgemm_sq_tim,ftrmm_mflops(i,i,i),j);
		Data.setEntry(series_nb,l,mflops);
	}
	std::string nam = "\"";
	if (TriSup)
		nam += "upper " ;
	else
		nam += "lower " ;
	if (LeftSide)
		nam += "left " ;
	else
		nam += "right " ;
	if (UnitDiag)
		nam += "non-" ;
	nam += "unit\"" ;

	Data.setSerieName(series_nb,nam);

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
	TimeWatcher TW(10,series_nb);
	typedef LinBox::Modular<T> Field ;
	Field F((int)charact);
	index_t l = 0 ;
	Timer ftrmm_blas_tim ;
	Timer chrono ;
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
	typedef typename LinBox::UnparametricField<T> FloatingDomain ;
	FloatingDomain G ;

	for ( index_t i = min ; i < max ; i += step , ++l ) {
	double mflops ;
		int ii = i ; // sinon, le constructeur le plus proche serait (_Matrix,_Field)... n'impnawak...
		index_t mimi = (index_t) ii*ii ;
		index_t j = 0 ;
		ftrmm_blas_tim.clear() ;
		while(TW.keepon(j,ftrmm_blas_tim)) {
			for (index_t jl = 0 ; jl < mimi ; ++jl) R.random(A[jl]);
			for (index_t jl = 0 ; jl < mimi ; ++jl) R.random(B[jl]);
			for (index_t jl = 0 ; jl < mimi ; ++jl) R.random(C[jl]);
			chrono.clear(); chrono.start() ;
			FFLAS::ftrmm(G,FFLAS::FflasLeft,FFLAS::FflasUpper,
					     FFLAS::FflasNoTrans,
					     FFLAS::FflasUnit,
					     ii,ii,
					     1.,
					     A,ii, B,ii) ;
			chrono.stop() ;
			ftrmm_blas_tim += chrono ;
		}
		mflops = compute_mflops(ftrmm_blas_tim,ftrmm_mflops(i,i,i),j);
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
 * @param n cols in A
 * @param Data where data is stored
 * @param point_nb point to be computed
 */
template<class Field, bool LeftSide, bool UnitDiag, bool TriSup>
void launch_bench_rectangular(Field & F // const problem
			      , index_t m, index_t n
			      , LinBox::PlotData<std::string> & Data
			      , index_t point_nb)
{
	TimeWatcher TW(10,point_nb);
	Timer ftrmm_rect_tim ;
	Timer chrono ; chrono.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	typedef typename LinBox::TriangularBlasMatrix<Field> TriangularMatrix ;
	typedef typename LinBox::BlasMatrix<Field>  Matrix ;

	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);

	index_t k = (LeftSide?m:n);

	TriangularMatrix A (F,k,k,
			    (TriSup?LinBox::LinBoxTag::Upper:LinBox::LinBoxTag::Lower),
			    (UnitDiag?LinBox::LinBoxTag::Unit:LinBox::LinBoxTag::NonUnit));
	Matrix B (F,(int)m,(int)n);

	index_t j = 0 ;
	ftrmm_rect_tim.clear() ;
	while (TW.keepon(j,ftrmm_rect_tim)) {
		RandMat.random(A);
		RandMat.random(B);
		chrono.clear() ; chrono.start() ;

		if(LeftSide)
			BMD.mulin_right(A,B) ; // B <- AB
		else
			BMD.mulin_left(B,A) ;  // B <- BA

		chrono.stop();
		ftrmm_rect_tim += chrono ;
	}
	if (!j) {
		std::cout << "multiplication did not happen" << std::endl;
	}
#ifdef _LB_DEBUG
	else {
		std::cout << point_nb << std::endl;
	}
#endif
	mflops = compute_mflops(ftrmm_rect_tim,ftrmm_mflops(m,k,n),(int)j);
	Data.setEntry(0,point_nb,mflops);
	std::ostringstream nam ;
	if (LeftSide)
		nam << "\"(" << k << ':' << m << ',' << n << ")" ;
	else
		nam << "\"(" << m << ',' << n << ':' << k << ")" ;
	nam << " (" << (UnitDiag?"":"non") << "unit) on ";
	F.write(nam);
	nam << '\"' ;


	Data.setAbsciName(point_nb,nam.str());
	return ;
}

/*! @internal
 * @brief launches the benchmarks for various parameters of a, b.
 * B = aAB or B=aBA with A triangular.
 * @todo We test various parameters \f$alpha\f$.
 * We test A being transposed or not.
 * @param F field
 * @param m rows in B
 * @param n cols in B
 * @param alpha alpha (not zero)
 * @param Data where data is stored
 * @param point_nb point to be computed
 */
template<class Field, bool LeftSide, bool UnitDiag, bool TriSup, bool tA>
void launch_bench_scalar(Field & F // const problem
			 , int m, int n
			 , const typename Field::Element & alpha //!@warning not used yet.
			 , LinBox::PlotData<std::string> & Data
			 , index_t point_nb)
{
	TimeWatcher TW(10,point_nb);
	Timer ftrmm_scal_tim ;
	Timer chrono ;
	ftrmm_scal_tim.clear();
	double mflops ;
	typedef typename Field::Element  Element;
	typedef typename Field::RandIter Randiter ;
	typedef typename LinBox::BlasMatrix<Field>                    Matrix ;
	typedef typename LinBox::TriangularBlasMatrix<Field>           TriangularMatrix ;
	typedef typename LinBox::TransposedBlasMatrix<TriangularMatrix > TransposedTriangular ;


	Randiter R(F) ;
	LinBox::BlasMatrixDomain<Field> BMD(F) ;
	LinBox::RandomDenseMatrix<Randiter,Field> RandMat(F,R);

	index_t k = (LeftSide?m:n);

	TriangularMatrix A (F,k,k,
			    (TriSup?LinBox::LinBoxTag::Upper:LinBox::LinBoxTag::Lower),
			    (UnitDiag?LinBox::LinBoxTag::Unit:LinBox::LinBoxTag::NonUnit));
	TransposedTriangular At(A);
	Matrix B (F,m,n);

	index_t j = 0 ;
	while (TW.keepon(j,ftrmm_scal_tim)) {
		RandMat.random(A);
		RandMat.random(B);

		chrono.clear() ; chrono.start();

		if (tA) {
			if(LeftSide)
				BMD.mulin_right(At,B) ; // B <- AB
			else
				BMD.mulin_left(B,At) ;  // B <- BA
		}
		else  {
			if(LeftSide)
				BMD.mulin_right(A,B) ; // B <- AB
			else
				BMD.mulin_left(B,A) ;  // B <- BA
		}


		chrono.stop() ; ftrmm_scal_tim += chrono ;
	}
	if (!j) {
		std::cout << "multiplication did not happen" << std::endl;
	}
#ifdef _LB_DEBUG
	else {
		std::cout << point_nb << std::endl;
	}
#endif
	mflops = compute_mflops(ftrmm_scal_tim,ftrmm_mflops(m,k,n),(int)j);
	Data.setEntry(0,point_nb,mflops);
	std::ostringstream nam ;
	nam << "\"(" << m << 'x' << n << ") ";
	nam << "B=" << alpha << " " ;
	if (LeftSide)
		nam <<  (tA?"t":"")<< "A.B" ;
	else
		nam <<  "B." << (tA?"t":"")<< "A" ;
	nam << " (" << (UnitDiag?"":"non") << "unit) on ";
	F.write(nam);
	nam << '\"' ;

	Data.setAbsciName(point_nb,nam.str());
}

/* ********************** */
/*        Tests           */
/* ********************** */

/*! @brief Benchmark square ftrmm for differenct parameters.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 * @param charac characteristic of the field.
 */
void bench_square( index_t min, index_t max, int step, int charac )
{

	int nb = 9 ;// une col de plus (la première)
	typedef LinBox::Modular<double>          Field ; ++nb ;

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	LinBox::PlotStyle Style;
	int it = 0 ;
	Field F(charac) ;
	std::cout << "0.." << std::flush;
	/* _LB_LEFT,_LB_UNIT,_LB_TSUP */
	launch_bench_square<Field,true,true,true>(F,min,max,step,Data,it++);
	std::cout << "1.." << std::flush;
	launch_bench_square<Field,true,true,false>(F,min,max,step,Data,it++);
	std::cout << "2.." << std::flush;
	launch_bench_square<Field,true,false,true>(F,min,max,step,Data,it++);
	std::cout << "3.." << std::flush;
	launch_bench_square<Field,true,false,false>(F,min,max,step,Data,it++);
	std::cout << "4.." << std::flush;
	launch_bench_square<Field,false,true,true>(F,min,max,step,Data,it++);
	std::cout << "5.." << std::flush;
	launch_bench_square<Field,false,true,false>(F,min,max,step,Data,it++);
	std::cout << "6.." << std::flush;
	launch_bench_square<Field,false,false,true>(F,min,max,step,Data,it++);
	std::cout << "7.." << std::flush;
	launch_bench_square<Field,false,false,false>(F,min,max,step,Data,it++);
	std::cout << "9!!" << std::endl;



	linbox_check(it <= nb);

	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("ftrmm","x","y");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,nb));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("ftrmm_square");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark square ftrmm for several fields.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
 */
void bench_fields( index_t min, index_t max, int step )
{

	int nb = 1 ;// une col de plus (la première)
	typedef LinBox::Modular<float>                 Field0 ; ++nb ;
	typedef LinBox::Modular<double>                Field1 ; ++nb ;
	typedef LinBox::ModularBalanced<float>         Field2 ; ++nb ;
	typedef LinBox::ModularBalanced<double>        Field3 ; ++nb ;

	Field0 F0(2011);
	Field1 F1(65537);
	Field2 F2(2011);
	Field3 F3(65537);

	int nb_pts = (int) std::ceil((double)(max-min)/(double)step) ;
	LinBox::PlotData<index_t>  Data(nb_pts,nb);
	LinBox::PlotStyle Style;
	int it = 0 ;
	std::cout << "0.." << std::flush;
	/* _LB_LEFT,_LB_UNIT,_LB_TSUP */
	launch_bench_square<Field0,true,true,true>(F0,min,max,step,Data,it++);
	std::cout << "1.." << std::flush;
	launch_bench_square<Field1,true,true,true>(F1,min,max,step,Data,it++);
	std::cout << "2.." << std::flush;
	launch_bench_square<Field2,true,true,true>(F2,min,max,step,Data,it++);
	std::cout << "3.." << std::flush;
	launch_bench_square<Field3,true,true,true>(F3,min,max,step,Data,it++);
	std::cout << "4!!" << std::endl;



	linbox_check(it <= nb);

	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("ftrmm for various fields","x","y");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,nb));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("ftrmm_fields");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}


/*! Benchmark fgemm Y=AX for several sizes of sqare matrices.
 * @param min min size
 * @param max max size
 * @param step step of the size between 2 benchmarks
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
	/* _LB_LEFT,_LB_UNIT,_LB_TSUP */
	launch_bench_square<FloatField,true,true,true>(F0,min,max,step,Data,it++);
	launch_bench_blas<float>(min,max,step,Data,it++,13);

	DoubleField F1(65537) ;
	/* _LB_LEFT,_LB_UNIT,_LB_TSUP */
	launch_bench_square<DoubleField,true,true,true>(F1,min,max,step,Data,it++);
	launch_bench_blas<double>(min,max,step,Data,it++,65537);

	linbox_check(it+1==nb);


	LinBox::PlotStyle Style;
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("ftrmm","x","y");

	Style.setPlotType(LinBox::PlotStyle::Plot::graph);
	Style.setLineType(LinBox::PlotStyle::Line::linespoints);
	Style.setUsingSeries(std::pair<index_t,index_t>(2,nb));

	LinBox::PlotGraph<index_t> Graph(Data,Style);
	Graph.setOutFilename("ftrmm_blas");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}


/*! @brief Benchmark fgemm Y=AX for several shapes.
 * we test the following shapes :
 *   - (2k,2k):2k, 2k:(2k,2k)
 *   - (k,2k):2k, 2k:(2k,k)
 *   - (4k,2k):2k , 2k:(2k,4k)
 *   .
 * @param k parameter.
 * @param charac characteristic of the field.
 */
void bench_rectangular( index_t k, int charac )
{
	typedef LinBox::Modular<double> Field ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 6 ;
	LinBox::PlotData<std::string>  Data(nb,1);
	Data.setSerieName(0,"mflops");

	launch_bench_rectangular<Field,_LB_LEFT,_LB_UNIT,_LB_TSUP>(F,2*k,2*k,Data,it++);
	launch_bench_rectangular<Field,_LB_RITE,_LB_UNIT,_LB_TSUP>(F,2*k,2*k,Data,it++);

	launch_bench_rectangular<Field,_LB_LEFT,_LB_UNIT,_LB_TSUP>(F,2*k,k,Data,it++);
	launch_bench_rectangular<Field,_LB_RITE,_LB_UNIT,_LB_TSUP>(F,k,2*k,Data,it++);

	launch_bench_rectangular<Field,_LB_LEFT,_LB_UNIT,_LB_TSUP>(F,2*k,4*k,Data,it++);
	launch_bench_rectangular<Field,_LB_RITE,_LB_UNIT,_LB_TSUP>(F,4*k,2*k,Data,it++);

	//!@todo resize if it>nb !!

	linbox_check(it==nb);

	LinBox::PlotStyle Style;
	Style.setTerm(LinBox::PlotStyle::Term::eps);
	Style.setTitle("ftrmm","x","y");
	Style.setPlotType(LinBox::PlotStyle::Plot::histo);
	Style.setXtics(LinBox::PlotStyle::Options::oblique);// make long legends oblique.
	Style.addPlotType("set style histogram cluster gap 1");
	Style.addPlotType("set style fill solid border -1");
	Style.addPlotType("set boxwidth 0.9");
	//! @todo make long legends oblique.

	// Style.setType(LinBox::PlotStyle::lines);
	Style.setUsingSeries(2);

	LinBox::PlotGraph<std::string> Graph(Data,Style);
	Graph.setOutFilename("ftrmm_rect");

	// Graph.plot();

	Graph.print_gnuplot();

	Graph.print_latex();

	return ;

}

/*! @brief Benchmark ftrmm with transpose or alpha parameters on.
 * @param k parameter.
 * @param charac characteristic of the field.
 */
void bench_transpose( index_t k, int charac)
{
	typedef LinBox::Modular<double> Field ;
	typedef Field::Element Element;
	typedef Field::RandIter Randiter ;
	Field F(charac) ;

	index_t it = 0 ; index_t nb = 2 ;
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
	launch_bench_scalar<Field,true,true,true,true>(F,k,k,one,Data,it++);
	launch_bench_scalar<Field,true,true,true,false>(F,k,k,one,Data,it++);

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

/*  main */

int main( int ac, char ** av)
{
	/*  Argument parsing/setting */

	static index_t  min = 100;     /*  min size */
	static index_t  max = 2000;    /*  max size (not included) */
	static index_t step = 200;    /*  step between 2 sizes */
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

	bench_square(min,max,step,2011);
	bench_fields(min,max,step);
	bench_blas(min,max,step);
	bench_rectangular(max/2,2011);
	bench_transpose(max,2011);

#if 0 /*  to be uncommented */
	if (lst.empty()) {
		std::cerr << "Warning, you are not benchmarking anything. Please check the -l value." << std::endl;
	}
	std::list<int>::iterator it = lst.begin();

	/* against pure blas routine */
	if (*it == 1) {
		std::cout << "bench 1 : blas" << std::endl;
		bench_blas(min,max,step);
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
#endif
	return EXIT_SUCCESS ;
}

#undef _LB_RITE
#undef _LB_LEFT
#undef _LB_TSUP
#undef _LB_TLOW
#undef _LB_UNIT
#undef _LB_DIAG
#undef _LB_TRANS
#undef _LB_NOTRS

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

