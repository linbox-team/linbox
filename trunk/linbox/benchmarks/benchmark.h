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

/*! @file   benchmarks/benchmark.h
 * @ingroup benchmarks
 * @brief   Common header to ease benchmarking routines.
 * We provide a class to easily populate and plot files that represent
 * benchmarks.
 *
 * We use <a href="http://www.gnuplot.info/">gnuplot</a> for the plotting part
 * or LaTeX to provide tables.  A minimum knowledge of \c Gnuplot is a plus but
 * the \c benchmark-* files should provide enough examples for creating
 * standard (not too fancy) plots.
 *
 * We fall back to plain latex tabulars when gnuplot is not available. We plot
 * graphs in files whose name is extended by a random 8-char string to avoid
 * écraser files.
 *
 */

#ifndef __LINBOX_benchmarks_benchmark_H
#define __LINBOX_benchmarks_benchmark_H

#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include "tests/test-common.h"
#include "linbox/util/timer.h"
#include <string>
#include <fstream>
#include <iomanip> // setprecision

typedef uint32_t index_t ;

/* ********************** */
/*    Plot structures     */
/* ********************** */

namespace LinBox
// namespace benchmark ?
{

	/*! @brief Represents a table of values to plot (2D).
	 * list of values are reprensented by vectors.  the table is a vector
	 * of these vectors.
	 *
	 * @warning NaN, inf are used as missing data. More genenally
	 * we could store data in strings.
         * @todo Allow for 'speed up against col X' style
	 */
	class PlotStyle ;

	/*! @brief The raw data to plot (2D).
	 * Represents the series of points and the labels for the points (X)
	 * axis and the series (Y) axis.  The internal representation is a
	 * vector of vector, each series of point being a vector of double.
	 * @tparam NAM the X axis is parametered by NAM (string, int, double...)
	 * @todo write members that permute, add, scale,... data.
	 */
	template<class NAM>
	class PlotData ;

	/*! @brief The graph (2D).
	 * This class joins a PlotStyle and a PlotData to build up a graph.  A
	 * filename should be provided as well, indicating where the output
	 * graph and scripts will be generated.
	 *
	 * @warning the filename will get a random suffix before the extension
	 * so as not to overwrite files "par inadvertance".
	 */
	template<class NAM>
	class PlotGraph ;

}

/* ********************** */
/*        Outils          */
/* ********************** */

using Givaro::Timer; // template by timer_t ?

namespace LinBox {

	/** Helper.
	 * This helper has several functions :
	 *   - Records the timings
	 *   - predict the execution time for the next experiment
	 *   - helps producing enough experiments (but not too much and not too time consuming) for producing a valid measure.
	 *   .
	 *   See member function help for more information.
	 */
	class TimeWatcher  ;

}// LinBox
	void showAdvanceLinear(int curr, int min, int max)
	{
		std::cout << std::setprecision(4) << "\033[2K" << "\033[30D" << min <<std::flush;
		std::cout << '<' << curr << '<' << max << " (" << std::flush;
		std::cout << double(curr-min)/double(max-min)*100 << "%)" << std::flush;
	}
	void showFinish(int curr, int all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "finished : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}
	void showSkip(int curr, int all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "skipped : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}


// double compute_mflops(const Timer & t, const double mflo, const int rpt = 1)

#include "benchmark.inl"

#endif // __LINBOX_benchmarks_benchmark_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
