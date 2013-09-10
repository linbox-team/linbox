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

	/*! @brief The raw data to plot.
	 * Represents  the labels for the points (X axis) and the values for
	 * each series of measures (Y axis).
	 *
	 * @internal The internal representation is a
	 * vector of vector, each series of point being a vector of double.
	 *
	 * @tparam Xkind the X axis is parametrised by \p Xkind (string, int, double...)
	 * The Y axis is always represented by double.
	 * @todo replace Xkind by std::string everywhere.
	 */
	class PlotData ;

	/*! @brief The graph (2D).
	 * This class joins a PlotStyle and a PlotData to build up a graph.  A
	 * filename should be provided as well, indicating where the output
	 * graph and scripts will be generated.
	 *
	 * @warning the filename will get a random suffix before the extension
	 * so as not to overwrite files "par inadvertance".
	 */
	class PlotGraph ;

} // LinBox


/* ********************** */
/*        Tags            */
/* ********************** */

// Tags (cf linbox/linbox-tags.h )
namespace LinBox { namespace Tag {

//! selection of best times in a series.
#if HAVE_CXX11
	enum struct TimeSelect : int32_t
#else
	struct TimeSelect { enum enum_t
#endif
		{
			average   = 1, //!< select the average
			bestThree = 2, //!< select the average among the three best (supposes order is <)
			bestOne   = 3, //!< select the best one (supposes order is <)
			median    = 4, //!< median
			medmean   = 5  //!< interquartile mean (remove 25% extremes and average the rest)
		};
#if HAVE_CXX11
#else
	};
#endif

} // Tag
} // LinBox



/* ********************** */
/*        Outils          */
/* ********************** */


namespace LinBox {

	//! vector of double
	typedef std::vector<double>       dvector_t;
	typedef std::vector<std::string>  svector_t;

	//! matrix of double
	typedef std::vector<dvector_t>    dmatrix_t;

	/** @brief this structure holds a bunch of timings.
	 * It collects the points, the time spent at each point and a measure
	 * (for instance mflops).
	 * @todo Times and Values could be dmatrix_t (and mergeable)
	 */
	struct DataSeries {
		svector_t   PointLabels ; //!< points abscisa, values for the x axis. Used in legend for the X axis.
		dvector_t   Points      ; //!< points abscisa, values for the x axis. Used in TimeWatcher (for instance, if PointLabels are the names of sparse matrices, Points would be their number of non zeros, or 1,2,3,... or whatever relevant for predicting time)
		dvector_t   Times       ; //!< actual computation times.
		dvector_t   Values      ; //!< actual data to be plotted (for instance mflops)

		//! Constructor
		DataSeries() ;

		~DataSeries() ;

		/** @brief resize
		 * @param n new size
		 * @pre the size before was n-1
		 */
		void resize(const index_t & n);

		//! Size of the series of measurements.
		index_t size() const;

		void push_back(const std::string & nam, const double & val, const double & x = NAN, const double &y = NAN)
		{
			linbox_check(PointLabels.size() == Values.size());

			PointLabels.push_back(nam);

			Values.push_back(val);

			if ( std::isnan(x) )
				Points.push_back((double)Points.size());
			else
				Points.push_back(x);

			if ( std::isnan(y))
				Times.push_back(val);
			else
				Times.push_back(y);

			return;
		}


	}; // DataSeries

	/*! Helper.
	 * This helper has several functions :
	 *   - Records the timings
	 *   - predict the execution time for the next experiment
	 *   - helps producing enough experiments (but not too much and not too time consuming) for producing a valid measure.
	 *   .
	 *   @warning if the timings are too short, this may not be accurate.
	 *
	 *   See member function help for more information.
	 */
	class TimeWatcher  {
	private :
		dvector_t  *    Points_; //!< Points data. If <code>Points_[i] = x</code>, then <code>Values_[i]=f(x)<code>.
		dvector_t  *    Values_; //!< Time data. See \p Points_ .

		index_t       MaxRepet_; //!< Maximum number of repetitions of timings
		index_t       MinRepet_; //!< Minimum number of repetitions of timings
		double         MaxTime_; //!< Maximum time to be spent on repetitions (after MinRepet_ iters have been done)
		double       AbortTime_; //!< Time to abort a series of computation.
		bool           aborted_; //!< abort any subsequent computation


	public:
		/** Constructor.
		 * @param size number of experiments expected
		 * @param Line current line in PlotData
		 */

		TimeWatcher (dvector_t & pts, dvector_t & vals) ;
		TimeWatcher () ;

		void init(dvector_t & pts, dvector_t & vals);

		//! returns the vector of abscissa (points)
		dvector_t & refX() ;

		//! returns the vector of ordiantes (values)
		dvector_t & refY() ;

		/** Prediction for the next experiment time.
		 * It is assumed that \c predict(0)=0. If Curent_<3, a linear,
		 * then quadratic fit is done. Other wise, a cubic fit is
		 * performed.
		 * @param x the next evaluation point.
		 * @return f(x) where f tries to fit the points : \f$ f(\mathtt{Data\_}[0][0..\mathtt{Current\_}-1]) \approx  refY()[0..\mathtt{Current\_}-1]\f$
		 */
		double predict(double x) ;

		/*! @brief Watches a timer and a number and repet and signals if over.
		 *
		 * We want at least 2 repetions but not more than maxtime spent on timing.
		 *
		 * @param repet number of previous repetitions. Should be 0 on the first time \c keepon is called.
		 * @param tim timer to watch
		 * @param maxtime maximum time (in seconds) until \c watchon tells stop.
		 * @return \c true if we conditions are not met to stop, \c false otherwise.
		 * @pre \c tim should have been started previously !
		 *
		 */
		bool keepon(index_t & repet, double tim, bool usePrediction = false) ;

		//! size
		index_t size() const
		{
			if (Points_ == NULL || Values_ == NULL) {
				linbox_check(Values_ == NULL && Points_ == NULL);
				return  0 ;
			}
			linbox_check(Points_->size() == Values_->size());
			return (index_t)Points_->size();
		}
	} ; // TimeWatcher


	template<class MyTimer>
	class Chrono  {
	private :
		MyTimer   _chrono_ ;
		dvector_t _times_  ;
		double    _total_  ;
	public:
		Chrono()  :
			_chrono_(),
			_times_(0),
			_total_(0)
		{
			_chrono_.clear();
		}

		~Chrono() {} ;

		void clear()
		{
			_chrono_.clear();
			_times_.resize(0);
			_total_ = 0;
		}

		void start()
		{
			_chrono_.start();
		}

		void stop()
		{
			_chrono_.stop();
			double split = _chrono_.userElapsedTime();
			_times_.push_back( split );
			_total_ += split ;
			_chrono_.clear();
		}

		double time() const
		{
			return _total_ ;
		}

		dvector_t times() const
		{
			return _times_ ;
		}


	};

	bool isDigit (const std::string & s)
	{
		std::istringstream ss(s);
		double d = 0.0;
		ss >> d ; // try to read.
		ss >> std::ws;  // suppress whitespace

		return (!ss.fail() && ss.eof()) ;
	}

	bool fortifiedString(const std::string & s)
	{
		if (isDigit(s))
			return true ;
		linbox_check(!s.empty());
		return s.front() == '\"' && s.back() ==  '\"' ;
	}

	std::string unfortifyString(const std::string &s)
	{
		std::string t = s ;
		if (fortifiedString(s)) {
			t.erase(t.begin());
			t.pop_back();
		}
		return t;
	}
	std::string fortifyString(const std::string & s)
	{
		if (fortifiedString(s))
			return s ;
		string r = "\"" ;
		return r + s + "\"";
	}

	template<class T>
	std::string toString(T & nam)
	{
		std::ostringstream nam_ss ;
		nam_ss << nam ;
		return nam_ss.str();
	}

}// LinBox

// advancement printing on terminal
namespace LinBox {

	/*! show the advancement (on the terminal)
	 * suppose linear advancement
	 * @param curr current iteration
	 * @param min starting iteration
	 * @param max terminal iteration
	 */
	void showAdvanceLinear(index_t curr, index_t min, index_t max);

	/*! tells the current series of measure has completed (on the terminal)
	 * @param curr current iteration
	 * @param all number of iterations
	 */
	void showFinish(index_t curr, index_t all);

	/*! tells the current series of measure was skipped (on the terminal)
	 * @param curr current iteration
	 * @param all number of iterations
	 */
	void showSkip(index_t curr, index_t all);


	class showProgression {
	private :
		index_t _cur_ ;
		index_t _tot_ ;
	public :
		showProgression (index_t tot) :
			_cur_(0)
			,_tot_(tot)
		{}

		void FinishIter()
		{
			++_cur_ ;
			showFinish(_cur_,_tot_);
		}

		void SkipIter()
		{
			++_cur_ ;
			showSkip(_cur_,_tot_);
		}
	};
} // LinBox

// computing megaflops from time (series)/number ops
namespace LinBox {

	/**
	 * @brief computes the number of megaflops.
	 * @param tim timer (seconds)
	 * @param mflo number of operations (1e6 operations)
	 * @param rpt number of experiences
	 * @return mflo/(tim*rpt)
	 */
	double computeMFLOPS(const double & tim, const double mflo, const index_t rpt = 1);

	/*! @internal @brief inserts a time in a vector of 3 best times.
	 * @param tim3 ordered form min to max vector of 3 best times
	 * @param tps inserts that time in \p tim3 if better.
	 * @return a reference to \p tim3
	 */
	dvector_t &
	insertTime(dvector_t & tim3, const double & tps);

	/**
	 * @brief computes the number of megaflops.
	 * @param tim timer (seconds)
	 * @param mflo number of operations (1e6 operations)
	 * @param ts number of experiences to select. @see TimeSelect. Default to the best three
	 * @return mflo/(tim*rpt)
	 */
	double computeMFLOPS(const dvector_t & tim, const double mflo, LINBOX_enum(Tag::TimeSelect) ts = Tag::TimeSelect::bestThree);

}

#include "benchmarks/benchmark.inl"

#ifdef LinBoxSrcOnly
#include "benchmarks/benchmark.C"
#endif

#endif // __LINBOX_benchmarks_benchmark_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
