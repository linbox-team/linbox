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

#include <sys/utsname.h>
#include <ctime>
#ifdef __LINBOX_HAVE_TINYXML2
#include <tinyxml2.h>
#endif


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

	//! selection of printers
#if HAVE_CXX11
	enum struct Printer : int32_t
#else
			      struct Printer { enum enum_t
#endif
				      {
					      dat       = 1, //!< print data  raw (format chosen by implementation)
					      tex       = 2, //!< print latex useable
					      xml       = 3, //!< print using xml syntax
					      gnuplot   = 4, //!< print using gnuplot
					      csv       = 5, //!< coma separated values
					      html      = 6  //!< coma separated values
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

typedef uint32_t index_t ;

//
// typedefs data formats
//

namespace LinBox {

	//! vector of double
	typedef std::vector<double>       dvector_t;
	typedef std::vector<std::string>  svector_t;

	//! matrix of double
	typedef std::vector<dvector_t>    dmatrix_t;
	typedef std::vector<svector_t>    smatrix_t;

} // LinBox

//
// TimeWatcher
//

namespace LinBox {
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
		/*! constructor.
		 * Inits the time watcher with a pair of points/values
		 * @param pts vector of points
		 * @param vals vector of times
		 */
		TimeWatcher (dvector_t & pts, dvector_t & vals) ;

		//! Null Constructor. The pointers are intialised to NULL
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
		index_t size() const;

		//! clear the pointers (not the settings)
		void clear();

	} ; // TimeWatcher
} // LinBox


//
// Chrono
//

namespace LinBox {
	template<class MyTimer>
	class Chrono  {
	private :
		MyTimer   _chrono_ ;
		dvector_t _times_  ;
		double    _total_  ;
	public:
		Chrono()  ;

		~Chrono()  ;

		void clear();

		void start();

		void stop();

		double time() const;

		dvector_t times() const;



	}; // Chrono
} // LinBox

//
// DataSeries
//

namespace LinBox {
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

		//! add some new data.
		void push_back(const std::string & nam, const double & val, const double & x = NAN, const double &y = NAN);

	}; // DataSeries

} // LinBox

//
//String processing
//

namespace LinBox {

	/** Check if a string is actually a double.
	 * @param s string to check
	 * @return true/false
	 */
	bool isDigit (const std::string & s);

	/** Tells is a string has double quotes around.
	 * @param s string to test
	 * @return true if s[0] == s[last] == '"'
	 */
	bool fortifiedString(const std::string & s);

	/** removes the surrounding quotes.
	 * @param s removes quotes around if necessary
	 * @return s without surrounding double quotes if necessary
	 */
	std::string unfortifyString(const std::string &s);

	/** adds surrounding quotes.
	 * @param s add quotes around if necessary
	 * @return s with surrounding double quotes if necessary
	 */
	std::string fortifyString(const std::string & s);

	/** Converts anything to a string
	 * @param nam to be put in a string.
	 */
	template<class T>
	std::string toString(T & nam)
	{
		std::ostringstream nam_ss ;
		nam_ss << nam ;
		return nam_ss.str();
	}

}// LinBox

/* ********************** */
/*    Plot structures     */
/* ********************** */

//
// PlotStyle
//

namespace LinBox {

	/*! @brief Represents a table of values to plot (2D).
	 * list of values are reprensented by vectors.  the table is a vector
	 * of these vectors.
	 *
	 * @warning NaN, inf are used as missing data. More genenally
	 * we could store data in strings.
	 * @todo Allow for 'speed up against col X' style
	 * @todo make depend on PlotData
	 */
	//!@todo setUsingSeries(const svector_t &)
	class PlotStyle {
	public:
		//! What format the plot should be in?
		struct Term {
			//! Term type
			enum Type {
				png  = 100, //!< png. Portable Network Graphics file.
				pdf  = 101, //!< pdf. Portable Document Format actually, this is eps converted to pdf.
				eps  = 102, //!< eps. Encapsulated PostScript. Cool for inclusion in LaTex files. This is the default.
				svg  = 103, //!< sgv. Scalable Vector Graphics.
				tex  = 104, //!< tex. Simple tabular.
				html = 105, //!< html. HTML table.
				other= 106 //!< don't know yet...
			} ;
		};

		// enum NoType { other = 0 } ;
		//! What style of graphic : histogram ? graph ?
		struct Plot {
			//! Plot type
			enum Type {
				histo = 200, //! histogram plot. This is the default. x ticks are evenly spaced, whatever there value and are labelled with their value.
				graph = 201, //! standard plot. Plots y_i=f(x) with x in the first colum and y_i in ith column. x-ticks are well spaced. This will not work if the X are not numbers (but strings).
				other = 202  //! other (ie user supplied).
			} ;
		};

		struct Line {
			enum Type {
				lines      = 300,  //! lines.
				histogram  = 301,  //! histogram (boxes).
				linespoints= 302,  //! lines with points. (default)
				points     = 303,  //! only the points.
				other      = 304   //! rien.
			} ;
		};

		struct Options {
			enum Type {
				oblique = 400,
				other   = 401
			};
		};

		/*! @brief Constructor.
		 * By default, creates an histogram representing the data in an eps plot.
		 */
		PlotStyle() ;

		/*! @brief sets the titles in the graph.
		 * @param titre Title of the graph
		 * @param titre_y Title of the y-axis (series)
		 * @param titre_x Title of the x-axis (data points)
		 */
		void setTitle ( const std::string  &  titre
				, const std::string  & titre_y
				, const std::string  & titre_x);

		/*! @brief Gets the title of the graph.
		 * @return a gnuplot command to set the title of the graph.
		 */
		std::string  getTitle() const ;

		/*! @brief Gets the title of points abscissa.
		 * @return a gnuplot command to set the title of the abscissa.
		 */
		std::string getTitleX() const ;

		/*! @brief Gets the title of the series.
		 * @return a gnuplot command to set the title of the ordinate (series).
		 */
		std::string getTitleY() const ;

		/*! @brief get the title string.
		 * @param index can be (0,1,2)
		 */
		std::string getRawTitle(int index=0) const ;

		/*! @brief Sets the output format.
		 * @sa TermType
		 * @param term type
		 */
		void setTerm( enum Term::Type term) ;


		/*! @brief Gets the output format.
		 * @return string for setting the expected output format in gnuplot.
		 * @warning noenhanced allows underscores while enhanced does subscripts.
		 * if we add a (no) enhanced option, we'll have to add a safeFormat(std::string) that replaces \c _ with <code>\_</code> .
		 * This is tricky and can be done at "post production" stage :-)
		 */
		std::string getTerm() const ;

		/*! @brief Gets the graph output extension.
		 * By default, this is ".eps".
		 * @return a string for this extension, including the sepatating dot.
		 *
		 */
		std::string getExt() const ;

#if 0
		/*! @brief gets the style of the graph.
		 * This is very user-tweakable !!
		 * @return the style for gnuplot.
		 */
		const std::string & getStyle()
		{
			return "#style\n"+_styleopts_ ;
		}

		/*! @brief sets the style of the graph.
		 * This is very user-tweakable !!
		 * @param style the style for gnuplot as a gnuplot command.
		 */
		void setStyle(const std::string & style)
		{
			_styleopts_ = style ;
		}
#endif

		/*! @brief sets the legend position.
		 * @param keypos the arguments to key (where the legend should be put)
		 * can be :
		 * <code>
		 *      set key {on|off} {default}
		 *              {{inside | outside} | {lmargin | rmargin | tmargin | bmargin} | {at <position>}}
		 *              {left | right | center} {top | bottom | center}
		 *              {vertical | horizontal} {Left | Right}
		 *              {{no}reverse} {{no}invert}
		 *              {samplen <sample_length>} {spacing <vertical_spacing>}
		 *              {width <width_increment>}
		 *              {height <height_increment>}
		 *              {{no}autotitle {columnheader}}
		 *              {title "<text>"} {{no}enhanced}
		 *              {{no}box { {linestyle | ls <line_style>} | {linetype | lt <line_type>} {linewidth | lw <line_width>}}}
		 * </code>
		 */
		void setKeyPos(const std::string & keypos) ;


		/*! @brief Gets the legend position.
		 * by default, it is "under".
		 */
		std::string getKeyPos() const ;

		/*! @brief sets the  position of the labels on the X absciss.
		 * @param opt
		 * @param more more stuff
		 */
		void setXtics ( enum Options::Type opt, const std::string & more="") ;

		/*! @brief Gets the legend position.
		 * by default, it is 45° inclined (use in on long tics legends).
		 */
		const std::string & getXtics() const ;

		/*! @brief Gets the name of the output graph.
		 * @param basnam the raw name for the output.
		 * @return basnam+extenstion.
		 */
		std::string getOutput(const std::string & basnam) const ;

		/*! @brief Sets the type of plot.
		 * @param type the type.
		 * @sa PlotType
		 *
		 */
		void setPlotType(enum Plot::Type type) ;

		/*! @brief Sets the way dots are linked.
		 * @sa LineType
		 * @param type type
		 */
		void setLineType( enum Line::Type type) ;

		/*! @brief Gets the type of plot.
		 * default is histogram, or if graph is supplied, then the default is linespoints.
		 * Can be totally customized.
		 * @return a string for gnuplot to set the plot type.
		 * @sa PlotType
		 *
		 */
		std::string getPlotType(const std::string & extraargs ="") ;// const

		/*! @brief adds some style line to the graph.
		 * This is very user-tweakable !!
		 * @param style a style line for gnuplot as a gnuplot command.
		 */
		void addPlotType(const std::string & style) ;

		/*! @brief tells which columns to use.
		 * @param col a column to use.
		 * @param moreargs more stuff
		 */
		void setUsingSeries(index_t col, const std::string &  moreargs= "") ;

		/*! @brief adds a column to use
		 * @param col a  column to use.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(index_t col, const std::string &  moreargs= "") ;

		/*! @brief tells which columns to use.
		 * @param cols a list of column to use.
		 * @param moreargs more stuff
		 */
		void setUsingSeries(std::list<index_t> cols, const std::string & moreargs= "") ;

		/*! @brief adds a set of columns to use.
		 * @param cols a list of column to use.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(std::list<index_t> cols, const std::string & moreargs= "") ;

		/*! @brief tells which columns to use.
		 * @param cols all colums between \c cols.first and \c cols.second (included)
		 * will be used.
		 * @param moreargs more stuff
		 *
		 */
		void setUsingSeries(std::pair<index_t,index_t> cols, const std::string & moreargs= "") ;

		/*! @brief adds contiguous columns to use.
		 * @param cols all colums between \c cols.first and \c
		 * cols.second will be used.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has
		 * already been called.
		 *
		 */
		void addUsingSeries(std::pair<index_t,index_t> cols, const std::string & moreargs= "") ;

		const std::string & getUsingSeries() const ;

	private :
		// int                                 _precision_ ;   //!< precision of the output. by default 2.
		/* Legend. */
		std::string                         _legend_pos_;   //!< legend position
		/*  titles */
		std::string                         _title_     ;   //!< name of the graph
		std::string                         _title_x_   ;   //!< title for the points
		std::string                         _title_y_   ;   //!< title for the series
		std::string                         _xtics_     ;   //!< format for the x tics.
		/*  units */
		// std::string                      _unit_      ;
		/*  terminal output */
		enum Term::Type                     _term_      ; //!< output data format.
		// std::string                      _termopts_  ;
		/*  plotting style */
		enum Plot::Type                     _plot_type_ ; //!< histogram/graph style
		// std::string                         _plot_extra_; //!< extra specification for the plot style. default empty.
		enum Line::Type                     _line_type_ ; //!< style for the representation of points
		std::string                         _styleopts_ ; //!< gp style command.
		/*  columns to use */
		std::string                         _usingcols_ ; //!< columns to be used (gp command)


	} ; // PlotStyle

} // LinBox


namespace LinBox {
	/*! @brief The raw data to plot.
	 * Represents  the labels for the points (X axis) and the values for
	 * each series of measures (Y axis).
	 *
	 * @internal The internal representation is a
	 * vector of vector, each series of point being a vector of double.
	 *
	 * @tparam Xkind the X axis is parametrised by \p Xkind (string, int, double...)
	 * The Y axis is always represented by double.
	 * @todo put the legend (title, x, y) in there
	 */
	class PlotData {
	private :
		std::vector<DataSeries > _tableau_     ;   //!< data. \c _tableau_[i] represents a series of measurements. A data series is resized only when a new element comes in.
		std::vector< std::string >      _serie_label_ ;   //!< label for each serie of measures. Used in the legend of the plots/tables of points.
		index_t                         _curr_serie_  ;   //!< index of the current series of measurements.
		TimeWatcher                     _time_watch_  ;   //!< time predictor, helper. See \c TimeWatcher.
	private:

#ifdef __LINBOX_HAVE_TINYXML2
		tinyxml2::XMLElement * saveData(tinyxml2::XMLDocument & doc) ;
#endif
	public :

		/*! Inits a plot with series of data.
		 * @param nb_pts number of points in each serie.
		 * @param nb_srs number of series of points. Default is 1.
		 */
		PlotData() ;
		/*! destructor.
		*/
		~PlotData()  ;

		/*! copy constructor.
		 * @param PD a PlotData to copy.
		 */
		PlotData(const PlotData & PD);

		/** Return the index of a series according to its name.
		 * Creates one if necessary.
		 * the current serie is set to this one.
		 * @param name of the series
		 * @return index of the series
		 */
		index_t getIndex(const std::string & nom) ;

		/** @brief initialize to empty
		*/
		void clear() ;

		void merge(const PlotData &PD) ;

		/*! @brief get the number of series.
		 * @return number of series.
		 */
		index_t size() const ;

		/** @brief gets the current series number.
		*/
		index_t getCurrentSerieNumber() const ;

		/*! @brief Sets the name of a serie.
		 * @param i index of the serie
		 * @param nom name of the serie
		 */
		void setSerieName(const index_t & i, const std::string & nom) ;

		/** Gets the name of a serie.
		 * @param i index of the serie
		 */
		const std::string &  getSerieName(const index_t & i) const ;

		/** Gets the name of the current serie.
		*/
		const std::string &  getCurrentSerieName() const ;

		/*! @brief Sets the name of the current serie.
		 * @param nom name of the serie
		 */
		void setCurrentSerieName(const std::string & nom) ;


		/** Inits the watch on a series
		 * @param i index of a series
		 */
		void initWatch ( const size_t & i) ;

		/** Inits the watch to current series
		*/
		void initCurrentSeriesWatch () ;

		/** Creates a new series.
		 * It is created after the last series.
		 * \c getCurrentSeries() points to it.
		 * @param nom name of the new series
		 */
		void newSerie(const std::string & nom = "") ;

		/** Finish a serie of measurements.
		 * Nothing is done for the moment.
		 */
		void finishSerie() ;

		/*! Returns the ith series of measurements.
		 * @param i ith series to be returned
		 */
		const DataSeries & getSeries(const index_t  &i) const ;

		/*! Returns the ith series of measurements.
		 * @param i ith series to be returned
		 */
		DataSeries & refSeries(const index_t  &i) ;


		/*! Returns the current series of measurements.
		*/
		const DataSeries & getCurrentSeries() const ;

		/*! Returns the current series of measurements.
		*/
		DataSeries & refCurrentSeries() ;

		/** @brief size of a series.
		 * @param i index of the series
		 */
		index_t getSerieSize(const index_t & i) const ;

		/** @brief size of the current series.
		*/
		index_t getCurrentSerieSize() const ;

		/** Sets the current series of measurements.
		 * @param s index for the series we wish to add measures to.
		 */
		void selectSeries(const index_t & s) ;

		/*! goes to the next series of points
		*/
		bool selectNextSeries() ;

		/*! @brief Sets the name of a point.
		 * @param i series number
		 * @param j index for the the point
		 * @param nom name of the point
		 */
		template<class T>
		void setSeriesPointLabel(const index_t & i, const index_t & j, const T & nom)
	{
		std::string nom_s = fortifyString(toString(nom));
		linbox_check(j<getSerieSize(i) );
		refSeries(i).PointLabels[j] = nom_s ;

		return;
	}
		/*! @brief Sets the name of a point.
		 * @param j index for the the point
		 * @param nom name of the point
		 */
		void setCurrentSeriesPointLabel(const index_t & j, const std::string & nom) ;

		/*! @brief gets the name of a point.
		 * @param i series number
		 * @param j its index.
		 * @return its name.
		 * @warning no default. \c setPointXLabel has to be used beforehands.
		 */
		const std::string & getSeriesPointLabel(const index_t &i, const index_t & j) const ;

		/*! @brief gets the name of a point.
		 * @param j its index.
		 * @return its name.
		 * @warning no default. \c setPointXLabel has to be used beforehands.
		 */
		const std::string & getCurrentSeriesPointLabel(const index_t & j) const ;

		/*! @brief gets the name of a serie.
		 * Defaults to \c "serie.i"
		 * @param i its index.
		 * @return its name.
		 */
		std::string getSerieLabel(const index_t & i) const ;

		/*! @brief gets the name of a serie.
		 * Defaults to \c "serie.i"
		 * @param i its index.
		 * @return its name.
		 */
		std::string  getCurrentSerieLabel() const;

		/*! @brief gets all the names in the series.
		 * @return a vector of names.
		 */
		const std::vector<std::string > & getSerieLabels() const ;

		/*! @brief gets all the names in the points of a series
		 * @param  i index of the series
		 * @return a vector of names.
		 */
		const svector_t & getSeriePointLabel( const index_t & i) const ;

		/*! @brief gets all the names in the points.
		 * @return a vector of names.
		 */
		const svector_t & getCurrentSeriesPointLabel() const ;

		/*! @brief gets all the names in the points of a series.
		 * @param i index of the series
		 * @return a vector of names.
		 */
		const dvector_t & getSeriesValues(const index_t & i) const ;

		/*! @brief gets all the names in the points.
		 * @return a vector of names.
		 */
		const dvector_t & getCurrentSeriesValues() const ;

		/*! @brief sets a new entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 * @param xval x value of the point (eg size of the matrix, of a sparse matrix,...)
		 * @param yval time for this computation (seconds)
		 */
		void setSeriesEntry(const index_t &i, const std::string & nam, const double & val
				    , const double & xval = NAN, const double & yval = NAN) ;

		/*! @brief sets a new entry.
		 * @param name name of the series
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 * @param xval x value of the point (eg size of the matrix, of a sparse matrix,...)
		 * @param yval time for this computation (seconds)
		 */
		void setEntry(const std::string & nom, const std::string & nam, const double & val
			      , const double & xval = NAN, const double & yval = NAN) ;


		/*! @brief sets a new entry.
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 */
		template<class T>
		void setCurrentSeriesEntry(const T & nam, const double & val
					   , const double & xval = NAN, const double & yval = NAN)
	{
		std::string nam_s = fortifyString(toString(nam));
		return setSeriesEntry(_curr_serie_,nam_s,val,xval,yval) ;
	}

		/*! @brief gets a value for an entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @return val value of point j in ith serie.
		 */
		double getSeriesEntry(const index_t & i, const index_t & j) const ;

		/*! @brief gets a value for an entry.
		 * @param j index of the point
		 * @return val value of point j in current serie.
		 */
		double getCurrentSeriesEntry(const index_t & j) const ;


		/*! @brief gets a time spent on an entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getSeriesEntryTime(const index_t &i, const index_t & j) const ;

		/*! @brief gets a time spent on an entry.
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getCurrentSeriesEntryTime(const index_t & j) const ;


		/*! @brief gets the point corresponding to an entry.
		 * @warning, this is not the label, but the value associated to the point
		 * @param i index of the series
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getSeriesEntryPoint(const index_t & i, const index_t & j) const ;

		/*! @brief gets the point corresponding to an entry.
		 * @warning, this is not the label, but the value associated to the point
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getCurrentSeriesEntryPoint(const index_t & j) const ;


		/*! gets a reference to the array of data.
		 * @return a reference to the member \c _tableau_ representing the data.
		 */
		const std::vector<DataSeries > & getTable() const ;

		/*! gets a reference to the array of data.
		 * @return a reference to the member \c _tableau_ representing the data.
		 */
		std::vector<DataSeries > & refTable() ;


		/** @brief Continue for another time measure ?
		 * @see TimeWatcher::keepon
		 * @param repet  current number of repetitions for this new measure
		 * @param tim    time previously spent on the measures.
		 * @return true if one more measure can be done
		 */
		bool keepon(index_t & repet, double tim, bool usePrediction=false) ;

		void load( const std::string & filename) ;

		/** @brief saves the data in XML format.
		 * @param filename file name
		 * @param title titles of the data
		 * @param xtitle legend of the X axis
		 * @param ytitle legend of the Y axis.
		 */
		void save( const std::string & filename
			   , const std::string & title = ""
			   , const std::string & xtitle = ""
			   , const std::string & ytitle = "") ;


	}; // PlotData

} // LinBox

//
// PlotGraph
//

namespace LinBox {

	/*! @brief The graph (2D).
	 * This class joins a PlotStyle and a PlotData to build up a graph.  A
	 * filename should be provided as well, indicating where the output
	 * graph and scripts will be generated.
	 *
	 * @warning the filename will get a random suffix before the extension
	 * so as not to overwrite files "par inadvertance".
	 * @warning don't name anything else than a folder "data" in your working directory. You've been warned.
	 * @todo make depend on PlotStyle (that owns data)
	 */
	//!@todo use getUsingSeries in latex/html/csv/xml
	class PlotGraph {
	private :
		PlotData              & _data_ ;   //!< reference to the data points
		PlotStyle            & _style_ ;   //!< reference to a plotting style
		std::string         _filename_ ;   //!< name for the output file (without extension). a random \c _XXXXXX suffix will be added to make it unique.
		std::string        _printname_ ;   //!< name for the output file (without extension) to be printed. a random \c _XXXXXX suffix makes it unique.
		dmatrix_t         _merge_data_ ;
		svector_t       _merge_points_ ;

	private :
		/*! @internal
		 * @brief random <code>:alnum:</code> \c char.
		 * [[:alnum:]] characters are in range
		 *  - num : 48-57
		 *  - AL  : 67-90
		 *  - al  : 97-122
		 *  .
		 * @return a random alphabetic or numeric char.
		 */
		char _randomAlNum();

		/*! @internal
		 * @brief Appends random suffix.
		 * Appends to \p _filename_ a random string constituted of an
		 * underscore followed by 8 random alnum chars.
		 * @return the concatenation of \c _filename_ and this suffix.
		 */
		void _randomName();

		const std::string & getFileName()  ;

		//! @bug this supposes the two series have unique measurements for one point.
		void mergeTwoSeries( svector_t & merge_points
				     , dmatrix_t & merge_data
				     , const svector_t & pts
				     , const dvector_t & dat
				     , const index_t & idx) const ;

		//! merge all series of points into a vector of absissa points  and a vector of vector of data points
		void mergeSeries();

		void print_csv();

		void print_dat();

		void print_xml();

		void print_html() ;

		/*! @brief Prints data in a latex tabular.
		*/
		void print_latex();

		/*!@brief Plots the data with gnuplot.
		 * Produces data in a .dat file, creates a .gp gnuplot script and
		 * outputs a graph calling gnuplot.
		 * @warning If gnuplot is not available, fall back to the latex method.
		 */
		void print_gnuplot(bool only_data=false);

	public :

		/*! @brief Sets a new data structure.
		 * @param data a reference to a PlotData class.
		 */
		void setData( PlotData & data );

		/*! @brief Gets the data.
		 * @param[in,out] data a reference to a PlotData class.
		 */
		PlotData & refData( PlotData & data);

		/*! @brief Sets a new style structure.
		 * @param style a reference to a PlotStyle class.
		 */
		void setStyle( PlotStyle & style ) ;

		/*! @brief Gets the style.
		 * @param[in,out] style a reference to a PlotStyle class.
		 */
		PlotStyle & refStyle( PlotStyle & style) ;


		// not implemented yet
		void sortSeries()  ;

		// not implemented yet
		void unique() ;

		/*! @brief Constructor for the PlotGraph class.
		 * Plots a series of data according to a style.
		 * @param data data to be plot, will be processed by the style
		 * @param style sets parameters to gnuplot to achieve a nice
		 * plot.
		 */
		PlotGraph( PlotData & data, PlotStyle & style );

		/*! @brief sets the ouput file name.
		 * All output is put in a "data" subfolder.
		 * @warning Since no file is overwritten, this
		 * directory can rapidly get very populated.
		 */
		void setOutFilename( const std::string & filename ) ;

		const std::string & getUsingSeries() ; // const

		/*! @brief Gets the plot command line.
		 * @param File the name of/path to the data file (with extension)
		 * @return a gnuplot "plot" command stream.
		 */
		std::string getPlotCommand(const std::string & File) ; //const

		void print( LINBOX_enum(Tag::Printer) pt = Tag::Printer::xml) ;

		void save() ;

		void load(const std::string & filename) ;

	}; // PlotGraph

} // LinBox



//
// advancement printing on terminal
//

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

	//! Show progression on the terminal (helper)
	class showProgression {
	private :
		index_t _cur_ ; //!< current iter
		index_t _tot_ ; //!< max iter
	public :
		//! constructor
		showProgression (index_t tot) ;

		//! show an inter has finished
		void FinishIter();

		//! show an inter has been skipped.
		void SkipIter();
	}; // showProgression
} // LinBox

//
// computing megaflops from time (series)/number ops
//

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

//
// Machine info
//


namespace LinBox {

	//! get ISO time and date
	std::string getDateTime();

	//! get some machine information (not cpu yet)
	smatrix_t getMachineInformation();

} // LinBox

//
// Least Squares
//
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

//
// Curve fitting
//

namespace LinBox {

	//! fit X[n-1,n],Y[n-1,n] and return evaluation at x.
	double fit2(const dvector_t & X, const dvector_t & Y, int n, double x);

#ifdef __LINBOX_HAVE_LAPACK
	//! fits with a degree 3 polynomial and return evaluation at x.
	double fit_lapack3(const dvector_t &X, const dvector_t &Z, double x);
#endif // __LINBOX_HAVE_LAPACK


	//! fit X[n-2,n],Y[n-2,n] and return evaluation at x.
	double fit3(const dvector_t & X, const dvector_t & Y,int n, double x);

} // LinBox




#ifdef LinBoxSrcOnly
#include "benchmarks/benchmark.C"
#endif

#include "benchmarks/benchmark.inl"

#endif // __LINBOX_benchmarks_benchmark_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
