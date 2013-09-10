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

/*! @file   benchmarks/benchmark.inl
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

#ifndef __LINBOX_benchmarks_benchmark_INL
#define __LINBOX_benchmarks_benchmark_INL

#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include "tests/test-common.h"
#include "linbox/util/timer.h"
#include <string>
#include <fstream>


// Plot structures (Style, Data, Graph)
namespace LinBox {


	// template by string *always*
	template<class Xkind>
	class PlotData {
	private :
		std::vector<DataSeries<Xkind> > _tableau_     ;   //!< data. \c _tableau_[i] represents a series of measurements. A data series is resized only when a new element comes in.
		std::vector< std::string >      _serie_label_ ;   //!< label for each serie of measures. Used in the legend of the plots/tables of points.
		TimeWatcher                     _time_watch_  ;   //!< time predictor, helper. See \c TimeWatcher.
		index_t                         _curr_serie_  ;   //!< index of the current series of measurements.
	private :
		bool fortifiedString(const std::string & s)
		{
			linbox_check(!s.empty());
			return s.front() == '\"' && s.back() ==  '\"' ;
		}

		std::string fortifyString(const std::string & s)
		{
			if (fortifiedString(s))
				return s ;
			string r = "\"" ;
			return r + s + "\"";
		}

	public :

		/*! Inits a plot with series of data.
		 * @param nb_pts number of points in each serie.
		 * @param nb_srs number of series of points. Default is 1.
		 */
		PlotData() :
			_tableau_      (0)
			,_serie_label_ (0)
			,_time_watch_  ()
			,_curr_serie_()
		{
		}

		/*! destructor.
		 */
		~PlotData() {} ;

		/*! copy constructor.
		 * @param PD a PlotData to copy.
		 */
		PlotData(const PlotData<Xkind> & PD):
			_tableau_(PD.getTable())
			,_serie_label_(PD.getSerieLabels())
			,_curr_serie_(PD.getSerieNumber())
			,_time_watch_  (_tableau_[_curr_serie_].Points,_tableau_[_curr_serie_].Times)
		{
		}

		/** Return the index of a series according to its name.
		 * Creates one if necessary.
		 * the current serie is set to this one.
		 * @param name of the series
		 * @return index of the series
		 */
		index_t getIndex(const std::string & nom)
		{
			std::vector<std::string>:: iterator it ;
			std::string nomf = fortifyString(nom);
			it = std::find(_serie_label_.begin() , _serie_label_.end() , nomf);
			if ( it != _serie_label_.end() ) {
				return (index_t)std::distance(_serie_label_.begin(),it);
			}
			else {
				index_t sz = (index_t)_serie_label_.size() ;
				_serie_label_.push_back(fortifyString(nom));
				_tableau_.resize(sz+1);
				_curr_serie_ = sz ;
				initWatch(sz);
				return sz ;
			}
		}

		void merge(const PlotData<Xkind> &PD)
		{
			for (size_t i = 0 ; i < PD.size() ; ++i) {
				_tableau_.push_back(PD.getSeries(i));
				_serie_label_.push_back(fortifyString(PD.getSerieName(i)));
			}
			return ;
		}

		/*! @brief get the number of series.
		 * @return number of series.
		 */
		index_t size() const
		{
			linbox_check(_tableau_.size() == _serie_label_.size());
			return (index_t)_tableau_.size() ;
		}

		/** @brief gets the current series number.
		 */
		index_t getCurrentSerieNumber() const
		{
			return _curr_serie_ ;
		}

		/*! @brief Sets the name of a serie.
		 * @param i index of the serie
		 * @param nom name of the serie
		 */
		void setSerieName(const index_t & i, const std::string & nom)
		{
			linbox_check(i<size());
			_serie_label_[i] = fortifyString(nom) ;
		}

		/** Gets the name of a serie.
		 * @param i index of the serie
		 */
		std::string  getSerieName(const index_t & i) const
		{
			linbox_check(i<size());
			return _serie_label_[i] ;
		}

		/** Gets the name of the current serie.
		 */
		std::string  getCurrentSerieName() const
		{
			return getSerieName(_curr_serie_) ;
		}

		/*! @brief Sets the name of the current serie.
		 * @param nom name of the serie
		 */
		void setCurrentSerieName(const std::string & nom)
		{
			return setSerieName(_curr_serie_,nom);
		}


		/** Inits the watch on a series
		 * @param i index of a series
		 */
		void initWatch ( const size_t & i)
		{
			linbox_check(i < size());
			_time_watch_.init(_tableau_[i].Points,_tableau_[i].Times);
		}

		/** Inits the watch to current series
		 */
		void initCurrentSeriesWatch ()
		{
			initWatch(_curr_serie_);
		}

		/** Creates a new series.
		 * It is created after the last series.
		 * \c getCurrentSeries() points to it.
		 * @param nom name of the new series
		 */
		void newSerie(const std::string & nom = "")
		{

			index_t old_size = size() ;
			_curr_serie_ = old_size ;

			_tableau_.resize(old_size+1);
			_serie_label_.resize(old_size+1);

			if (nom.empty()) {
				std::ostringstream emptytitle ;
				emptytitle << "serie." << _curr_serie_ ;
				setCurrentSerieName(emptytitle.str());
			}
			else
				setCurrentSerieName(nom);

			initCurrentSeriesWatch();

			return;
		}

		/** Finish a serie of measurements.
		 * Nothing is done for the moment.
		 */
		void finishSerie()
		{
		}

		/*! Returns the ith series of measurements.
		 * @param i ith series to be returned
		 */
		const DataSeries<Xkind> & getSeries(const index_t  &i) const
		{
			linbox_check(i < size());
			return _tableau_[i] ;
		}

		/*! Returns the ith series of measurements.
		 * @param i ith series to be returned
		 */
		DataSeries<Xkind> & refSeries(const index_t  &i)
		{
			linbox_check(i < size());
			return _tableau_[i] ;
		}


		/*! Returns the current series of measurements.
		 */
		const DataSeries<Xkind> & getCurrentSeries() const
		{
			return getSeries(_curr_serie_);
		}

		/*! Returns the current series of measurements.
		 */
		DataSeries<Xkind> & refCurrentSeries()
		{
			return refSeries(_curr_serie_);
		}

		/** @brief size of a series.
		 * @param i index of the series
		 */
		index_t getSerieSize(const index_t & i) const
		{
			return getSeries(i).size();
		}

		/** @brief size of the current series.
		 */
		index_t getCurrentSerieSize() const
		{
			return getSerieSize(_curr_serie_);
		}

		/** Sets the current series of measurements.
		 * @param s index for the series we wish to add measures to.
		 */
		void selectSeries(const index_t & s)
		{
			linbox_check(s < size());
			_curr_serie_ = s ;
			initCurrentSeriesWatch();

			return;
		}

		/*! goes to the next series of points
		 */
		void selectNextSeries()
		{
			++_curr_serie_ ;
			selectSeries(_curr_serie_);
		}

		/*! @brief Sets the name of a point.
		 * @param i series number
		 * @param j index for the the point
		 * @param nom name of the point
		 */
		void setSeriesPointLabel(const index_t & i, const index_t & j, const Xkind & nom)
		{
			linbox_check(j<getSerieSize(i) );
			getSeries(i).PointLabels[j] = nom ;

			return;
		}

		/*! @brief Sets the name of a point.
		 * @param j index for the the point
		 * @param nom name of the point
		 */
		void setCurrentSeriesPointLabel(const index_t & j, const Xkind & nom)
		{
			return setSeriesPointLabel(_curr_serie_,j,nom);
		}

		/*! @brief gets the name of a point.
		 * @param i series number
		 * @param j its index.
		 * @return its name.
		 * @warning no default. \c setPointXLabel has to be used beforehands.
		 */
		Xkind getSeriesPointLabel(const index_t &i, const index_t & j) const
		{
			linbox_check(j<getSerieSize(i));
			return(getSeries(i).PointLabels[j]) ;
		}

		/*! @brief gets the name of a point.
		 * @param j its index.
		 * @return its name.
		 * @warning no default. \c setPointXLabel has to be used beforehands.
		 */
		Xkind getCurrentSeriesPointLabel(const index_t & j) const
		{
			return getSeriesPointLabel(_curr_serie_,j);
		}

		/*! @brief gets the name of a serie.
		 * Defaults to \c "serie.i"
		 * @param i its index.
		 * @return its name.
		 */
		std::string getSerieLabel(const index_t & i) const
		{
			linbox_check(i<size());
			linbox_check(i<_serie_label_.size() );
			if (_serie_label_[i].empty()) {
				std::ostringstream emptytitle ;
				emptytitle << "serie." << i ;
				return emptytitle.str();
			}
			return(_serie_label_[i]);
		}

		/*! @brief gets the name of a serie.
		 * Defaults to \c "serie.i"
		 * @param i its index.
		 * @return its name.
		 */
		std::string getCurrentSerieLabel() const
		{
			return getSerieLabel(_curr_serie_);
		}

		/*! @brief gets all the names in the series.
		 * @return a vector of names.
		 */
		std::vector<std::string > getSerieLabels() const
		{
			return _serie_label_ ;
		}

		/*! @brief gets all the names in the points of a series
		 * @param  i index of the series
		 * @return a vector of names.
		 */
		const svector_t & getSeriePointLabel( const index_t & i) const
		{
			return(getSeries(i).PointLabels) ;
		}

		/*! @brief gets all the names in the points.
		 * @return a vector of names.
		 */
		const svector_t & getCurrentSeriesPointLabel() const
		{
			return getSeriePointLabel(_curr_serie_);
		}

		/*! @brief gets all the names in the points of a series.
		 * @param i index of the series
		 * @return a vector of names.
		 */
		const dvector_t & getSeriesValues(const index_t & i) const
		{
			return(getSeries(i).Values) ;
		}

		/*! @brief gets all the names in the points.
		 * @return a vector of names.
		 */
		const dvector_t & getCurrentSeriesValues() const
		{
			return getSeriesValues(_curr_serie_);
		}

		/*! @brief sets a new entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 * @param xval x value of the point (eg size of the matrix, of a sparse matrix,...)
		 * @param yval time for this computation (seconds)
		 */
		void setSeriesEntry(const index_t &i, const Xkind & nam, const double & val
			      , const double & xval = NAN, const double & yval = NAN)
		{
			refSeries(i).push_back(nam,val,xval,yval);
			// std::cout << "points : " << refSeries(i).Points << std::endl;
			initWatch(i); //  in case series has changed
			return ;
		}

		/*! @brief sets a new entry.
		 * @param name name of the series
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 * @param xval x value of the point (eg size of the matrix, of a sparse matrix,...)
		 * @param yval time for this computation (seconds)
		 */
		void setEntry(const std::string & nom, const Xkind & nam, const double & val
			      , const double & xval = NAN, const double & yval = NAN)
		{
			std::cout << nom << " has index : " << getIndex(nom) << std::endl;
			return setSeriesEntry(getIndex(nom),nam,val,xval,yval);
		}


		/*! @brief sets a new entry.
		 * @param j index of the point
		 * @param nam name of the point (eg size of the matrix, name of a sparse matrix,...)
		 * @param val value to be inserted (eg mflops, sec,...).
		 */
		void setCurrentSeriesEntry(const Xkind & nam, const double & val
					   , const double & xval = NAN, const double & yval = NAN)
		{
			return setSeriesEntry(_curr_serie_,nam,val,xval,yval) ;
		}


		/*! @brief gets a value for an entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @return val value of point j in ith serie.
		 */
		double getSeriesEntry(const index_t & i, const index_t & j) const
		{
			linbox_check(j<getSerieSize());
			linbox_check(i<size());
			return getSeries(i).Values[j] ;
		}

		/*! @brief gets a value for an entry.
		 * @param j index of the point
		 * @return val value of point j in current serie.
		 */
		double getCurrentSeriesEntry(const index_t & j) const
		{
			return getSeries(_curr_serie_,j);
		}


		/*! @brief gets a time spent on an entry.
		 * @param i index of the series
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getSeriesEntryTime(const index_t &i, const index_t & j) const
		{
			linbox_check(j<getSerieSize(i));
			return getSeries(i).Time[j] ;
		}

		/*! @brief gets a time spent on an entry.
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getCurrentSeriesEntryTime(const index_t & j) const
		{
			return getSeriesEntryTime(_curr_serie_,j);
		}


		/*! @brief gets the point corresponding to an entry.
		 * @warning, this is not the label, but the value associated to the point
		 * @param i index of the series
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getSeriesEntryPoint(const index_t & i, const index_t & j) const
		{
			linbox_check(j<getSerieSize(i));
			return getSeries(i).Point[j] ;
		}

		/*! @brief gets the point corresponding to an entry.
		 * @warning, this is not the label, but the value associated to the point
		 * @param j index of the point
		 * @return time for the point j in current serie.
		 */
		double getCurrentSeriesEntryPoint(const index_t & j) const
		{
			return getSeriesEntryPoint(_curr_serie_,j);
		}


		/*! gets a reference to the array of data.
		 * @return a reference to the member \c _tableau_ representing the data.
		 */
		const std::vector<DataSeries<Xkind > > & getTable() const
		{
			return _tableau_ ;
		}

		/*! gets a reference to the array of data.
		 * @return a reference to the member \c _tableau_ representing the data.
		 */
		std::vector<DataSeries<Xkind > > & refTable()
		{
			return _tableau_ ;
		}


		/** @brief Continue for another time measure ?
		 * @see TimeWatcher::keepon
		 * @param repet  current number of repetitions for this new measure
		 * @param tim    time previously spent on the measures.
		 * @return true if one more measure can be done
		 */
		bool keepon(index_t & repet, double tim, bool usePrediction=false)
		{
			return _time_watch_.keepon(repet,tim, usePrediction);
		}

		// tinyxml
		void exporte( const std::string & filename) const ;
		// tinyxml
		void importe( const std::string & filename) ;



	}; // PlotData


	class PlotStyle {
	public:
		//! What format the plot should be in?
		struct Term
		{
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
		struct Plot
		{
			//! Plot type
			enum Type {
				histo = 200, //! histogram plot. This is the default. x ticks are evenly spaced, whatever there value and are labelled with their value.
				graph = 201, //! standard plot. Plots y_i=f(x) with x in the first colum and y_i in ith column. x-ticks are well spaced.
				other = 202  //! other (ie user supplied).
			} ;
		};

		struct Line
		{
			enum Type {
				lines      = 300,  //! lines.
				histogram  = 301,  //! histogram (boxes).
				linespoints= 302,  //! lines with points. (default)
				points     = 303,  //! only the points.
				other      = 304   //! rien.
			} ;
		};

		struct Options
		{
			enum Type {
				oblique = 400,
				other   = 401
			};
		};

		/*! @brief Constructor.
		 * By default, creates an histogram representing the data in an eps plot.
		 */
		PlotStyle() :
			_term_(Term::eps),_plot_type_(Plot::histo),_line_type_(Line::histogram)
		{

		}

		/*! @brief sets the titles in the graph.
		 * @param titre Title of the graph
		 * @param titre_y Title of the y-axis (series)
		 * @param titre_x Title of the x-axis (data points)
		 */
		void setTitle ( std::string    titre
				, std::string  titre_y
				, std::string  titre_x)
		{
			_title_   = titre ;
			_title_x_ = titre_x ;
			_title_y_ = titre_y ;
		}

		/*! @brief Gets the title of the graph.
		 * @return a gnuplot command to set the title of the graph.
		 */
		std::string getTitle() const
		{
			std::string title = "#title\nset title \""  + _title_ + '\"';
			if (!_title_x_.empty())
				title +="\nset xlabel \"" + _title_x_ +'\"' ;
			if (!_title_y_.empty())
				title +="\nset ylabel \"" + _title_y_ +'\"' ;
			return title ;
		}

		/*! @brief Gets the title of points abscissa.
		 * @return a gnuplot command to set the title of the abscissa.
		 */
		std::string getTitleX() const
		{
			return "\nset xlabel \"" + _title_x_ + '\"' ;
		}

		/*! @brief Gets the title of the series.
		 * @return a gnuplot command to set the title of the ordinate (series).
		 */
		std::string getTitleY() const
		{
			return "\nset ylabel \"" + _title_y_ + '\"' ;
		}

		/*! @brief get the title string.
		 * @param index can be (0,1,2)
		 */
		std::string getRawTitle(int index=0) const
		{
			switch (index) {
			case 0 :
				return _title_ ;
			case 1 :
				return _title_x_ ;
			case 2 :
				return _title_y_ ;
			default :
				return "bad index" ;
			}
		}

		/*! @brief Sets the output format.
		 * @sa TermType
		 * @param term type
		 */
		void setTerm( enum Term::Type term)
		{
			_term_ = term ;
		}

#if 0
		/*! Set additionnal term features
		*/
		void setTermOption(std::string & opts)
		{
			_termopts_ = opts;
		}
#endif

		/*! @brief Gets the output format.
		 * @return string for setting the expected output format in gnuplot.
		 */
		std::string getTerm() const
		{
			std::string term = "#term\nset term " ;
			switch(_term_) {
			case (Term::png) :
				term += "png enhanced" ;
				break;
			case (Term::pdf) :
				std::cerr << "warning, pdf not really working for now" << std::endl;
				term += "postscript eps enhanced color" ;
				break;
			case (Term::eps) :
				term += "postscript eps enhanced color" ;
				break;
			case (Term::svg) :
				term += "svg" ;
				break;
			case (Term::other) :
			default :
				std::cerr  << " *** error ***" << std::endl << "No supported term set" << std::endl;
				term += "unknown" ;
			}
			return term ;
		}

		/*! @brief Gets the graph output extension.
		 * By default, this is ".eps".
		 * @return a string for this extension, including the sepatating dot.
		 *
		 */
		std::string getExt() const
		{
			switch(_term_) {
			case (Term::png) :
				return ".png" ;
			case (Term::pdf) :
#ifndef __LINBOX_HAVE_GHOSTSCRIPT
				std::cerr << "warning, pdf not available. falling back to eps" << std::endl;
#endif
				return ".pdf" ;
			case (Term::eps) :
				return ".eps" ;
			case (Term::svg) :
				return ".svg" ;
			default :
				std::cerr << "unknown extension set" << std::endl;
				return ".xxx" ;
			}
		}

#if 0
		/*! @brief gets the style of the graph.
		 * This is very user-tweakable !!
		 * @return the style for gnuplot.
		 */
		std::string getStyle()
		{
			return "#style\n"+_styleopts_ ;
		}

		/*! @brief sets the style of the graph.
		 * This is very user-tweakable !!
		 * @param style the style for gnuplot as a gnuplot command.
		 */
		void setStyle(std::string style)
		{
			_styleopts_ = style ;
		}
#endif

		/*! @brief sets the legend position.
		 * @param keypos the arguments to key (where the legend should be put)
		 * can be : inside, outside,...
		 */
		void setKeyPos(std::string keypos)
		{
			_legend_pos_ = keypos ;
		}


		/*! @brief Gets the legend position.
		 * by default, it is "under".
		 */
		std::string getKeyPos() const
		{
			std::string lgd ="#legend\nset key " ;
			if (!_legend_pos_.empty())
				lgd +=  _legend_pos_ ;
			else
				lgd += " under" ;
			return lgd;
		}

		/*! @brief sets the  position of the labels on the X absciss.
		 * @param opt
		 * @param more more stuff
		 */
		void setXtics ( enum Options::Type opt, std::string more="")
		{
			_xtics_ =  "#xtics\nset xtics ";
			if (opt == Options::oblique)
				_xtics_ +=  "nomirror rotate by -45 scale 0 ";
			else {
				linbox_check(opt == Options::other);
				_xtics_ += more ;
			}
		}

		/*! @brief Gets the legend position.
		 * by default, it is 45° inclined (use in on long tics legends).
		 */
		std::string getXtics() const
		{
			return _xtics_ ;
		}

		/*! @brief Gets the name of the output graph.
		 * @param basnam the raw name for the output.
		 * @return basnam+extenstion.
		 */
		std::string getOutput(std::string basnam) const
		{
			std::string setout = "#output\nset output \'" ;
#ifdef __LINBOX_HAVE_GHOSTSCRIPT
			if (_term_ == Term::pdf)
				setout += "| ps2pdf - " ;
			setout += basnam + getExt() + '\'' ;
#else
			setout += basnam + ".eps\'" ;
#endif

			return setout ;
		}

		/*! @brief Sets the type of plot.
		 * @param type the type.
		 * @sa PlotType
		 *
		 */
		void setPlotType(enum Plot::Type type)
		{
			_plot_type_ = type ;
			// _plot_extra_ = moreargs ;
		}

		/*! @brief Sets the way dots are linked.
		 * @sa LineType
		 * @param type type
		 */
		void setLineType( enum Line::Type type)
		{
			_line_type_ = type ;
		}

		/*! @brief Gets the type of plot.
		 * default is histogram, or if graph is supplied, then the default is linespoints.
		 * Can be totally customized.
		 * @return a string for gnuplot to set the plot type.
		 * @sa PlotType
		 *
		 */
		std::string getPlotType() // const
		{
			_styleopts_ += "\nset datafile missing \"inf\"" ;
			std::string mystyle = "#style\nset style data " ;
			if (_line_type_ != Line::other) {
				switch (_line_type_) {
				case (Line::lines) :
					mystyle += "lines" ;
					break;
				case (Line::histogram) :
					mystyle += "histogram" ;
					break;
				case (Line::points) :
					mystyle += "points" ;
					break;
				case (Line::linespoints) :
					mystyle += "linespoints" ;
					break;
				default :
					std::cout << __func__ << " : you should have set the LineType when ploting PlotType::graph !" << std::endl;
					mystyle += "other" ;
				}
			}
			else { // userd defined datastyle
				return _styleopts_  ;
			}
			// some more style args :
			mystyle += "\n" + _styleopts_ + "\n";
			return mystyle ;
		}

		/*! @brief adds some style line to the graph.
		 * This is very user-tweakable !!
		 * @param style a style line for gnuplot as a gnuplot command.
		 */
		void addPlotType(std::string style)
		{
			_styleopts_ += "\n" + style ;
		}

		/*! @brief tells which columns to use.
		 * @param col a column to use.
		 * @param moreargs more stuff
		 */
		void setUsingSeries(index_t col, std::string moreargs= "")
		{
			linbox_check(col>1);
			std::ostringstream usingcols ;
			if (_plot_type_ == Plot::histo) {
				usingcols << " using " << col << ":xtic(1) title columnheader(" << col << ") "  << moreargs << " ";
			}
			else {
				linbox_check(_plot_type_ == Plot::graph);
				usingcols << " using 1:" << col << " title columnheader(" << col << ") "  << moreargs << " ";
			}
			_usingcols_ = usingcols.str();
		}

		/*! @brief adds a column to use
		 * @param col a  column to use.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(index_t col, std::string moreargs= "")
		{
			linbox_check(col>2);
			linbox_check(!_usingcols_.empty());
			std::ostringstream usingcols ;
			usingcols << ", \'\' using " ;
			if (_plot_type_ == Plot::graph)
				usingcols << "1:" ;
			usingcols << col << " ti col "  << moreargs << " ";
			_usingcols_ += usingcols.str();
		}

		/*! @brief tells which columns to use.
		 * @param cols a list of column to use.
		 * @param moreargs more stuff
		 */
		void setUsingSeries(std::list<index_t> cols, std::string moreargs= "")
		{
			linbox_check(!cols.empty());
			std::list<index_t>::iterator it = cols.begin();
			// no way to check *it< coldim...
			std::ostringstream usingcols ;
			if (_plot_type_ == Plot::histo) {
				usingcols << " using " << *it << ":xtic(1) title columnheader(" << *it << ") " << moreargs << " " ;
				++it ;
				for (;it != cols.end();++it) {
					usingcols << ", \'\' using " << *it << " ti col "  << moreargs << " ";
				}
			}
			else {
				linbox_check(_plot_type_ == Plot::graph);
				usingcols << " using 1:" << *it << " title columnheader(" << *it << ") "  << moreargs << " ";
				++it ;
				for (;it != cols.end();++it) {
					usingcols << ", \'\' using 1:" << *it << " ti col "  << moreargs << " ";
				}
			}

			_usingcols_ = usingcols.str();
			return;
		}

		/*! @brief adds a set of columns to use.
		 * @param cols a list of column to use.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(std::list<index_t> cols, std::string moreargs= "")
		{
			linbox_check(!cols.empty());
			linbox_check(!_usingcols_.empty());
			std::list<index_t>::iterator it = cols.begin();
			std::ostringstream usingcols ;
			if (_plot_type_ == Plot::histo) {
				for (;it != cols.end();++it) {
					usingcols << ", \'\' using " << *it << " ti col "  << moreargs << " ";

				}
			}
			else {
				linbox_check(_plot_type_ == Plot::graph);
				for (;it != cols.end();++it) {
					usingcols << ", \'\' using 1:" << *it << " ti col "  << moreargs << " ";

				}
			}
			_usingcols_ += usingcols.str();
			return;
		}

		/*! @brief tells which columns to use.
		 * @param cols all colums between \c cols.first and \c cols.second
		 * will be used.
		 * @param moreargs more stuff
		 *
		 */
		void setUsingSeries(std::pair<index_t,index_t> cols, std::string moreargs= "")
		{
			std::ostringstream usingcols ;
			if (_plot_type_ == Plot::histo) {
				usingcols << " using " << cols.first << ":xtic(1) title columnheader(" << cols.first << ") "  << moreargs << " ";
				usingcols << ", for [i=" << cols.first+1 << ":" << cols.second << "] \'\' using i title columnheader(i) "  << moreargs << " ";
			}
			else {
				linbox_check(_plot_type_ == Plot::graph);
				usingcols << " using 1:" << cols.first << " title columnheader(" << cols.first << ") "  << moreargs << " ";
				usingcols << ", for [i=" << cols.first+1 << ":" << cols.second << "] \'\' using 1:i title columnheader(i) "  << moreargs << " ";
			}

			_usingcols_ = usingcols.str();
			return;

		}

		/*! @brief adds contiguous columns to use.
		 * @param cols all colums between \c cols.first and \c
		 * cols.second will be used.
		 * @param moreargs more stuff
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has
		 * already been called.
		 *
		 */
		void addUsingSeries(std::pair<index_t,index_t> cols, std::string moreargs= "")
		{
			linbox_check(!_usingcols_.empty());
			std::ostringstream usingcols ;
			if (_plot_type_ == Plot::histo) {
				usingcols << ", for i=[" << cols.first << ":" << cols.second << "] \'\' using i title columnheader(i) " << moreargs << " ";
			}
			else {
				usingcols << ", for i=[" << cols.first << ":" << cols.second << "] \'\' using 1:i title columnheader(i) "  << moreargs << " ";
				linbox_check(_plot_type_ == Plot::graph);
			}

			_usingcols_ += usingcols.str();

		}

		/*! @brief Gets the plot command line.
		 * @param File the name of/path to the data file (with extension)
		 * @return a gnuplot "plot" command stream.
		 */
		std::string getPlotCommand(std::string File) const
		{
			std::string PC = "#plot\nplot \'" + File + "\' ";
			PC += _usingcols_ ;
			return PC ;
		}

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


	template<class Xkind>
	class PlotGraph {
	private :
		PlotData<Xkind>          & _data_ ;   //!< reference to the data points
		PlotStyle             & _style_ ;   //!< reference to a plotting style
		std::string         _filename_  ;   //!< name for the output file (without extension). a random \c _XXXXXX suffix will be added to make it unique.
		dmatrix_t             _merge_data_   ;
		svector_t     _merge_points_ ;

		/*! @internal
		 * @brief random <code>:alnum:</code> \c char.
		 * [[:alnum:]] characters are in range
		 *  - num : 48-57
		 *  - AL  : 67-90
		 *  - al  : 97-122
		 *  .
		 * @return a random alphabetic or numeric char.
		 */
		char _randomAlNum()
		{
			int c = rand()%62 ;
			c += 48 ; // c entre 48 et 109
			if (c < 58) {
				return (char) c;
			}
			else {
				c += 7 ;
				if (c < 91) {
					return (char)c ;
				}
				else {
					c += 6 ;
					return (char)c;
				}
			}


		}

		/*! @internal
		 * @brief Appends random suffix.
		 * Appends to \p _filename_ a random string constituted of an
		 * underscore followed by 8 random alnum chars.
		 * @return the concatenation of \c _filename_ and this suffix.
		 */
		std::string _randomName()
		{
			std::ostringstream unique_filename ;
			unique_filename << _filename_ << '_' ;
			for (index_t i = 8 ; i-- ; ) {
				unique_filename << _randomAlNum() ;
			}
			// std::cout << unique_filename.str() << std::endl;
			return unique_filename.str() ;

		}


		//! @bug this supposes the two series have unique measurements for one point.
		void mergeTwoSeries( svector_t & merge_points
				     , dmatrix_t & merge_data
				     , const svector_t & pts
				     , const dvector_t & dat
				     , const index_t & idx) const
		{
			index_t data_size = (index_t)merge_points.size();
			linbox_check(data_size == (index_t)merge_data[0].size());

			merge_data[idx].resize(data_size,NAN);
			typename svector_t::iterator it ;

			for (index_t i = 0 ; i < pts.size() ; ++i) {
				// iterators change because of push_back, so they are here :
				typename svector_t::iterator beg = merge_points.begin() ;
				typename svector_t::iterator end = merge_points.begin()+data_size ;

				// std::cout << "inserting "<< pts[i] << std::endl;
				it = std::find( beg, end, pts[i] ) ;
				if (it != end){
					index_t j = (index_t) std::distance(beg,it);
					merge_data[idx][j] = dat[i] ;
				}
				else {
					for (index_t j = 0 ; j < idx ; ++j) {
						merge_data[j].push_back(NAN);
					}

					merge_data[idx].push_back(dat[i]) ;
					merge_points.push_back(pts[i]) ;
				}
				// std::cout << "..." << std::endl;
				// std::cout << merge_points << std::endl;
				// std::cout << merge_data << std::endl;
				// std::cout << "..." << std::endl;
			}

			return;

		}

		//! merge all series of points into a vector of absissa points  and a vector of vector of data points
		void mergeSeries()
		{
			_data_. selectSeries(0);
			_merge_points_ = _data_.getCurrentSeriesPointLabel() ;
			_merge_data_[0] = _data_.getCurrentSeriesValues() ;

			std::cout << "merge points " << _merge_points_ << std::endl;
			std::cout << "merge data   " << _merge_data_ << std::endl;

			for (index_t i = 1 ; i < _data_.size() ; ++i) {
				_data_. selectNextSeries() ;
				std::cout << "to be merged "  << i << " : "  << std::endl;
				std::cout << "new points " << _data_.getCurrentSeriesPointLabel() << std::endl;
				std::cout << "new data   " << _data_.getCurrentSeriesValues() << std::endl;

				mergeTwoSeries(_merge_points_,_merge_data_,
					       _data_. getCurrentSeriesPointLabel(), _data_. getCurrentSeriesValues(),i);

				std::cout << "result : " << std::endl;
				std::cout << "merge points " << _merge_points_ << std::endl;
				std::cout << "merge data   " << _merge_data_ << std::endl;

			}

			return ;
		}

		public :

		/*! @brief Sets a new data structure.
		 * @param data a reference to a PlotData class.
		 */
		void setData( PlotData<Xkind> & data )
		{
			_data_ = data ;
		}

		/*! @brief Gets the data.
		 * @param[in,out] data a reference to a PlotData class.
		 */
		PlotData<Xkind> & refData( PlotData<Xkind> & data)
		{
			return data = _data_ ;
		}

		/*! @brief Sets a new style structure.
		 * @param style a reference to a PlotStyle class.
		 */
		void setStyle( PlotStyle & style )
		{
			_style_ = style ;
		}

		/*! @brief Gets the style.
		 * @param[in,out] style a reference to a PlotStyle class.
		 */
		PlotStyle & refStyle( PlotStyle & style)
		{
			return style = _style_ ;
		}


		// not implemented yet
		void sortSeries()  ;

		// not implemented yet
		void unique() ;

		// metadata : machine (uname -a) ; date (date --rfc-3339=seconds) ; program name ; other
		void addMetadata() ;

			/*! @brief Constructor for the PlotGraph class.
		 * Plots a series of data according to a style.
		 * @param data data to be plot, will be processed by the style
		 * @param style sets parameters to gnuplot to achieve a nice
		 * plot.
		 */
		PlotGraph( PlotData<Xkind> & data, PlotStyle & style ) :
			_data_(data)
			,_style_(style)
			,_merge_data_(data.size())
			,_merge_points_(data.getSeries(0).size())
		{
			srand((unsigned)time(NULL));
			mergeSeries();
		}

		/*! @brief sets the ouput file name.
		 * All output is put in a "data" subfolder.
		 * @warning Since no file is overwritten, this
		 * directory can rapidly get very populated.
		 */
		void setOutFilename( std::string filename )
		{
			int err = system( "test -d data || ( rm -rf data && mkdir data )" ) ;
			if (err) {
				throw LinBoxError("could not create directory data");
			}

			if ( filename.empty() ) {
				_filename_ = "./data/plotdata" ;
				std::cerr << "you should provide a filename. Using " << _filename_ << " as default ."<<std::endl;
			}
			else {
				_filename_ = "./data/" + filename;
			}

		} ;

		//! @todo
		 void print_csv() ;

 		 //! @todo
		 void print_xml() ;

		//! @todo
		 void print_html() ;

		/*! @brief Prints data in a latex tabular.
		*/
		void print_latex()
		{
			index_t nb_points = (index_t)_merge_points_.size();
			index_t nb_series = _data_.size();

			linbox_check(nb_points);
			linbox_check(nb_series);
			// srand(time(NULL));
			// std::ostringstream unique_filename  ;
			std::string unique_filename = _randomName();
			unique_filename += ".tex" ;
			// std::cout << _filename_ << " plot in " << unique_filename << '.'<< std::endl;
			std::ofstream FN(unique_filename.c_str());
			//!@todo check FN opened.
			// begin
			FN << "%\\usepackage{slashbox}" << std::endl;
			FN << "\\begin{table}" << std::endl;
			FN << "\\centering"    << std::endl;
			// format
			FN << "\\begin{tabular}{c||" ;
			for (index_t j = nb_points ; j-- ; )
				FN << 'c' ;
			FN << "|}" << std::endl;
			// top left case
			std::string series = _style_.getRawTitle(2);
			std::string points = _style_.getRawTitle(1);
			if (!points.empty()) {
				FN << "\\backslashbox{" << points << "}{" << series << "}" ;
			}
			else {
				FN << series ;
			}
			// first line
			for (index_t j = 0 ; j < nb_points ; ++j ) {
				FN << " & " <<  _merge_points_[j] ;
			}
			// lines of data
			FN << std::endl << "\\hline" << std::endl;
			FN.precision(2);
			for (index_t i = 0 ; i < nb_series ; ++i) {
				FN << _data_.getSerieLabel(i) ;
				for (index_t j =  0 ; j < nb_points ; ++j )
					FN << " & " << _merge_data_[i][j] ;
				if (i+1 < nb_series )
					FN << "\\\\" ;
				FN << std::endl;
			}
			// end
			FN << "\\end{tabular}" << std::endl;
			FN << "\\caption{" << _style_.getRawTitle() << "}" << std::endl;
			FN << "\\label{tab:<+" << "label+>}" << std::endl;
			FN << "\\end{table}" << std::endl ;

			FN.close();

			std::cout << "latex table in " << unique_filename << '.' << std::endl;
			return ;

		}

		/*!@brief Plots the data with gnuplot.
		 * Produces data in a .dat file, creates a .gp gnuplot script and
		 * outputs a graph calling gnuplot.
		 * @warning If gnuplot is not available, fall back to the latex method.
		 */
		void print_gnuplot()
		{
#ifndef __LINBOX_HAVE_GNUPLOT
			std::cout << "gnuplot is not available on your system. using latex table as fallback" << std::endl;
			print_latex();
#else
			// srand(time(NULL));
			index_t nb_points = (index_t)_merge_points_.size() ;
			index_t nb_series = (index_t)_data_.size() ;

			std::string unique_filename  = _randomName();
			std::string DataFileName = unique_filename + ".dat" ;
			std::string PlotFileName = unique_filename + ".gp" ;
			std::ofstream DF(DataFileName.c_str());
			std::ofstream PF(PlotFileName.c_str());
			/*  Data file to be plot */
			// DF.precision(_style_.getPrecision());
			DF.precision(2);
			DF << "legend " ;
			for (index_t i = 0 ; i < nb_series ; ++i) {
				DF << _data_.getSerieLabel(i) << ' ' ;
			}
			DF << std::endl;

			for (index_t j = 0 ; j < nb_points ; ++j) {
				DF << _merge_points_[j] ;
				for (index_t i = 0 ; i < nb_series ; ++i) {
					DF << " " << _merge_data_[i][j] ;
				}
				DF << std::endl;
			}

			/*  Ploting script */
			PF << "#" << _filename_                    << std::endl;
			PF << _style_.getTerm()                    << std::endl;
			PF << _style_.getOutput(unique_filename)   << std::endl;
			PF << _style_.getTitle()                   << std::endl;
			PF << _style_.getKeyPos()                  << std::endl;
			PF << _style_.getXtics()                   << std::endl;
			PF << _style_.getPlotType()                << std::endl;
#if 0
			for (index_t i = 0 ; i < nb_series ; ++i) {
				PF << "set style line " << _style_.getLineStyle() << std::endl;
			}
#endif
			PF << _style_.getPlotCommand(DataFileName) << std::endl;
#if 0
			for (index_t i = 0 ; i < nb_series ; ++i) {
				PF << '\"' << DataFileName << "\" using 1:" << i+2 << " with lines " ;
				PF << " ls " << i+1 ;
				PF << " title '" << names[i] << "'";
				if (i < (nb_series-1))
					PF << ",\\" << std::endl;
				else
					PF << ';' ;
			}
			PF << std::endl;
#endif
			PF.close();

			std::string command( "gnuplot " ) ;
			command += PlotFileName ;
			int err = system( command.c_str() ) ;
			if (err) {
				std::cout << "errors have occured. Look at gnuplot output." << std::endl;
			}
			else {
				std::cout << "Output generated as " << unique_filename  + _style_.getExt() << std::endl;
			}


#endif
		}

	}; // PlotGraph

} // LinBox

//
// DataSeries
//
namespace LinBox {

	template<class Xkind>
	DataSeries<Xkind>:: DataSeries() :
		PointLabels(0)
		, Points(0)
		, Times(0)
		, Values(0)
	{}

	template<class Xkind>
	DataSeries<Xkind>::~DataSeries() {}

	template<class Xkind>
	void
	DataSeries<Xkind>::resize(const index_t & n)
	{
		linbox_check(n == Values.size()+1);
		PointLabels.resize(n);
		Times.resize(n);
		Points.resize(n);
		Values.resize(n);

		return;
	}

	template<class Xkind>
	index_t
	DataSeries<Xkind>::size() const
	{
		linbox_check(PointLabels.size() == Points.size())
		linbox_check(Times.size() == Points.size())
		linbox_check(Times.size() == Values.size())

		return (index_t)Values.size();
	}


} // LinBox




#endif // __LINBOX_benchmarks_benchmark_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
