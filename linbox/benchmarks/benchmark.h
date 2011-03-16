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

#ifndef __LINBOX_benchmark_H
#define __LINBOX_benchmark_H

#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/util/debug.h"
#include <string>
#include <fstream>

typedef uint32_t index_t ;

namespace LinBox
{


	/*! @brief Represents a table of values to plot.
	 * list of values are reprensented by vectors.  the table is a vector
	 * of these vectors.
	 *
	 * @warning NaN, inf are used as missing data. More genenally
	 * we could store data in strings.
	 */
	class PlotStyle {
	public:
		//! What format the plot should be in?
		enum TermType {
			png, //!< png. Portable Network Graphics file.
			pdf, //!< pdf. Portable Document Format actually, this is eps converted to pdf.
			eps, //!< eps. Encapsulated PostScript. Cool for inclusion in LaTex files. This is the default.
		       	svg, //!< sgv. Scalable Vector Graphics.
			none //!< don't know yet...
		} ;

		//! What style of graphic : histogram (bars), lines between points ?
		enum StyleType {
			histogram, //! histogram plot. bars
			lines,     //! points with lines. this is the default.
			points     //! only points
		} ;

		/*! @brief Constructor.
		 * By default, creates an eps plot with lines joining data points.
		 */
		PlotStyle() :
			_term_(eps),_style_(lines)
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
		std::string getTitle()
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
		std::string getTitleX()
		{
			return "\nset xlabel \"" + _title_x_ + '\"' ;
		}

		/*! @brief Gets the title of the series.
		 * @return a gnuplot command to set the title of the ordinate (series).
		 */
		std::string getTitleY()
		{
			return "\nset ylabel \"" + _title_y_ + '\"' ;
		}

		/*! @brief get the title string.
		 * @param index can be (0,1,2)
		 * @param out the corresponding string ( title , x title , y title)
		 */
		std::string getRawTitle(int index=0)
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
		void setTerm( enum TermType term)
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
		std::string getTerm()
		{
			std::string term = "#term\nset term " ;
			switch(_term_) {
			case (png) :
				term += "png enhanced" ;
				break;
			case (pdf) :
				std::cerr << "warning, pdf not really working for now" << std::endl;
				term += "postscript eps enhanced color" ;
				break;
			case (eps) :
				term += "postscript eps enhanced color" ;
				break;
			case (svg) :
				term += "svg" ;
				break;
			case (none) :
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
		std::string getExt()
		{
			switch(_term_) {
			case (png) :
				return ".png" ;
			case (pdf) :
#ifndef __LINBOX_HAVE_GHOSTSCRIPT
				std::cerr << "warning, pdf not available. falling back to eps" << std::endl;
#endif
				return ".pdf" ;
			case (eps) :
				return ".eps" ;
			case (svg) :
				return ".svg" ;
			default :
				std::cerr << "unknown extension set" << std::endl;
				return ".xxx" ;
			}
		}

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

		/*! @brief adds some style line to the graph.
		 * This is very user-tweakable !!
		 * @param style a style line for gnuplot as a gnuplot command.
		 */
		void addStyle(std::string style)
		{
			_styleopts_ += "\n" + style ;
		}

		/*! @brief sets the legend position.
		 * @param keypos the arguments to key (where the legend should be put)
		 * can be : inside, outside,...
		 */
		void setKeyPos(std::string keypos)
		{
			_legend_pos_ = keypos ;
		}

		/*! @brief Gets the legend position.
		 * by default, it is "on".
		 */
		std::string getKeyPos()
		{
			std::string lgd ="#legend\nset key " ;
			if (!_legend_pos_.empty())
				lgd +=  _legend_pos_ ;
			else
				lgd += " on default" ;
			return lgd;
		}

		/*! @brief Gets the name of the output graph.
		 * @param basnam the raw name for the output.
		 * @return basnam+extenstion.
		 */
		std::string getOutput(std::string basnam)
		{
			std::string setout = "#output\nset output \'" ;
#ifdef __LINBOX_HAVE_GHOSTSCRIPT
			if (_term_ == pdf)
				setout += "| ps2pdf - " ;
			setout += basnam + getExt() + '\'' ;
#else
			setout += basnam + ".eps\'" ;
#endif

			return setout ;
		}

		/*! @brief Sets the type of plot.
		 * @param type the type.
		 * @sa StyleType
		 *
		 */
		void setType(enum StyleType type)
		{
			_style_ = type ;
		}

		/*! @brief Sets the type of plot.
		 * default is points joined with colored lines.
		 * @return a string for gnuplot to set the plot type.
		 * @sa StyleType
		 *
		 */
		std::string getType()
		{
			std::string mystyle = "#style\n" ;
			switch(_style_) {
			case (histogram) :
				mystyle += "set style data histogram " ;
				break;
			case (lines) :
				mystyle += "set style data linespoints " ;
				break;
			case (points) :
				break;
			default :
				mystyle += "#style unknown." ;
				break;
			}
			mystyle += "\nset datafile missing \"inf\"" ;
			return mystyle ;
		}

		/*! @brief tells which columns to use.
		 * @param col a column to use.
		 */
		void setUsingSeries(index_t col)
		{
			linbox_check(col>1);
			std::ostringstream usingcols ;
			usingcols << " using " << col << ":xtic(1) title columnheader(" << col << ")" ;
			_usingcols_ = usingcols.str();
		}

		/*! @brief adds a column to use
		 * @param col a  column to use.
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(index_t col)
		{
			linbox_check(col>2);
			linbox_check(!_usingcols_.empty());
			std::ostringstream usingcols ;
			usingcols << ", \'\' u " << col << " ti col " ;
			_usingcols_ += usingcols.str();
		}

		/*! @brief tells which columns to use.
		 * @param cols a list of column to use.
		 */
		void setUsingSeries(std::list<index_t> cols)
		{
			linbox_check(!cols.empty());
			std::list<index_t>::iterator it = cols.begin();
			// no way to check *it< coldim...
			std::ostringstream usingcols ;
			usingcols << " using " << *it << ":xtic(1) title columnheader(" << *it << ")" ;
			++it ;
			for (;it != cols.end();++it) {
				usingcols << ", \'\' u " << *it << " ti col " ;

			}
			_usingcols_ = usingcols.str();


		}

		/*! @brief adds a set of columns to use.
		 * @param cols a list of column to use.
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 */
		void addUsingSeries(std::list<index_t> cols)
		{
			linbox_check(!cols.empty());
			linbox_check(!_usingcols_.empty());
			std::list<index_t>::iterator it = cols.begin();
			std::ostringstream usingcols ;
			for (;it != cols.end();++it) {
				usingcols << ", \'\' u " << *it << " ti col " ;

			}
			_usingcols_ += usingcols.str();
		}

		/*! @brief tells which columns to use.
		 * @param cols all colums between \c cols.first and \c cols.second
		 * will be used.
		 *
		 */
		void setUsingSeries(std::pair<index_t,index_t> cols)
		{
			std::ostringstream usingcols ;
			usingcols << " using " << cols.first << ":xtic(1) title columnheader(" << cols.first << ")" ;
			usingcols << ", for [i=" << cols.first+1 << ":" << cols.second << "] \'\' using i title columnheader(i) " ;

			_usingcols_ = usingcols.str();

		}

		/*! @brief adds contiguous columns to use.
		 * @param cols all colums between \c cols.first and \c cols.second
		 * will be used.
		 * @pre \p _usingcols_ is not empty, ie \c setUsingSeries has already been called.
		 *
		 */
		void addUsingSeries(std::pair<index_t,index_t> cols)
		{
			linbox_check(!_usingcols_.empty());
			std::ostringstream usingcols ;
			usingcols << ", for i=[" << cols.first << ":" << cols.second << "] \'\' using i title columnheader(i) " ;

			_usingcols_ += usingcols.str();

		}

		/*! @brief Gets the plot command line.
		 * @param File the name of/path to the data file (with extension)
		 * @return a gnuplot "plot" command stream.
		 */
		std::string getPlotCommand(std::string File)
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
		/*  units */
		// std::string                         _unit_      ;
		/*  terminal output */
		enum TermType                       _term_      ; //!< output data format.
		// std::string                         _termopts_  ;
		/*  plotting style */
		enum StyleType                      _style_     ; //!< line/bar style
		std::string                         _styleopts_ ; //!< gp style command.
		/*  columns to use */
		std::string                         _usingcols_ ; //!< columns to be used (gp command)


	} ;

	// the minimal data to plot.
	template<class NAM>
	class PlotData {
	private :
		std::vector< std::vector< float > >  _tableau_   ;   //!< data. \c _tableau_[i][j] is the \e jth value of serie \e i.
		index_t                              _nb_points_ ;   //!< number of points in each series. (size of \c _tableau_[0])
		index_t                              _nb_series_ ;   //!< number of series. (size of \c _tableau_)
		std::vector< std::string >           _serie_name_;   //!< name for each serie of points
		std::vector< NAM >                   _absci_name_;   //!< values of the x axis.
	public :
		/*! Inits a plot with series of data.
		 * @param nb_pts number of points in each serie.
		 * @param nb_srs number of series of points. Default is 1.
		 */
		PlotData(index_t nb_pts, index_t nb_srs=1) :
			_tableau_(nb_srs),_nb_points_(nb_pts),_nb_series_(nb_srs),_serie_name_(nb_srs),_absci_name_(nb_pts)
		{
			linbox_check(nb_srs);
			linbox_check(nb_pts);
			for (index_t i = 0 ; i < nb_srs ; ++i)
				_tableau_[i].resize(nb_pts);
		}

		//! destructor.
		~PlotData() {} ;

		//! copy constructor.
		//! @param PD a PlotData to copy.
		PlotData(const PlotData<NAM> & PD):
			_tableau_(PD.getTable()),_nb_points_(PD.getPointsDim()),_nb_series_(PD.getSeriesDim()),_serie_name_(PD.getSerieNames()),_absci_name_(PD.getAbsciNames())
		{

		}

		/*! @brief get the number of series.
		 * @return number of series.
		 */
		index_t getSeriesDim()
		{
			return _nb_series_ ;
		}

		/*! @brief get the common number of points in each serie.
		 * @return number of points.
		 */
		index_t getPointsDim()
		{
			return _nb_points_ ;
		}

		/*! @brief Sets the name of a serie.
		 * @param i index of the serie
		 * @param nom name of the serie
		 */
		void setSerieName(index_t i, std::string nom)
		{
			linbox_check(i<_nb_series_);
			_serie_name_[i] = nom ;
		}

		/*! @brief Sets the name of a point.
		 * @param i index for the the point
		 * @param nom name of the point
		 */
		void setAbsciName(index_t j, NAM nom)
		{
			linbox_check(j<_nb_points_);
			_absci_name_[j] = nom ;
		}

		/*! @brief gets the name of a serie.
		 * Defaults to \c "serie.i"
		 * @param i its index.
		 * @return its name.
		 */
		std::string getSerieName(index_t i)
		{
			linbox_check(i<_nb_series_);
			if (_serie_name_[i].empty()) {
				std::ostringstream emptytitle ;
				emptytitle << "serie." << i ;
				return emptytitle.str();
			}
			return(_serie_name_[i]);
		}

		/*! @brief gets the name of a point.
		 * @param j its index.
		 * @return its name.
		 * @warning no default. \c setAbsciName has to be used beforehands.
		 */
		NAM getAbsciName(index_t j)
		{
			linbox_check(j<_nb_points_);
			return(_absci_name_[j]) ;
		}

		/*! @brief gets all the names in the series.
		 * @return a vector of names.
		 */
		std::vector<std::string > getSerieNames()
		{
			return _serie_name_ ;
		}

		/*! @brief gets all the names in the points.
		 * @return a vector of names.
		 */
		std::vector<NAM > getAbsciNames()
		{
			return _absci_name_ ;
		}

		/*! @brief modifies the number of series.
		 * @param n the new number of series. Some data will be lost if n is smaller than
		 * the current size.
		 */
		void resizeSeries( index_t & n)
		{
			if (n<_nb_series_) {
				std::cerr  << "warning, you are dropping series" << std::endl;
			}
			_tableau_.resize(n);
			return;
		}

		/*! @brief modifies the number of points.
		 * @param n the new number of points in every series. Some data
		 * will be lost if n is smaller than the current size.
		 */
		void resizePoints( index_t & n)
		{
			if (n<_nb_points_) {
				std::cerr  << "warning, you are dropping points in the series" << std::endl;
			}
			for (index_t i = 0 ; i < _nb_series_ ; ++i)
				_tableau_[i].resize(n);
		}

		/*! @brief sets a new entry.
		 * @param i index of the serie
		 * @param j index of the point
		 * @param val value to be inserted.
		 */
		void setEntry(index_t i, index_t j, float val)
		{
			linbox_check(i<_nb_series_);
			linbox_check(j<_nb_points_);
			_tableau_[i][j] = val ;
			return ;
		}

		/*! @brief gets a value for an entry.
		 * @param i index of the serie
		 * @param j index of the point
		 * @return val value of point j in serie j.
		 */
		float getEntry(index_t i, index_t j)
		{
			linbox_check(i<_nb_series_);
			linbox_check(j<_nb_points_);
			return _tableau_[i][j] ;
		}

		/*! gets a reference to the array of data.
		 * @return a reference to the member \c _tableau_ representing the data.
		 */
		std::vector<std::vector< float > > & getTable()
		{
			return _tableau_ ;
		}

	};

	template<class NAM>
	class PlotGraph {
	private :
		PlotData<NAM>          & _data_ ;   //!< reference to the data points
		PlotStyle             & _style_ ;   //!< reference to a plotting style
		std::string         _filename_  ;   //!< name for the output file (without extension). a random \c _XXXXXX suffix will be added to make it unique.

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

		/*! @brief Sets a new data structure.
		 * @param data a reference to a PlotData class.
		 */
		void setData( PlotData<NAM> & data )
		{
			_data_ = data ;
		}

		/*! @brief Gets the data.
		 * @param[in,out] data a reference to a PlotData class.
		 */
		PlotData<NAM> & refData( PlotData<NAM> & data)
		{
			return data = _data_ ;
		}

		/*! @brief Sets a new style structure.
		 * @param style a reference to a PlotData class.
		 */
		void setStyle( PlotStyle & style )
		{
			_style_ = style ;
		}

		/*! @brief Gets the style.
		 * @param[in,out] style a reference to a PlotData class.
		 */
		PlotStyle & refStyle( PlotStyle & style)
		{
			return style = _style_ ;
		}

	public :

		/*! @brief Constructor for the PlotGraph class.
		 * Plots a series of data according to a style.
		 * @param data data to be plot, will be processed by the style
		 * @param style sets parameters to gnuplot to achieve a nice
		 * plot.
		 */
		PlotGraph( PlotData<NAM> & data, PlotStyle & style ) :
			_data_(data),_style_(style)
		{
			srand(time(NULL));
		}

		/*! @brief sets the ouput file name.
		 * All output is put in a "data" subfolder.
		 * @warning Since no file is overwritten, this
		 * directory can rapidly get very populated.
		 */
		void setOutFilename( std::string filename )
		{
			if ( filename.empty() ) {
				_filename_ = "./data/plotdata" ;
				std::cerr << "you should provide a filename. Using " << _filename_ << " as default ."<<std::endl;
			}
			else {
				_filename_ = "./data/" + filename;
			}

		} ;

		/*! @brief Prints data in a latex tabular.
		 */
		void print_latex()
		{
			index_t nb_points = _data_.getPointsDim();
			index_t nb_series = _data_.getSeriesDim();

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
			FN << "\%\\usepackage{slashbox}" << std::endl;
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
				FN << " & " <<  _data_.getAbsciName(j) ;
			}
			// lines of data
			FN << std::endl << "\\hline" << std::endl;
			FN.precision(2);
			for (index_t i = 0 ; i < nb_series ; ++i) {
				FN << _data_.getSerieName(i) ;
				for (index_t j =  0 ; j < nb_points ; ++j )
					FN << " & " << _data_.getEntry(i,j) ;
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
			index_t nb_points = _data_.getPointsDim() ;
			index_t nb_series = _data_.getSeriesDim() ;

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
				DF << _data_.getSerieName(i) << ' ' ;
			}
			DF << std::endl;

			for (index_t j = 0 ; j < nb_points ; ++j) {
				DF << _data_.getAbsciName(j) ;
				for (index_t i = 0 ; i < nb_series ; ++i) {
					DF << " " << _data_.getEntry(i,j) ;
				}
				DF << std::endl;
			}

			/*  Ploting script */
			PF << "#" << _filename_ << std::endl;
			PF << _style_.getTerm() << std::endl;
			PF << _style_.getOutput(unique_filename)  << std::endl;
			PF << _style_.getTitle() << std::endl;
			PF << _style_.getKeyPos() << std::endl ;
			PF << _style_.getType() << std::endl;
			PF << _style_.getStyle() << std::endl;
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

	};

}
#endif // __LINBOX_benchmark_H
