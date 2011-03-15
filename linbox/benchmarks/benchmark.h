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


	/*! @brief represents a table of values to plot.
	 * list of values are reprensented by vectors.  the table is a vector
	 * of these vectors.
	 *
	 * @warning NaN, inf are used as missing data instead of a more genenal
	 * <code> union { double ; char } </code> which is either
	 * { <double>, 'o'} (a float) or {0, 'x'} (a missing value). Could be a
	 * string as well.
	 */

	class PlotStyle {
	public:
		enum TermType { png, pdf, eps, svg } ;
		enum StyleType {histogram,lines} ;

		PlotStyle()
		{

		}

		void setTitle ( std::string  titre
				, std::string  titre_y
				, std::string  titre_x)
		{
			_title_ = titre ;
			_title_x_ = titre_x ;
			_title_y_ = titre_y ;
		}

		std::string getTitle()
		{
			std::string title = "#title\nset title \"" ;
			title += _title_ ;
			title += "\"";
			return title ;
		}
		std::string getTitleX()
		{
			return _title_x_ ;
		}
		std::string getTitleY()
		{
			return _title_y_ ;
		}

		std::string getUnit()
		{
			return _unit_ ;
		}

		void setUnit(std::string unit)
		{
			_unit_ = unit ;
		}

		void setTerm( enum TermType term)
		{
			_term_ = term ;
		}

		void setTermOption(std::string & opts)
		{
			_termopts_ = opts;
		}

		std::string getTerm()
		{
			switch(_term_) {
			case (png) :
				return "#term\nset term png enhanced" ;
			case (pdf) :
				std::cerr << "warning, pdf not really working for now" << std::endl;
				return "#term\nset term postscript eps enhanced color" ;
			case (eps) :
				return "#term\nset term postscript eps enhanced color" ;
			case (svg) :
				return "#term\nset term svg" ;
			default :
				std::cerr  << "unknown term set" << std::endl;
				return "set unknown term" ;
			}
		}

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

		std::string getStyle()
		{
			return "#style\n"+_styleopts_ ;
		}

		void setStyle(std::string style)
		{
			_styleopts_ = style ;
		}
		void addStyle(std::string style)
		{
			_styleopts_ += "\n" + style ;
		}

		void setKeyPos(std::string keypos)
		{
			_legend_pos_ = keypos ;
		}

		std::string getKeyPos()
		{
			std::string lgd ="#legend\nset key " ;
			if (!_legend_pos_.empty())
				lgd +=  _legend_pos_ ;
			else
				lgd += " on default" ;
			return lgd;
		}

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

		void setType(enum StyleType type)
		{
			_style_ = type ;
		}

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
			default :
				mystyle += "#style unknown." ;
				break;
			}
			mystyle += "\nset datafile missing \"inf\"" ;
			return mystyle ;
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
			usingcols << " u" << *it << ":xtic(1) title columnheader(" << *it << ")" ;
			++it ;
			for (;it != cols.end();++it) {
				usingcols << ", \'\' u " << *it << " ti col " ;

			}
			_usingcols_ = usingcols.str();


		}

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

		void addUsingSeries(std::pair<index_t,index_t> cols)
		{
			linbox_check(!_usingcols_.empty());
			std::ostringstream usingcols ;
			usingcols << ", for i=[" << cols.first << ":" << cols.second << "] \'\' using i title columnheader(i) " ;

			_usingcols_ += usingcols.str();

		}



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
		std::string                         _unit_      ;
		/*  terminal output */
		enum TermType                       _term_      ;
		std::string                         _termopts_  ;
		/*  plotting style */
		enum StyleType                      _style_     ;
		std::string                         _styleopts_ ;
		/*  columns to use */
		std::string                         _usingcols_ ;


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
			for (index_t i = 0 ; i < nb_srs ; ++i)
				_tableau_[i].resize(nb_pts);
		}

		//! destructor.
		~PlotData() {} ;

		//! copy constructor.
		PlotData(const PlotData & PD):
			_tableau_(PD.getTable()),_nb_points_(PD.getPointsDim()),_nb_series_(PD.getSeriesDim()),_serie_name_(PD.getSerieNames()),_absci_name_(PD.getAbsciNames())
		{

		}

		index_t getSeriesDim()
		{
			return _nb_series_ ;
		}

		index_t getPointsDim()
		{
			return _nb_points_ ;
		}

		void setSerieName(index_t i, std::string nom)
		{
			linbox_check(i<_nb_series_);
			_serie_name_[i] = nom ;
		}

		void setAbsciName(index_t j, NAM nom)
		{
			linbox_check(j<_nb_points_);
			_absci_name_[j] = nom ;
		}

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

		NAM getAbsciName(index_t j)
		{
			linbox_check(j<_nb_points_);
			return(_absci_name_[j]) ;
		}

		std::vector<std::string > getSerieNames()
		{
			return _serie_name_ ;
		}

		std::vector<NAM > getAbsciNames()
		{
			return _absci_name_ ;
		}

		void resizeSeries( index_t & n)
		{
			if (n<_nb_series_) {
				std::cerr  << "warning, you are dropping series" << std::endl;
			}
			_tableau_.resize(n);
			return;
		}

		void resizePoints( index_t & n)
		{
			if (n<_nb_points_) {
				std::cerr  << "warning, you are dropping points in the series" << std::endl;
			}
			for (index_t i = 0 ; i < _nb_series_ ; ++i)
				_tableau_[i].resize(n);
		}

		void setEntry(index_t i, index_t j, float val)
		{
			linbox_check(i<_nb_series_);
			linbox_check(j<_nb_points_);
			_tableau_[i][j] = val ;
			return ;
		}

		float getEntry(index_t i, index_t j)
		{
			linbox_check(i<_nb_series_);
			linbox_check(j<_nb_points_);
			return _tableau_[i][j] ;
		}

		std::vector<std::vector< float > > & getTable()
		{
			return _tableau_ ;
		}

	};

	template<class NAM>
	class PlotGraph {
	private :
		PlotData<NAM>          & _data_ ;
		PlotStyle             & _style_ ;
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

	public :

		PlotGraph( PlotData<NAM> & data, PlotStyle & style ) :
			_data_(data),_style_(style)
		{
			srand(time(NULL));
		}

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

		// ~PlotData() {} ;

		// no copy constructor.

		/*! Prints data in a latex tabular.
		 * @todo 1ère case : titre_y,titre_y ;
		 * @todo unités.
		 */
		void print_latex()
		{
			index_t nb_points = _data_.getPointsDim();
			index_t nb_series = _data_.getSeriesDim();

			linbox_check(nb_points);
			linbox_check(nb_series_);
			// srand(time(NULL));
			// std::ostringstream unique_filename  ;
			std::string unique_filename = _randomName();
			unique_filename += ".tex" ;
			// std::cout << _filename_ << " plot in " << unique_filename << '.'<< std::endl;
			std::ofstream FN(unique_filename.c_str());
			//!@todo check FN opened.
			// begin
			FN << "\\begin{table}" << std::endl;
			FN << "\\centering"    << std::endl;
			// format
			FN << "\\begin{tabular}{c||" ;
			for (index_t j = nb_points ; j-- ; )
				FN << 'c' ;
			FN << "|}" << std::endl;
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
			// commenting out what is not printed yet.
			FN << "\% unit is : "      <<  _style_.getUnit() << std::endl;
			FN << "\% title is : "     <<  _style_.getTitle() << std::endl;
			FN << "\% row title is : " <<  _style_.getTitleX() << std::endl;
			FN << "\% col title is : " <<  _style_.getTitleY() << std::endl;
			// end
			FN << "\\end{tabular}" << std::endl;
			FN << "\\caption{<+" << "Caption text+>}" << std::endl;
			FN << "\\label{tab:<+" << "label+>}" << std::endl;
			FN << "\\end{table}" << std::endl ;

			FN.close();

			std::cout << "latex table in " << unique_filename << '.' << std::endl;
			return ;

		}

		/*!Plots the data with gnuplot.
		 * Produces data in a .dat file, creates a .gp gnuplot script and
		 * outputs a graph.
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
