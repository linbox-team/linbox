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

/*!@internal
 * @file   benchmarks/benchmark.C
 * @ingroup benchmarks
 * @brief   utils
  */

#ifndef __LINBOX_benchmarks_benchmark_C
#define __LINBOX_benchmarks_benchmark_C

#include "linbox-config.h"

#include "benchmark.h"

//
// PlotData
//

namespace LinBox {
#ifdef __LINBOX_HAVE_TINYXML2
	tinyxml2::XMLElement * PlotData::saveData(tinyxml2::XMLDocument & doc)
	{
		using namespace tinyxml2;
		XMLElement * data = doc.NewElement( "data" );

		selectFirstSeries();
		for (index_t i = 0 ; i < size() ; ++i  ) {

			XMLElement * serie = doc.NewElement ( "serie" );
			serie->SetAttribute("name",unfortifyString(getCurrentSerieName()).c_str());

			for (index_t j = 0 ;  j < getCurrentSerieSize() ; ++j)
			{
				XMLElement * point = doc.NewElement ( "point" );
				point->SetAttribute("x",unfortifyString(getCurrentSeriesEntry(j,Point::Labels() )).c_str());
				point->SetAttribute("y",getCurrentSeriesEntry(j,Point::Values() ));
				point->SetAttribute("time",getCurrentSeriesEntry(j,Point::Times() ));
				point->SetAttribute("xval",getCurrentSeriesEntry(j,Point::Points() ));
				point->SetAttribute("id",getCurrentSeriesId(j).c_str());

				serie->InsertEndChild( point );
			}

			data->InsertEndChild( serie );

			selectNextSeries();
		}

		return data ;
	}
#endif


	PlotData::PlotData() :
		_tableau_      (0)
		,_serie_label_ (0)
		,_curr_serie_  ( )
		,_time_watch_  ( )
	{
	}

	PlotData::~PlotData() {}

	PlotData::PlotData(const PlotData & PD):
		_tableau_(PD.getTable())
		,_serie_label_(PD.getSerieLabels())
		,_curr_serie_(PD.getCurrentSerieNumber())
		,_time_watch_  (_tableau_[_curr_serie_].Points,_tableau_[_curr_serie_].Times)
	{
	}

	index_t PlotData::selectIndex(const std::string & nom)
	{

		std::string nomf = fortifyString(nom);
		index_t j ;
		bool ok = findKeyword(j,_serie_label_.begin() , _serie_label_.end() , nomf);


#if 1
		if ( ! ok ) {
			linbox_check(j ==(index_t)_serie_label_.size() );
			_serie_label_.push_back(fortifyString(nom));
			_tableau_.resize(j+1);
		}
#else
		linbox_check(ok);
#endif

		initWatch(j);

		return j ;
	}

	index_t PlotData::getIndex(const std::string & nom) const
	{

		std::string nomf = fortifyString(nom);
		index_t j ;
		bool ok = findKeyword(j,_serie_label_.begin() , _serie_label_.end() , nomf);

		linbox_check(ok);


		return j ;
	}

	void PlotData::clear()
	{
		_tableau_.resize(0);
		_serie_label_.resize(0);
		_curr_serie_ =  0;
		_time_watch_.clear();
	}

	void PlotData::merge(const PlotData &PD)
	{
		for (index_t i = 0 ; i < (index_t)PD.size() ; ++i) {
			_tableau_.push_back(PD.getSeries(i));
			_serie_label_.push_back(fortifyString(PD.getSerieName(i)));
		}
		return ;
	}

	index_t PlotData::size() const
	{
		linbox_check(_tableau_.size() == _serie_label_.size());
		return (index_t)_tableau_.size() ;
	}

	index_t PlotData::getCurrentSerieNumber() const
	{
		return _curr_serie_ ;
	}

	void PlotData::setSerieName(const index_t & i, const std::string & nom)
	{
		linbox_check(i<size());
		_serie_label_[i] = fortifyString(nom) ;
	}

	const std::string & PlotData::getSerieName(const index_t & i) const
	{
		linbox_check(i<size());
		return _serie_label_[i] ;
	}

	const std::string & PlotData::getCurrentSerieName() const
	{
		return getSerieName(_curr_serie_) ;
	}

	void PlotData::setCurrentSerieName(const std::string & nom)
	{
		return setSerieName(_curr_serie_,nom);
	}

	void PlotData::initWatch ( const index_t & i)
	{
		linbox_check(i < size());
		_curr_serie_ = i ;
		_time_watch_.init(_tableau_[i].Points,_tableau_[i].Times);
	}

	void PlotData::initCurrentSeriesWatch ()
	{
		initWatch(_curr_serie_);
	}

	void PlotData::newSerie(const std::string & nom )
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

	void PlotData::finishSerie()
	{
	}

	const DataSeries & PlotData::getSeries(const index_t  &i) const
	{
		linbox_check(i < size());
		return _tableau_[i] ;
	}

	const DataSeries & PlotData::selectSeries(const index_t  &i)
	{
		initWatch(i);
		return getSeries(i) ;
	}

	const DataSeries & PlotData::selectSeries(const std::string & name)
	{
		return selectSeries(selectIndex(name));
	}

	DataSeries & PlotData::refSeries(const index_t  &i)
	{
		linbox_check(i < size());
		return _tableau_[i] ;
	}

	DataSeries & PlotData::refSeries(const std::string &nom)
	{
		return refSeries(getIndex(nom));
	}

	const DataSeries & PlotData::getCurrentSeries() const
	{
		return getSeries(_curr_serie_);
	}

	DataSeries & PlotData::refCurrentSeries()
	{
		return refSeries(_curr_serie_);
	}

	index_t PlotData::getSerieSize(const index_t & i) const
	{
		return getSeries(i).size();
	}

	index_t PlotData::getCurrentSerieSize() const
	{
		return getSeries(_curr_serie_).size();
	}


	void PlotData::selectFirstSeries()
	{
		selectSeries(0);

		return;
	}

	bool PlotData::selectNextSeries()
	{
		linbox_check(_curr_serie_ < size());
		++_curr_serie_ ;
		if (_curr_serie_ < size()) {
			selectSeries(_curr_serie_);
			return true;
		}
		return false;
	}

	void PlotData::setCurrentSeriesPointLabel(const index_t & j, const std::string & nom)
	{
		return setSeriesPointLabel(_curr_serie_,j,nom);
	}


	std::string  PlotData::getSerieLabel(const index_t & i) const
	{
		linbox_check(i<size());
		linbox_check(i<_serie_label_.size() );
		linbox_check(! _serie_label_[i].empty()) ;

		return(_serie_label_[i]);
	}

	const svector_t &  PlotData::getSerieLabels() const
	{
		return _serie_label_ ;
	}


	void PlotData::setSeriesEntry(const index_t &i, const std::string & nam, const double & val
				      , const double & xval , const double & yval)
	{
		refSeries(i).push_back(fortifyString(nam),val,xval,yval);
		initWatch(i); //  in case series has changed
		return ;
	}

	void PlotData::setSeriesEntry(const std::string & nom, const std::string & nam, const double & val
				, const double & xval, const double & yval)
	{
		selectSeries(nom);
		return setCurrentSeriesEntry(nam,val,xval,yval);
	}

	const std::vector<DataSeries > & PlotData::getTable() const
	{
		return _tableau_ ;
	}

	std::vector<DataSeries > & PlotData::refTable()
	{
		return _tableau_ ;
	}

	bool PlotData::keepon(index_t & repet, double tim, bool usePrediction)
	{
		return _time_watch_.keepon(repet,tim, usePrediction);
	}

	void PlotData::load( const std::string & filename)
	{
#ifdef __LINBOX_HAVE_TINYXML2
		using namespace tinyxml2;
		XMLDocument doc;
		doc.LoadFile( filename.c_str() );
		// std::cout << "loaded " << filename << std::endl;
		linbox_check(!doc.ErrorID());

		XMLElement * bench = doc.FirstChildElement( "benchmark");
		linbox_check(bench);
		XMLElement* data = bench -> FirstChildElement( "data" ) ;
		linbox_check(data);
		XMLElement* series = data -> FirstChildElement( "serie" ) ;

		clear();

		while (series)
		{
			newSerie( series->Attribute( "name" ) );
			XMLElement * points = series->FirstChildElement() ;
			while (points) {
				std::string x = points->Attribute( "x" );
				double      y = points->DoubleAttribute( "y" );
				double      time = points->DoubleAttribute( "time" );
				double      xval = points->DoubleAttribute( "xval" );
				// std::string id   = points->Attribute( "id" );

				setCurrentSeriesEntry(x,y,xval,time);

				// setCurrentSeriesEntryId(id);

				points = points->NextSiblingElement();
			}

			series = series->NextSiblingElement();
		}
#else
		throw LinBoxError("You need tinyxml2 for loading data");
#endif
		// save("toto.xml");

	}

	void PlotData::save( const std::string & filename
			     , const std::string & title
			     , const std::string & xtitle
			     , const std::string & ytitle )
	{
#ifdef __LINBOX_HAVE_TINYXML2
		using namespace tinyxml2;
		XMLDocument doc;

		doc.InsertEndChild(doc.NewDeclaration());

		XMLElement * benchmark = doc.NewElement( "benchmark" );
		doc.InsertEndChild(benchmark);

		{ // Benchmark Metadata
			XMLElement * metadata = doc.NewElement( "metadata" );

			smatrix_t uname = getMachineInformation();

			for (size_t i = 0 ; i < uname[0].size() ; ++i) {
				metadata->SetAttribute(uname[0][i].c_str(),uname[1][i].c_str());
			}

			std::string myTime = getDateTime();

			metadata->SetAttribute("time",myTime.c_str());

			benchmark->InsertEndChild(metadata);
		}

		{ // Legende
			XMLElement * legende = doc.NewElement( "legende" );

			std::string mytitle = title;
			if (title.empty())
				mytitle =filename ;
			legende->SetAttribute("title",unfortifyString(mytitle).c_str());
			if (!xtitle.empty())
				legende->SetAttribute("X",unfortifyString(xtitle).c_str());
			if (!ytitle.empty())
				legende->SetAttribute("Y",unfortifyString(ytitle).c_str());


			benchmark->InsertEndChild(legende);
		}


		{ // series
			XMLElement * data = saveData(doc);

			benchmark->InsertEndChild(data);
		}

		{ // point metadata
			XMLElement * metapoint = doc.NewElement( "PointMetaData" );
			for (index_t i = 0 ; i < _meta_data_.MetaDataVec.size() ; ++i)
			{
				XMLElement * item = NULL ;
				_meta_data_.MetaDataVec[i].writeMetaData(&item,doc);
				std::string pts = "";
				for (index_t j = 0 ; j < _meta_data_.MetaDataIDs[i].size(); ++j){
					pts += _meta_data_.MetaDataIDs[i][j] ;
					if (j+1 < _meta_data_.MetaDataIDs[i].size()  )
						pts += ',';
				}
				item->SetAttribute("used_in",pts.c_str());
				metapoint->InsertEndChild(item);
			}

			benchmark->InsertEndChild(metapoint);
		}

		doc.SaveFile(filename.c_str());

		std::cout << "xml table saved in " << filename << std::endl;
#else
		std::cout << "tinyxml2 is not installed, could not save" << std::endl;
#endif

	}

} // LinBox


//
// PlotStyle
//

namespace LinBox {


	PlotStyle::PlotStyle() :
		_term_(Term::eps),_plot_type_(Plot::histo),_line_type_(Line::histogram)
	{

	}

	void PlotStyle::setTitle ( const std::string  &  titre
				   , const std::string  & titre_y
				   , const std::string  & titre_x)
	{
		_title_   = titre ;
		_title_x_ = titre_x ;
		_title_y_ = titre_y ;
	}

	std::string PlotStyle::getTitle() const
	{
		std::string title = "#title\nset title \""  + _title_ + '\"';
		if (!_title_x_.empty())
			title +="\nset xlabel \"" + _title_x_ +'\"' ;
		if (!_title_y_.empty())
			title +="\nset ylabel \"" + _title_y_ +'\"' ;
		return title ;
	}

	std::string PlotStyle::getTitleX() const
	{
		return "\nset xlabel \"" + _title_x_ + '\"' ;
	}

	std::string PlotStyle::PlotStyle::getTitleY() const
	{
		return "\nset ylabel \"" + _title_y_ + '\"' ;
	}

	std::string PlotStyle::getRawTitle(int index) const
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

	void PlotStyle::setTerm( enum Term::Type term)
	{
		_term_ = term ;
	}

	std::string PlotStyle::getTerm() const
	{
		std::string term = "#term\nset term " ;
		switch(_term_) {
		case (Term::png) :
			term += "png noenhanced" ;
			break;
		case (Term::pdf) :
			std::cerr << "warning, pdf not really working for now" << std::endl;
			term += "postscript eps noenhanced color" ;
			break;
		case (Term::eps) :
			term += "postscript eps noenhanced color" ;
			break;
		case (Term::epstex) :
			term += "epslatex color colortext" ;
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

	std::string PlotStyle::getExt() const
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
		case (Term::epstex) :
			return ".tex" ;
		case (Term::svg) :
			return ".svg" ;
		default :
			std::cerr << "unknown extension set" << std::endl;
			return ".xxx" ;
		}
	}

	void PlotStyle::setKeyPos(const std::string & keypos)
	{
		_legend_pos_ = keypos ;
	}

	std::string PlotStyle::getKeyPos() const
	{
		std::string lgd ="#legend\nset key " ;
		if (!_legend_pos_.empty())
			lgd +=  _legend_pos_ ;
		else
			lgd += " under" ;
		return lgd;
	}

	void PlotStyle::setXtics ( enum Options::Type opt, const std::string & more )
	{
		_xtics_ =  "#xtics\nset xtics ";
		if (opt == Options::oblique)
			_xtics_ +=  "nomirror rotate by -45 scale 0 ";
		else {
			linbox_check(opt == Options::other);
			_xtics_ += more ;
		}
	}

	const std::string & PlotStyle::getXtics() const
	{
		return _xtics_ ;
	}

	std::string PlotStyle::getOutput(const std::string & basnam) const
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

	void PlotStyle::setPlotType(enum Plot::Type type)
	{
		_plot_type_ = type ;
		// _plot_extra_ = moreargs ;
	}

	void PlotStyle::setLineType( enum Line::Type type)
	{
		_line_type_ = type ;
	}

	std::string PlotStyle::getPlotType(const std::string & extraargs) // const
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
				if (extraargs.empty()) // default style
					mystyle += "\nset style histogram cluster gap 1\nset style fill solid border rgb \"black\"";
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
		mystyle += "\n" + _styleopts_ + "\n" + extraargs + "\n";
		return mystyle ;
	}

	void PlotStyle::addPlotType(const std::string & style)
	{
		_styleopts_ += "\n" + style ;
	}

	void PlotStyle::setUsingSeries(index_t col, const std::string & moreargs)
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

	void PlotStyle::addUsingSeries(index_t col, const std::string & moreargs)
	{
		linbox_check(col>2);
		linbox_check(!_usingcols_.empty()); // we don't add if nothing was set
		std::ostringstream usingcols ;
		usingcols << ", \'\' using " ;
		if (_plot_type_ == Plot::graph)
			usingcols << "1:" ;
		usingcols << col << " ti col "  << moreargs << " ";
		_usingcols_ += usingcols.str();
	}

	void PlotStyle::setUsingSeries(std::list<index_t> cols, const std::string & moreargs)
	{
		linbox_check(!cols.empty());
		std::list<index_t>::iterator it = cols.begin();
		// no way to check *it< coldim...
		std::ostringstream usingcols ;
		if ( _plot_type_ == Plot::histo ) {
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

	void PlotStyle::addUsingSeries(std::list<index_t> cols, const std::string & moreargs)
	{
		linbox_check(!cols.empty());
		linbox_check(!_usingcols_.empty()); // we don't add if nothing was set
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

	void PlotStyle::setUsingSeries(std::pair<index_t,index_t> cols, const std::string & moreargs)
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

	void PlotStyle::addUsingSeries(std::pair<index_t,index_t> cols, const std::string & moreargs)
	{
		linbox_check(!_usingcols_.empty()); // we don't add if nothing was set
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

	const std::string & PlotStyle::getUsingSeries() const
	{
		return _usingcols_ ;
	}

} // LinBox

//
// Curve fitting
//

namespace LinBox {

	// fit X[nn-1,nn],Y[nn-1,nn] and return evaluation at x.
	double fit2(const dvector_t & X, const dvector_t & Y, int nn, double x)
	{
		index_t n = (index_t) nn;
		assert(n>0);
		if ( n==1 ) {
			if ( X[0]==X[1] ) {
				// std::cerr << "two of your evaluation points are at the same X" << std::endl;
				// this is NOT supposed to happen.
				return (Y[0]+Y[1])/2 ;
			}
		}
		if (X[n]==X[n-1]) // discard the last one.
			return fit2(X,Y,(int)n-1,x);

		double a = (Y[n-1]-Y[n])/(X[n-1]-X[n]) ;
		double b = (X[n-1]*Y[n]-X[n]*Y[n-1])/(X[n-1]-X[n]) ;
		return a*x+b ;
	}

#ifdef __LINBOX_HAVE_LAPACK

	double fit_lapack3(const dvector_t &X, const dvector_t &Z, double x)
	{

		dvector_t Y = Z ;

		int n = (int) Z.size();
		linbox_check((size_t)n == X.size());

		int deg = (int) std::min((int)4,n);
		dvector_t V((size_t)(deg*n));

		int ldv = deg ;

#if 0 // Clapack (not working)
		for(index_t i = 0 ; i < (index_t)n; ++i) {
			for (index_t j = 0 ; j < (index_t)ldv; ++j) {
				V[i*ldv+j] = std::pow(X[i],j);
			}
		}

		clapack_dgels(CblasRowMajor, CblasNoTrans, n, deg, 1, &V[0],
			      deg, &Y[0], 1);

#endif

#if 1 /* basic least squares */
		// std::cout << V.size() << std::endl;
		for(index_t i = 0 ; i < (index_t)n; ++i) {
			for (index_t j = 0 ; j < (index_t)ldv; ++j) {
				V[i+j*(index_t)n] = std::pow(X[i],j);
			}
		}

		// std::cout << V << std::endl;

		int info;
		int ldun = 1 ;
		{

			int lwork = 2*n*deg*4 ;
			dvector_t work((size_t)lwork);
			char N[] = "N";
			dgels_(N, &n, &ldv, &ldun, &(V[0]) , &n, &(Y[0]), &n, &work[0], &lwork, &info);
		}

#endif

#if 0 /* least squares using SVN and V not nec. full rank */
		{
			int lwork = 2*deg+std::max(2*deg,n)*4 ;
			dvector_t work(lwork);
			dvector_t s(deg);
			double rcond = 1e-8 ;
			int rank ;
			dgelss_( &n, &ldv, &ldun, &(V[0]) , &n, &(Y[0]), &n, &s[0], &rcond, &rank, &work[0], &lwork, &info);
		}
#endif

#if 0 /* weighted least squares */
		//DGGGLM

		// TODO

#endif


		// std::cout << Y << std::endl;
		// horner eval the poly
		double res = 0.0;

		for(int i=deg-1; i >= 0; i--) {
			res = res * x + Y[(index_t)i];
		}
		return res;
	}
#endif // __LINBOX_HAVE_LAPACK


	double fit3(const dvector_t & X, const dvector_t & Y,int n, double x)
	{
#ifndef __LINBOX_HAVE_LAPACK /* à la main */
		linbox_check(n>1);
		linbox_check((size_t)n< X.size());
		linbox_check((size_t)n< Y.size());
		if (n==2) {
			if (X[1]==X[2])
				return fit2(X,Y,1,x) ;
			if (X[0]==X[2]) {
				return fit2(X,Y,1,x) ;
			}
			if (X[0]==X[1]) {
				dvector_t X1(2); X1[0]=X[1]; X1[2]=X[2];
				dvector_t Y1(2); Y1[0]=Y[1]; Y1[2]=Y[2];
				return fit2(X1,Y1,1,x) ;
			}
		}
		if (X[n]==X[n-1]) { // discard last
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));

			return fit3(X1,Y1,n-1,x) ;
		}
		if (X[n]==X[n-2]) { // discard last
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));

			return fit3(X1,Y1,n-1,x) ;
		}
		if (X[n-1]==X[n-2]) { // discard last but one
			dvector_t X1(X.begin(),X.begin()+(n-1));
			dvector_t Y1(Y.begin(),Y.begin()+(n-1));
			X1[n-1]=X[n];
			Y1[n-1]=Y[n];

			return fit3(X1,Y1,n-1,x) ;
		}

		// todo: use Lagrange ?
		// std::cout << X[n-2] << ',' << X[n-1] << ',' << X[n] << std::endl;
		double d  = (-X[n]+X[n-1])*(-X[n]+X[n-2])*(X[n-2]-X[n-1]) ;
		double a1 = -X[n]*Y[n-2]+X[n-2]*Y[n]+X[n-1]*Y[n-2]-X[n-1]*Y[n]+X[n]*Y[n-1]-X[n-2]*Y[n-1];
		double a2 = -X[n-2]*X[n-2]*Y[n]+X[n-2]*X[n-2]*Y[n-1]+X[n-1]*X[n-1]*Y[n]-Y[n-2]*X[n-1]*X[n-1]+Y[n-2]*X[n]*X[n]-Y[n-1]*X[n]*X[n];
		double a3 = X[n-2]*X[n-2]*X[n-1]*Y[n]-X[n-2]*X[n-2]*X[n]*Y[n-1]-X[n-1]*X[n-1]*X[n-2]*Y[n]+Y[n-1]*X[n-2]*X[n]*X[n]+X[n-1]*X[n-1]*X[n]*Y[n-2]-Y[n-2]*X[n-1]*X[n]*X[n];

		// std::cout <<" (("<<a1<<"*x+"<<a2<<")*x+"<<a3<<")/"<<d << std::endl;
		return ((a1*x+a2)*x+a3)/d ;
#else // __LINBOX_HAVE_LAPACK
		int m = min(n,5);
		dvector_t X1((index_t)m) ;
		dvector_t Y1((index_t)m) ;
		for (int i = 0 ; i < m ; ++i) X1[(index_t)i] = X[(index_t)(n-m+i)] ;
		for (int i = 0 ; i < m ; ++i) Y1[(index_t)i] = Y[(index_t)(n-m+i)] ;
		return fit_lapack3(X1,Y1,x);

#endif // __LINBOX_HAVE_LAPACK
	}

} // LinBox


//
//  TimeWatcher
//

namespace LinBox {
	TimeWatcher::TimeWatcher (dvector_t & pts, dvector_t & vals) :
		Points_(&pts)
		,Values_(&vals)
		,MaxRepet_(12)
		,MinRepet_(2)
		,MaxTime_(0.5)
		,AbortTime_(10)
		,aborted_(false)
	{
		linbox_check(vals.size() == pts.size());
	}

	TimeWatcher::TimeWatcher () :
		Points_(NULL)
		,Values_(NULL)
		,MaxRepet_(12)
		,MinRepet_(2)
		,MaxTime_(0.5)
		,AbortTime_(10)
		,aborted_(false)
	{
	}


	void TimeWatcher::init(dvector_t & pts, dvector_t & vals)
	{
		linbox_check(vals.size() == pts.size());
		Points_ = &pts ;
		Values_ = &vals ;
	}


	dvector_t & TimeWatcher::refX()
	{
		return *Points_ ;
	}

	dvector_t & TimeWatcher::refY()
	{
		return *Values_ ;
	}

	// we don't assume that t(0sec) = 0 unless nothing has been computed yet.
	double TimeWatcher::predict(double x)
	{
		index_t Current_ = size();
		if (Current_ == 0)
			return 0 ;
		if (Current_ ==1 )
			return refX()[Current_-1] ; // unknown. could be known if we suppose t(0)=0
		if (Current_ == 2 ) {
			return fit2(refX(),refY(),1,x);
		}
		return fit3(refX(),refY(), (int) Current_-1,x);
	}

	bool TimeWatcher::keepon( index_t & repet, double tim, bool usePrediction )
	{

		if (aborted_)
			return false ;

		if (usePrediction) {
			if (predict(tim) < AbortTime_)
				return true ;
			else{
				aborted_ = true ;
				return false;
			}
		}

		if (repet<MinRepet_ || (tim < MaxTime_ && repet < MaxRepet_) ) {
			++repet ;
			return true;
		}
		return false ;
	}

	index_t TimeWatcher::size() const
	{
		if (Points_ == NULL || Values_ == NULL) {
			linbox_check(Values_ == NULL && Points_ == NULL);
			return  0 ;
		}
		linbox_check(Points_->size() == Values_->size());
		return (index_t)Points_->size();
	}

	void TimeWatcher::clear()
	{
		Points_ = NULL ;
		Values_ = NULL ;
	}

} // LinBox

//
// DataSeries
//

namespace LinBox {

	DataSeries:: DataSeries() :
		PointLabels(0)
		, Points(0)
		, Times(0)
		, Values(0)
		, UID(0)
	{}

	DataSeries::~DataSeries() {}

#if 0
	void
	DataSeries::resize(const index_t & n)
	{
		linbox_check(n == Values.size()+1);
		PointLabels.resize(n);
		Times.resize(n);
		Points.resize(n);
		Values.resize(n);
		UID.resize(n);

		return;
	}
#endif

	index_t
	DataSeries::size() const
	{
		linbox_check(PointLabels.size() == Points.size())
		linbox_check(Times.size() == Points.size())
		linbox_check(Times.size() == Values.size())
		linbox_check(Times.size() == UID.size())

		return (index_t)Values.size();
	}

	void DataSeries::push_back(const std::string & nam, const double & val, const double & x , const double &y )
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

		UID.push_back("point_" + randomAlNum(8));
		return;
	}

} // LinBox

//
// PlotGraph
//

namespace LinBox {

	//!@todo use getUsingSeries in latex/html/csv/xml


		void PlotGraph::_randomName()
		{
			std::ostringstream unique_filename ;
			unique_filename << _filename_ << '_' << getDateTime("_") << '_' << randomAlNum(4);

			// std::cout << unique_filename.str() << std::endl;

			_printname_ = unique_filename.str() ;

		}

		const std::string & PlotGraph::getFileName()
		{
			if (_printname_.empty())
				_randomName();
			return _printname_;
		}

		void PlotGraph::mergeTwoSeries( svector_t & merge_points
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

				index_t k ;
			       bool ok = findKeyword(k, merge_points.begin(), merge_points.begin()+data_size,pts[i]);

				if ( ok ){
					merge_data[idx][k] = dat[i] ;
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

		void PlotGraph::mergeSeries()
		{
			_data_. selectFirstSeries();
			_merge_points_ = _data_.getCurrentSeries( Point::Labels() ) ;
			_merge_data_[0] = _data_.getCurrentSeries( Point::Values() ) ;

			// std::cout << "merge points " << _merge_points_ << std::endl;
			// std::cout << "merge data   " << _merge_data_ << std::endl;

			for (index_t i = 1 ; i < _data_.size() ; ++i) {
				_data_. selectNextSeries() ;
				// std::cout << "to be merged "  << i << " : "  << std::endl;
				// std::cout << "new points " << _data_.getCurrentSeriesPointLabel() << std::endl;
				// std::cout << "new data   " << _data_.getCurrentSeriesValues() << std::endl;

				mergeTwoSeries(_merge_points_,_merge_data_,
					       _data_. getCurrentSeries( Point::Labels() ), _data_. getCurrentSeries( Point::Values() ),i);

				// std::cout << "result : " << std::endl;
				// std::cout << "merge points " << _merge_points_ << std::endl;
				// std::cout << "merge data   " << _merge_data_ << std::endl;

			}

			return ;
		}

		void PlotGraph::print_csv()
		{
			char comma   = ',';
			char comment = '#';
			index_t nb_points = (index_t)_merge_points_.size() ;
			index_t nb_series = (index_t)_data_.size() ;

			std::string unique_filename  = getFileName();
			std::string DataFileName = unique_filename + ".csv" ;
			std::ofstream DF(DataFileName.c_str());

			/*  Data file to be plot */
			DF.precision(2);
			// metadata
			DF << comment << fortifyString("title") << comma << fortifyString(_style_.getRawTitle()) << std::endl;
			DF << comment << fortifyString("date") << fortifyString(getDateTime()) << std::endl;
			smatrix_t uname = getMachineInformation();
			for (size_t i = 0 ; i < uname[0].size() ; ++i)
				DF << comment << fortifyString(uname[0][i]) << comma << fortifyString(uname[1][i]) << std::endl ;

			// data
			DF << fortifyString(_style_.getRawTitle(1)) << comma ;
			for (index_t i = 0 ; i < nb_series ; ++i) {
				DF << _data_.getSerieLabel(i) ;
				if (i != nb_series -1)
					DF << comma ;
			}
			DF << std::endl;

			for (index_t j = 0 ; j < nb_points ; ++j) {
				DF << _merge_points_[j] << comma;
				for (index_t i = 0 ; i < nb_series ; ++i) {
					DF  << _merge_data_[i][j] ;
					if (i != nb_series -1)
						DF << comma ;

				}
				DF << std::endl;
			}

			std::cout << "csv data in " << DataFileName << std::endl;

		}

		void PlotGraph::print_dat()
		{
			print_gnuplot(true);
		}

		void PlotGraph::print_xml()
		{
#ifdef __LINBOX_HAVE_TINYXML2
			std::string unique_filename = getFileName();
			unique_filename += ".xml" ;

			_data_.save(unique_filename,_style_.getRawTitle(),_style_.getRawTitle(1),_style_.getRawTitle(2));

#else
			std::cout << "tinyxml2 is not installed, could not print" << std::endl;
#endif

			load(unique_filename);
			return ;
		}

		void PlotGraph::print_html()
		{

			std::string comment_in = "<!--";
			std::string comment_out = "-->";
			index_t nb_points = (index_t)_merge_points_.size() ;
			index_t nb_series = (index_t)_data_.size() ;

			std::string unique_filename  = getFileName();
			std::string DataFileName = unique_filename + ".html" ;
			std::ofstream DF(DataFileName.c_str());

			/*  Data file to be plot */
			DF.precision(2);
			// metadata
			DF << comment_in << ("date") << (getDateTime()) << std::endl;
			smatrix_t uname = getMachineInformation();
			for (size_t i = 0 ; i < uname[0].size() ; ++i)
				DF << (uname[0][i]) << " : " << (uname[1][i]) << std::endl ;
			DF << comment_out << std::endl ;

			// data
			DF << "<table border=\"1\">" << std::endl;
			DF << "<caption> " <<  (_style_.getRawTitle()) << " (data in " << _style_.getRawTitle(2) << ')'  << " </caption>" << std::endl;
			DF << "<tr> " << std::endl;
			DF << "<th> " << _style_.getRawTitle(1) << " </th>";
			for (index_t i = 0 ; i < nb_series ; ++i) {
				DF << "<th> " << unfortifyString(_data_.getSerieLabel(i)) << " </th>";
			}
			DF << " </tr>" << std::endl;

			for (index_t j = 0 ; j < nb_points ; ++j) {
				DF << "<tr> " << std::endl;
				DF << "<th> " << unfortifyString(_merge_points_[j]) << " </th>";
				for (index_t i = 0 ; i < nb_series ; ++i) {
					DF  << "<td> " << _merge_data_[i][j]  << " </td>";

				}
				DF << std::endl;
				DF << "</tr>" << std::endl;
			}

			DF << "</table>" << std::endl;
			std::cout << "html data in " << DataFileName << std::endl;


		}

		void PlotGraph::print_latex()
		{
			index_t nb_points = (index_t)_merge_points_.size();
			index_t nb_series = _data_.size();

			linbox_check(nb_points);
			linbox_check(nb_series);
			// srand(time(NULL));
			// std::ostringstream unique_filename  ;
			std::string unique_filename = getFileName();
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

		void PlotGraph::print_gnuplot(bool only_data)
		{
#ifndef __LINBOX_HAVE_GNUPLOT
			std::cout << "gnuplot is not available on your system. only the data will be printed" << std::endl;
#endif
			index_t nb_points = (index_t)_merge_points_.size() ;
			index_t nb_series = (index_t)_data_.size() ;

			std::string unique_filename  = getFileName();
			std::string DataFileName = unique_filename + ".dat" ;
			std::ofstream DF(DataFileName.c_str());

			/*  Data file to be plot */
			// DF.precision(_style_.getPrecision());
			DF.precision(2);

			char comment = '#' ;
			char comma   = ' ';
			// metadata
			DF << comment << ("title") << comma << (_style_.getRawTitle()) << std::endl;
			DF << comment << ("date") << (getDateTime()) << std::endl;
			smatrix_t uname = getMachineInformation();
			for (size_t i = 0 ; i < uname[0].size() ; ++i)
				DF << comment << (uname[0][i]) << comma << (uname[1][i]) << std::endl ;


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

			if (only_data)
				std::cout << "data in " << DataFileName << std::endl;

#ifdef __LINBOX_HAVE_GNUPLOT
			if (!only_data) {
				std::string PlotFileName = unique_filename + ".gp" ;
				std::ofstream PF(PlotFileName.c_str());

				/*  Ploting script */
				PF << "#" << _filename_                    << std::endl;
				PF << _style_.getTerm()                    << std::endl;
				PF << _style_.getOutput(unique_filename)   << std::endl;
				PF << _style_.getTitle()                   << std::endl;
				PF << _style_.getKeyPos()                  << std::endl;
				PF << _style_.getXtics()                   << std::endl;
				PF << _style_.getPlotType()                << std::endl;

				PF << getPlotCommand(DataFileName) << std::endl;

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
			}
#endif


			return;
		}


		void PlotGraph::setData( PlotData & data )
		{
			_data_ = data ;
		}

		PlotData & PlotGraph::refData( PlotData & data)
		{
			return data = _data_ ;
		}

		void PlotGraph::setStyle( PlotStyle & style )
		{
			_style_ = style ;
		}

		PlotStyle & PlotGraph::refStyle( PlotStyle & style)
		{
			return style = _style_ ;
		}


		// not implemented yet
		void PlotGraph::sortSeries() {}

		// not implemented yet
		void PlotGraph::unique() {}

		PlotGraph::PlotGraph( PlotData & data, PlotStyle & style ) :
			_data_(data)
			,_style_(style)
			,_filename_("")
			,_printname_("")
			,_merge_data_(data.size())
			,_merge_points_(data.getSeries(0).size())
		{
			srand((unsigned)time(NULL));
			mergeSeries();
		}

		void PlotGraph::setOutFilename( const std::string & filename )
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

		}

		const std::string & PlotGraph::getUsingSeries() // const
		{
			// mutable _style_ ?
			linbox_check(_merge_points_.size());
			if (_style_.getUsingSeries().empty()) {
				_style_.setUsingSeries(std::pair<index_t,index_t>((index_t)2,(index_t)_merge_data_.size()+1));
			}
			return _style_.getUsingSeries();
		}

		std::string PlotGraph::getPlotCommand(const std::string & File) //const
		{
			std::string PC = "#plot\nplot \'" + File + "\' ";
			PC += getUsingSeries() ;
			return PC ;
		}

		void PlotGraph::print( LINBOX_enum(Tag::Printer) pt ) {
			switch (pt) {
			case (Tag::Printer::xml):
				{
					print_xml();
					break;
				}
			case (Tag::Printer::csv) :
				{
					print_csv();
					break;
				}
			case (Tag::Printer::dat) :
				{
					print_dat();
					break;
				}
			case (Tag::Printer::gnuplot) :
				{
					print_gnuplot();
					break;
				}
			case (Tag::Printer::tex) :
				{
					print_latex();
					break;
				}
			case (Tag::Printer::html) :
				{
					print_html();
					break;
				}
			default :
				{
					throw LinBoxError("printer unknown");
				}

			}

			return ;
		}

		void PlotGraph::save()
		{
			return print_xml();
		}

		void PlotGraph::load(const std::string & filename)
		{
			return _data_.load(filename);
		}

} // LinBox

#endif // __LINBOX_benchmarks_benchmark_C

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
