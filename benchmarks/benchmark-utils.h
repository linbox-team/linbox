/* Copyright (C) 2013 LinBox
 * Written by Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file   benchmarks/benchmark-utils.h
 * @ingroup benchmarks
 * @brief   utils
  */

#ifndef __LINBOX_benchmarks_benchmark_utils_H
#define __LINBOX_benchmarks_benchmark_utils_H


/* ********************** */
/*        Tags            */
/* ********************** */

// Tags (cf linbox/linbox-tags.h )
namespace LinBox { namespace Tag {

	//! selection of best times in a series.
#if HAVE_CXX11
	enum struct TimeSelect : int32_t
#else
				 enum TimeSelect 
#endif
					 {
						 average   = 1, //!< select the average
						 bestThree = 2, //!< select the average among the three best (supposes order is <)
						 bestOne   = 3, //!< select the best one (supposes order is <)
						 median    = 4, //!< median
						 medmean   = 5  //!< interquartile mean (remove 25% extremes and average the rest)
					 };

	//! selection of printers
#if HAVE_CXX11
	enum struct Printer : int32_t
#else
			      enum Printer 
#endif
				      {
					      dat       = 1, //!< print data  raw (format chosen by implementation)
					      tex       = 2, //!< print latex useable
					      xml       = 3, //!< print using xml syntax
					      gnuplot   = 4, //!< print using gnuplot
					      csv       = 5, //!< coma separated values
					      html      = 6  //!< coma separated values
				      };

	} // Tag
} // LinBox

/* ********************** */
/*        Outils          */
/* ********************** */

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
// String processing
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

	//! finds keyword betwen begin and end, return true if found and i is the index where it is (possibly correspondig to end)
	bool findKeyword(size_t & i, const svector_t::const_iterator & begin, const svector_t::const_iterator & end, const std::string & keyword)
		{
			svector_t::const_iterator it ;
			it = std::find(begin, end, keyword);
			i = (size_t)std::distance(begin, it);
			return (it != end) ;

		}

	/*! @internal
	 * @brief random <code>:alnum:</code> \c char.
	 * [[:alnum:]] characters are in range
	 *  - num : 48-57
	 *  - AL  : 67-90
	 *  - al  : 97-122
	 *  .
	 * @return a random alphabetic or numeric char.
	 */
	char randomAlNum() ;

	/*! @internal
	 * @brief random <code>:alnum:</code> \c string.
	 * @param m size of string
	 * [[:alnum:]] characters are in range
	 *  - num : 48-57
	 *  - AL  : 67-90
	 *  - al  : 97-122
	 *  .
	 * @return a random alphabetic or numeric char.
	 */
	std::string randomAlNum(const size_t & m) ;
}// LinBox

//
// Machine info
//

namespace LinBox {

	/*! get ISO time and date
	 * year-time YYYY-MM-DD 'sep' HH:MM:SS 'sep' (GMT)
	 * @param sep separation between
	 */
	std::string getDateTime( const std::string & sep = " ");

	//! get some machine information (not cpu yet)
	smatrix_t getMachineInformation();

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
	void showAdvanceLinear(size_t curr, size_t min, size_t max);

	/*! tells the current series of measure has completed (on the terminal)
	 * @param curr current iteration
	 * @param all number of iterations
	 */
	void showFinish(size_t curr, size_t all);

	/*! tells the current series of measure was skipped (on the terminal)
	 * @param curr current iteration
	 * @param all number of iterations
	 */
	void showSkip(size_t curr, size_t all);

	//! Show progression on the terminal (helper)
	class showProgression {
	private :
		size_t _cur_ ; //!< current iter
		size_t _tot_ ; //!< max iter
	public :
		//! constructor
		showProgression (size_t tot) ;

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
	double computeMFLOPS(const double & tim, const double mflo, const size_t rpt = 1);

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
	double computeMFLOPS(const dvector_t & tim, const double mflo, Tag::TimeSelect ts = Tag::TimeSelect::bestThree);

}

//
// Chrono
// warning : implem in benchmark.inl
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


#ifdef LinBoxSrcOnly
#include "benchmarks/benchmark-utils.C"
#endif

#endif // __LINBOX_benchmarks_benchmark_utils_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
