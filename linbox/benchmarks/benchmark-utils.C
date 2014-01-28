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
 * @file   benchmarks/benchmark-utils.C
 * @ingroup benchmarks
 * @brief   utils
  */

#ifndef __LINBOX_benchmarks_benchmark_utils_C
#define __LINBOX_benchmarks_benchmark_utils_C

#include "linbox/linbox-config.h"
#include "benchmark-utils.h"

//
// Terminal progression
//

namespace LinBox {

	void showAdvanceLinear(index_t curr, index_t min, index_t max)
	{
		std::cout << std::setprecision(4) << "\033[2K" << "\033[30D" << min <<std::flush;
		std::cout << '<' << curr << '<' << max << " (" << std::flush;
		std::cout << double(curr-min)/double(max-min)*100 << "%)" << std::flush;
	}

	void showFinish(index_t curr, index_t all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "finished : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}

	void showSkip(index_t curr, index_t all)
	{
		std::cout <<  "\033[2K" << "\033[30D" << "skipped : " << curr << std::flush;
		std::cout << '/' << all-1 << std::flush << std::endl;
	}

	showProgression::showProgression (index_t tot) :
		_cur_(0)
		,_tot_(tot)
	{}

	//! show an inter has finished
	void showProgression::FinishIter()
	{
		++_cur_ ;
		showFinish(_cur_,_tot_);
	}

	//! show an inter has been skipped.
	void showProgression::SkipIter()
	{
		++_cur_ ;
		showSkip(_cur_,_tot_);
	}


} // LinBox

//
// MFLOPS helper
//

namespace LinBox {

	double computeMFLOPS(const double & tim, const double mflo, const index_t rpt)
	{
		linbox_check(rpt);
		// linbox_check(tim != 0.);
		if (tim == 0.) return NAN ;
		return (double) ((mflo*rpt)/tim);
	}

	dvector_t &
	insertTime(dvector_t & tim3, const double & tps)
	{
		linbox_check(tim3.size() == 3);
		if (tim3[0] > tps) {
			tim3[2] = tim3[1] ;
			tim3[1] = tim3[0] ;
			tim3[0] = tps ;
		}
		else if (tim3[1] > tps) {
			tim3[2] = tim3[1] ;
			tim3[1] = tps ;
		}
		else if (tim3[2] > tps) {
			tim3[2] = tps ;
		}
		return tim3 ;
	}

	double computeMFLOPS(const dvector_t & tim, const double mflo, LINBOX_enum(Tag::TimeSelect) ts )
	{
		linbox_check(tim.size());
		switch (ts) {
		case (Tag::TimeSelect::average) :
			{
				double tps = 0 ;
				for (size_t i = 0 ; i < tim.size() ; ++i)
					tps += tim[i] ;
				return computeMFLOPS(tps,mflo,(index_t)tim.size());
			}
		case (Tag::TimeSelect::bestThree) :
			{
				if (tim.size() <4)
					return computeMFLOPS(tim,mflo,Tag::TimeSelect::average);

				dvector_t tps (3);
				double t1,t2 ;
				if (tim[0]<tim[1]) {
					t1 = tim[0];
					t2 = tim[1];
				}
				else {
					t1 = tim[1];
					t2 = tim[0];
				}
				if (tim[3] < t1) {
					tps[0] = tim[3] ;
					tps[1] = t1 ;
					tps[2] = t2 ;
				}
				else if (tim[2] < t1) {
					tps[0] = t1 ;
					tps[1] = tim[3] ;
					tps[2] = t2 ;
				}
				else {
					tps[0] = t1 ;
					tps[1] = t2;
					tps[2] = tim[3] ;
				}

				for (size_t i = 3 ; i < tim.size() ; ++i)
					insertTime(tps,tim[i]);

				return computeMFLOPS(tim,mflo,Tag::TimeSelect::average);

			}
		case (Tag::TimeSelect::bestOne) :
			{
				double t1 = tim[0] ;
				for (size_t i = 1 ; i < tim.size() ; ++i)
					if (tim[i] < t1)
						t1 = tim[i] ;
				return computeMFLOPS(t1,mflo,1);

			}
		case (Tag::TimeSelect::median) :
			{
				if (tim.size() == 1)
					return computeMFLOPS(tim[0],mflo,1) ;

				dvector_t tps (tim);
				std::sort(tps.begin(),tps.end());
				index_t mid = (index_t)tps.size()/2 ;
				double t1 ;
				if (Givaro::isOdd(tps.size()))
					t1 = tps[mid] ;
				else
					t1 = (tps[mid-1]+tps[mid])/2;
				return computeMFLOPS(t1,mflo,1);
			}
		case (Tag::TimeSelect::medmean) :
			{
				if (tim.size() < 3)
					return computeMFLOPS(tim,mflo,Tag::TimeSelect::median); ;

				index_t q1 = (index_t)std::round((double)tim.size()/(double)4) ;
				index_t q3 = (index_t)tim.size()-q1 ;
				dvector_t tps (tim);
				std::sort(tps.begin(),tps.end());
				dvector_t tps2 (tim.begin()+q1,tim.begin()+q3);
				return computeMFLOPS(tps2,mflo,Tag::TimeSelect::average);
			}

		default :
			{
				throw("not among tags");
			}

		} // switch(ts)
	}

}

//
// String processing
//

namespace LinBox {

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

	char randomAlNum()
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


	std::string randomAlNum(const size_t & m)
	{
		std::string r = "" ;
		for (size_t i = 0 ; i < m ; ++i)
			r += randomAlNum();
		return r ;
	}

}// LinBox

//
// Machine information
//

namespace LinBox {

	//! get ISO time and date
	std::string getDateTime(const std::string & sep)
	{
		std::time_t rawtime;
		std::tm* timeinfo;
		char buffer [80];

		std::time(&rawtime);
		timeinfo = std::gmtime(&rawtime);

		std::string fmt ;
		std::string date = "%Y-%m-%d" ;
		std::string time = "%H:%M:%S" ;
		std::string tz   = "GMT" ;

		fmt = date + sep + time + sep + tz ;

		std::strftime(buffer,80,fmt.c_str(),timeinfo);

		std::string mytime(buffer);

		return mytime;
	}

	//! get some machine information (not cpu yet)
	smatrix_t getMachineInformation()
	{
		smatrix_t Machine(2);
		Machine[0].resize(5);
		Machine[1].resize(5);
		struct utsname unameData;
		uname(&unameData);
		Machine[0][0] = "sysname";
		Machine[1][0] = unameData.sysname;
		Machine[0][1] = "nodename";
		Machine[1][1] = unameData.nodename;
		Machine[0][2] = "release";
		Machine[1][2] = unameData.release;
		Machine[0][3] = "version";
		Machine[1][3] = unameData.version;
		Machine[0][4] = "machine";
		Machine[1][4] = unameData.machine;
		// Machine[0][5] = "RAM (kb)";
		// system("cat /proc/meminfo  | grep MemTotal | awk '{print $2}'");
		// Machine[0][6] = "CPU name";
		// system("/proc/cpuinfo  | grep 'model name' | awk '{$1=$2=$3=""; print $0}'");
		// Machine[0][7] = "CPU nb";
		// system("/proc/cpuinfo  | grep 'model name' | wc -l");
		// Machine[0][8] = "CPU Ghz";
		// cpuid ?
		return Machine ;
	}

}

#endif // __LINBOX_benchmarks_benchmark_utils_C

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
