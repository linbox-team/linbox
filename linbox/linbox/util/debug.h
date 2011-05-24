/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/util/debug.h
 *
 * Copyright (C) 2001,2010 LinBox
 * Copyright (C) 2001 Bradford Hovinen
 * Copyright (C) 2010 LinBox
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by BB.
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 */

/*! @file util/debug.h
 * @ingroup util
 * Various utilities for debugging.
 * @todo we should put vector printing elsewhere.
 */

#ifndef __LINBOX_util_debug_H
#define __LINBOX_util_debug_H

#include <iostream>
#include <sstream>
#include "linbox/util/error.h"
#include <vector>
#include <list>

/*! Check an assertion (à la \c std::assert).
 * If in DEBUG mode, throws a \ref PreconditionFailed exception.
 * In REALEASE mode, nothing is checked.
 * @param check assertion to be checked.
 */
#ifndef DEBUG
#  define linbox_check(check)
#else
#  ifdef __GNUC__
#    define linbox_check(check) \
        if (!(check)) \
                 throw LinBox::PreconditionFailed (__func__, __FILE__, __LINE__, #check); //BB : should work on non gnu compilers too
#  else
#    define linbox_check(check) \
        if (!(check)) \
                 throw LinBox::PreconditionFailed (__FILE__, __LINE__, #check);
#  endif
#endif

namespace LinBox
{
	/*!  A precondtion failed.
	 * @ingroup util
	 * The \c throw mechanism is usually used here as in
	 \code
	 if (!check)
	 throw(PreconditionFailed(__func__,__LINE__,"this check just failed");
	 \endcode
	 * The parameters of the constructor help debugging.
	 */
	class PreconditionFailed {//: public LinboxError BB: otherwise,  error.h:39 segfaults
		static std::ostream *_errorStream;

	public:
		/*! @internal
		 * A precondtion failed.
		 * @param function usually \c __func__, the function that threw the error
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param check    a string telling what failed.
		 */
		PreconditionFailed (const char *function, int line, const char *check)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (" << function << ":" << line << "): ";
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		/*! @internal
		 * A precondtion failed.
		 * The parameter help debugging. This is not much different from the previous
		 * except we can digg faster in the file where the exception was triggered.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param check    a string telling what failed.
		 */
		PreconditionFailed (const char* function, const char *file, int line, const char *check)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		static void setErrorStream (std::ostream &stream);

		/*! @internal overload the virtual print of LinboxError.
		 * @param o output stream
		 */
		std::ostream &print (std::ostream &o) const
		{
			if (std::ostringstream * str = dynamic_cast<std::ostringstream*>(_errorStream))
				return o << str->str() ;
			else
				throw LinboxError("LinBox ERROR: PreconditionFailed exception is not initialized correctly");
		}
	};

	/*! @internal A function is "not implemented yet(tm)".
	 * where, why ?
	 */
	class NotImplementedYet {
	protected:
		static std::ostream *_errorStream;

	public:
		/*! @internal
		 * A precondtion failed.
		 * The parameter help debugging. This is not much different from the previous
		 * except we can digg faster in the file where the exception was triggered.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param why      by default, lazy people don't provide an explanation.
		 */
		NotImplementedYet() {}

		NotImplementedYet(const char * function,
				  const char* file,
				  int line,
				  const char * why='\0')
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << " This function is not implemented yet" ;
			if (why)
				(*_errorStream)	<< " (" << why << ")" <<std::endl;
			else
				(*_errorStream)	<<  "." << std::endl;

		}
	};

	/*! @internal Something went wrong.
	 * what ?
	 */
	class LinBoxFailure : public NotImplementedYet {
	public:
		/*! @internal
		 * LinBox failed.
		 * The parameter help debugging/explaining.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param what     what happened ? should not be NULL...
		 */
		LinBoxFailure(const char * function,
			      const char* file,
			      int line,
			      const char * what='\0')
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << " failure : " ;
			if (what)
				(*_errorStream)	<< what << "." <<std::endl;
			else
				(*_errorStream)	<<  "no explanation." << std::endl;

		}
	};

	/*! @internal Something is wrong.
	 * what ?
	 */
	class LinBoxError : public NotImplementedYet {
	public:
		/*! @internal
		 * User failed.
		 * The parameter help debugging/explaining.
		 * @param function usually \c __func__, the function that threw the error
		 * @param file     usually \c __FILE__, the file where this function is
		 * @param line     usually \c __LINE__, the line where it happened
		 * @param what     what happened ? should not be NULL...
		 */
		LinBoxError(const char * what,
			    const char * function='\0',
			    const char* file='\0',
			    int line=-1)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << " *** ERROR *** : " << what << std::endl;
			if (function) {
				(*_errorStream) << "(at " << function ;
				if (file) {
					(*_errorStream) << " in " << file ;
					if (line>=0)
						(*_errorStream) << ':' <<  line ;
					(*_errorStream) << ")" ;
				}
			}
			else if (file) { // extremely unlikely...
				(*_errorStream) << "(in " << file ;
				if (line>=0)
					(*_errorStream) << ':' <<  line ;
				(*_errorStream) << ")" ;
			}
			(*_errorStream) << std::endl;
		}
	};
}

#ifdef LinBoxSrcOnly
// for all-source compilation
#include <linbox/util/debug.C>
#endif

#include <fflas-ffpack/utils/print-utils.h>

#endif // __LINBOX_util_debug_H

