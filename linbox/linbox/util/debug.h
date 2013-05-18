/* linbox/util/debug.h
 *
 * Copyright (C) 2001,2010 LinBox
 * Copyright (C) 2001 Bradford Hovinen
 * Copyright (C) 2010 LinBox
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by BB.
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
#ifdef NDEBUG // à la assert.
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

#define THIS_CODE_COMPILES_BUT_IS_NOT_TESTED \
     std::cout << "*** Warning *** " << std::endl << __func__ << " in " << __FILE__ << ':' << __LINE__ << " is not tested" << std::endl;

#define THIS_CODE_MAY_NOT_COMPILE_AND_IS_NOT_TESTED \
	throw(" *** Warning ***  this piece of code is not compiled by default and may not work")

namespace LinBox
{ /*  Preconditions,Error,Failure,NotImplementedYet */
	/*!  A precondition failed.
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

		NotImplementedYet( const char * why)
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "*** ERROR ***"  << std::endl;
			(*_errorStream) << " This function is not implemented yet " ;
			(*_errorStream)	<< " (" << why << ")" <<std::endl;
		}

		NotImplementedYet(const char * function,
				  const char* file,
				  int line,
				  const char * why="\0")
		{
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << " *** ERROR *** (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
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
			      const char* file="\0",
			      int line=-1,
			      const char * what="\0")
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
			    const char * function="\0",
			    const char* file="\0",
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


namespace LinBox
{ /* Exceptions. */

	/*! @defgroup exceptions Exceptions.
	 * @brief Exceptions in LinBox (proposal, example in \c algorithms/hermite.h).
	 * If the algorithms cannot return as expected, then an exception is
	 * thrown. Hopefully it is catched later on.  Execptions, when thrown,
	 * don't write to any stream except in debug mode.  However, they can
	 * explain what they are for with the
	 * <code>const char *what(void)</code> member.
	 *
	 * Any exception derives from the \c LinBox::Exception class.
	 */

	/*! This is the exception class in LinBox.
	 * Any LinBox exception can derive from it.
	 */
	class Exception {
	// public:
		// Exception() {};
	} ;

	/*! Algorithmic exception.
	 */
	class algoException : public Exception {
	// public:
		// algoException() {};
	};

	/*! Not implemented yet.
	 * This piece of code is not fully implemented.
	 */
	class NotImplementedYetException : public Exception {
	};

	/*! Something bad an unexpected happened.
	 */
	class IrrecuperableException : public Exception {
	};

	/*! The input is not as expected.
	 */
	class BadInputException : public Exception {
	};
}

#define CONC(a,b) a ## b

#define LINBOX_SILENT_EXCEPTION(name) \
	   throw CONC(name,Exception) ()

#ifndef DEBUG
#define LINBOX_RAISE_EXCEPTION(name,why) \
	   throw CONC(name,Exception) ()
#else
#define LINBOX_RAISE_EXCEPTION(name,why) \
   do { \
	   std::cerr << " *** EXCEPTION *** (at " << __func__ << " in " << __FILE__ << ':' <<  __LINE__ << ") " << std::endl; \
	   std::cerr << "                   " << why << std::endl; \
	   throw CONC(name,Exception) () ; \
   } while(0)
#endif



#if defined(LinBoxSrcOnly) or defined(LinBoxTestOnly)
// for all-source compilation
#include "linbox/util/debug.C"
#endif

#include <fflas-ffpack/utils/print-utils.h>

#endif // __LINBOX_util_debug_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

