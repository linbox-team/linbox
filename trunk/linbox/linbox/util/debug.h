/* linbox/util/debug.h
 *
 * Copyright (C) 2001,2010 LinBox
 * Copyright (C) 2001 Bradford Hovinen
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
 *
 */

#ifndef __LINBOX_util_debug_H
#define __LINBOX_util_debug_H

#include <iostream>
#include <sstream>
#include "linbox/util/error.h"

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
	class PreconditionFailed //: public LinboxError BB: otherwise,  error.h:39 segfaults
	{
		static std::ostream *_errorStream;

	    public:
		PreconditionFailed (const char *function, int line, const char *check) { //BB : not so useful, is it ?
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (" << function << ":" << line << "): ";
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		PreconditionFailed (const char* function, const char *file, int line, const char *check) {
			if (_errorStream == (std::ostream *) 0)
				_errorStream = &std::cerr;

			(*_errorStream) << std::endl << std::endl;
			(*_errorStream) << "ERROR (at " << function << " in " << file << ':' <<  line << "): " << std::endl;
			(*_errorStream) << "Precondition not met:" << check << std::endl;
		}

		static void setErrorStream (std::ostream &stream);

		// -- overload the virtual print of LinboxError
		std::ostream &print (std::ostream &o) const { 
			if (std::ostringstream * str = dynamic_cast<std::ostringstream*>(_errorStream))
				return o << str->str() ; 
			else
				throw LinboxError("LinBox ERROR: PreconditionFailed exception is not initialized correctly");
		}
	};
}

#ifdef LinBoxSrcOnly
// for all-source compilation
#include <linbox/util/debug.C>
#endif

#endif // __LINBOX_util_debug_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
