/* -*- mode: c; style: linux -*- */

/* linbox/src/util/debug.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __DEBUG_H
#define __DEBUG_H

#include <iostream>

#ifdef NDEBUG
#  define linbox_check(check)
#else
#  define linbox_check(check) \
        if (!(check)) \
                 throw LinBox::PreconditionFailed (__FUNCTION__, __LINE__, #check);
#endif

// Namespace in which all LinBox library code resides
namespace LinBox
{

	class PreconditionFailed
	{
		static ostream *_errorStream;

	    public:
		PreconditionFailed (char *function, int line, char *check) {
			if (_errorStream == (ostream *) 0)
				_errorStream = &cerr;

			(*_errorStream) << "ERROR (" << function << ":" << line << "): ";
			(*_errorStream) << "Precondition " << check << " not met" << endl;
		}

		static void setErrorStream (ostream &stream);
	};

	void PreconditionFailed::setErrorStream (ostream &stream)
	{
		_errorStream = &stream;
	}

	ostream *PreconditionFailed::_errorStream;

} // namespace LinBox

#endif // __DEBUG_H
