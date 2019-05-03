/* linbox/util/error.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
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

#ifndef __LINBOX_util_error_H
#define __LINBOX_util_error_H

#include <cstring>
#include <iostream>

/**
 * I feel like we need:
 *  CheckFailed     (old PreconditionFailed, only throw by linbox_check()
 *  NotImplementedYet
 *  NotTested       (old THIS_CODE_COMPILES_BUT_IS_NOT_TESTED)
 *      ?? (NO PROBABLY NOT AN EXCEPTION)
 *  InternalError   (old LinBoxFailure)
 *  Error           (old LinBoxError/LinboxError/lb_runtime_error)
 *      ?? (OR RuntimeError)
 *  InconsistentSystem (old InconsistentSystem/MathInconsistentSystem)
 *      ?? (SHOULD WE KEEP CERTIFICATE INSIDE? SHOULD IT BE A CHILD OF SOLVEFAILED?)
 *  SolveFailed
 *  InvalidInput    (old InvalidMatrixInput)
 *
 * Also feels like we need a auto __LINE__ __FILE__,
 * so Error(msg) might be a macro expanding to ErrorException(msg, __LINE__, __FILE__) being a class.
 */

namespace LinBox { /*  Preconditions,Error,Failure,NotImplementedYet */
    /*!  A precondition failed.
     * @ingroup util
     * The \c throw mechanism is usually used here as in
     \code
     if (!check)
     throw(PreconditionFailed(__func__,__LINE__,"this check just failed");
     \endcode
     * The parameters of the constructor help debugging.
     */
    class PreconditionFailed { //: public LinboxError BB: otherwise,  error.h:39 segfaults
        static std::ostream* _errorStream;

    public:
        /*! @internal
         * A precondtion failed.
         * @param function usually \c __func__, the function that threw the error
         * @param line     usually \c __LINE__, the line where it happened
         * @param check    a string telling what failed.
         */
        PreconditionFailed(const char* function, int line, const char* check)
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

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
        PreconditionFailed(const char* function, const char* file, int line, const char* check)
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << "ERROR (at " << function << " in " << file << ':' << line << "): " << std::endl;
            (*_errorStream) << "Precondition not met:" << check << std::endl;
        }

        static void setErrorStream(std::ostream& stream);

        /*! @internal overload the virtual print of LinboxError.
         * @param o output stream
         */
        std::ostream& print(std::ostream& o) const
        {
            if (std::ostringstream* str = dynamic_cast<std::ostringstream*>(_errorStream))
                return o << str->str();
            else
                throw LinboxError("LinBox ERROR: PreconditionFailed exception is not initialized correctly");
        }
    };

    /*! @internal A function is "not implemented yet(tm)".
     * where, why ?
     */
    class NotImplementedYet {
    protected:
        static std::ostream* _errorStream;

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

        NotImplementedYet(const std::string& why)
            : NotImplementedYet(why.c_str())
        {
        }

        NotImplementedYet(const char* why)
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << "*** ERROR ***" << std::endl;
            (*_errorStream) << " This function is not implemented yet ";
            (*_errorStream) << " (" << why << ")" << std::endl;
        }

        NotImplementedYet(const char* function, const char* file, int line, const char* why = "\0")
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << " *** ERROR *** (at " << function << " in " << file << ':' << line << "): " << std::endl;
            (*_errorStream) << " This function is not implemented yet";
            if (why)
                (*_errorStream) << " (" << why << ")" << std::endl;
            else
                (*_errorStream) << "." << std::endl;
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
        LinBoxFailure(const char* function, const char* file = "\0", int line = -1, const char* what = "\0")
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << "ERROR (at " << function << " in " << file << ':' << line << "): " << std::endl;
            (*_errorStream) << " failure : ";
            if (what)
                (*_errorStream) << what << "." << std::endl;
            else
                (*_errorStream) << "no explanation." << std::endl;
        }
    };

    /*! @internal Something is wrong.
     * what ?
     */
    class LinBoxError : public NotImplementedYet {
    public:
        LinBoxError(const std::string& what, const char* function = "\0", const char* file = "\0", int line = -1)
            : LinBoxError(what.c_str(), function, file, line)
        {
        }

        /*! @internal
         * User failed.
         * The parameter help debugging/explaining.
         * @param function usually \c __func__, the function that threw the error
         * @param file     usually \c __FILE__, the file where this function is
         * @param line     usually \c __LINE__, the line where it happened
         * @param what     what happened ? should not be NULL...
         */
        LinBoxError(const char* what, const char* function = "\0", const char* file = "\0", int line = -1)
        {
            if (_errorStream == (std::ostream*)0) _errorStream = &std::cerr;

            (*_errorStream) << std::endl << std::endl;
            (*_errorStream) << " *** ERROR *** : " << what << std::endl;
            if (function) {
                (*_errorStream) << "(at " << function;
                if (file) {
                    (*_errorStream) << " in " << file;
                    if (line >= 0) (*_errorStream) << ':' << line;
                    (*_errorStream) << ")";
                }
            }
            else if (file) { // extremely unlikely...
                (*_errorStream) << "(in " << file;
                if (line >= 0) (*_errorStream) << ':' << line;
                (*_errorStream) << ")";
            }
            (*_errorStream) << std::endl;
        }
    };
}

namespace LinBox {

    // ------------------------------- LinboxError
    /** base class for execption handling in LinBox
      \ingroup util
      */
    class LinboxError {
        static const size_t max_error_string = 256;

    public:
        LinboxError(const std::string& msg)
            : LinboxError(msg.c_str())
        {
        }

        LinboxError(const char* msg = "\0")
        {
            std::strncpy(strg, msg, max_error_string);
            strg[max_error_string - 1] = 0;
        };

        const char* what() const { return strg; }

        // -- virtual print of the error message
        virtual std::ostream& print(std::ostream& o) const { return o << strg << std::endl; }

        // -- non virtual output operator
        friend std::ostream& operator<<(std::ostream& o, const LinboxError& E);

        // - useful to setup a break point on it
        static void throw_error(const LinboxError& err) { throw err; }

        virtual ~LinboxError() {}

    protected:
        char strg[max_error_string];
    };

    class LinboxMathError : public LinboxError {
    public:
        LinboxMathError(const char* msg)
            : LinboxError(msg){};
    };

    class LinboxMathInconsistentSystem : public LinboxMathError {
    public:
        LinboxMathInconsistentSystem(const char* msg)
            : LinboxMathError(msg){};
    };

    // -- Exception thrown when probabilistic solve fails
    class SolveFailed : public LinboxError {
    public:
        SolveFailed(const char* msg)
            : LinboxError(msg){};
        SolveFailed()
            : LinboxError(){};
    };

    /**
     * Exception thrown when the system to be solved is
     * inconsistent. Contains a certificate of inconsistency.
     */
    template <class Vector>
    class InconsistentSystem {
    public:
        InconsistentSystem(const Vector& u)
            : _u(u)
        {
        }

        const Vector& certificate() const { return _u; }

    private:
        Vector _u;
    };
}

/** Exception class for invalid matrix input
 */
namespace Exceptions {
    class InvalidMatrixInput {
    };
}

#ifdef LinBoxSrcOnly // for all-source compilation
#include "linbox/util/error.C"
#endif

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
