/* File: src/wrappers/by_library/ntl.h
 * Author: William J Turner for the LinBox group
 */

#ifndef _NTL_
#define _NTL_

#include "unparam_field.h"

/** @name NTL Fields.
 * Partial template instantiantions for wrapping objects from 
 * \URL[Victor Shoup]{http://www.shoup.net/index.html}'s
 * package \URL[NTL]{http://www.shoup.net/ntl/index.html} 5.0b
 * objects to comply with the common object interface for 
 * \Ref{LinBox} fields.
 *
 * See \URL{http://www.shoup.net/ntl/doc/tour.html} for a tour of the
 * \URL[NTL]{http://www.shoup.net/ntl/index.html} library. 
 */
//@{

/** @name class RR: arbitrary precision floating point numbers
 */
//@{

#include <NTL/RR.h>

// zero equality
template <> bool unparam_field<NTL::RR>::isZero(const NTL::RR& x) const
{ return static_cast<bool>(IsZero(x)); }

// one equality
template <> bool unparam_field<NTL::RR>::isOne(const NTL::RR& x) const
{ return static_cast<bool>(IsOne(x)); }




//@} // class RR







//@} NTL Fields

#endif // _NTL_
