/* Copyright (C) 2010,2011,2012 LinBox
 * Written by <brice.boyer@imag.fr>
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

/*! @file linbox-tags.h
 * @ingroup algorithms
 * @brief Provides tags for various algorithms/solutions, à la \c FFLAS.
 * This is a subset of Fflas* enums.
 */

#ifndef __LINBOX_linbox_tags_H
#define __LINBOX_linbox_tags_H

#include <fflas-ffpack/fflas/fflas.h>

#if HAVE_CXX11
#define LINBOX_enum(name) name
#else
#define LINBOX_enum(name) name::enum_t
#endif

namespace LinBox
{

	/*! Structure for tags.
	 * Tags are simple enums that set a choice in a routine.
	 * For instance, if the user wants a <i>right</i> nullspace,
	 * she will use a \c LinBoxTag::Right parameter.
	 *
	 * There it total compatiblity with \c FFLAS tags (cross link)
	 * For instance, in LinBox, it is similar to use \c LinBoxTag::Upper and
	 * <code>(LinBoxTag::Shape) FFLAS::FflasUpper</code>.
	 *
	 * @note Tags are not Methods.
	 */
	namespace LinBoxTag {
		//! Left/Right Tag
#if HAVE_CXX11
		enum struct Side : int32_t
#else
		struct Side { enum enum_t
#endif
					   {
						   Left  = FFLAS::FflasLeft, //!< Left
						   Right = FFLAS::FflasRight  //!< Right
					   };

					   //! (No)Transpose Tag
#if HAVE_CXX11

		enum struct Transpose : int32_t
#else
		};

		struct Transpose { enum enum_t
#endif
		{

			NoTrans = FFLAS::FflasNoTrans,
			Trans   = FFLAS::FflasTrans
		};

		//! (Upp/Low)er Tag
#if HAVE_CXX11
		enum struct Shape : int32_t
#else
		};

		struct Shape { enum  enum_t
#endif
		{
			Upper = FFLAS::FflasUpper,
			Lower = FFLAS::FflasLower
		} ;

		//! (Non)Unit Diagonal Tag
#if HAVE_CXX11
		enum struct Diag : int32_t
#else
		};

		struct Diag { enum enum_t
#endif
		{
			NonUnit = FFLAS::FflasNonUnit,
			Unit    = FFLAS::FflasUnit
		} ;

		//! Dense format (table) output Tag
#if HAVE_CXX11
		enum  class Format : int32_t
#else
		};

		struct Format { enum  enum_t
#endif
		{
			FormatPlain = 0,
			FormatMaple = 1,
			FormatHTML  = 2,
			FormatLaTeX = 3
		} ;

#if HAVE_CXX11
		enum struct Direction : int32_t
#else
		} ;

		struct Direction {  enum enum_t
#endif
		{
			Row = 10 ,
			Col = 11
		} ;

#if HAVE_CXX11
		enum struct Sign : int32_t
#else
		};

		struct Sign { enum enum_t
#endif
		{
			Positive    = 500 , //! is  >0 (strictement positif)
			Negative    = 501 , //! is  <0 (strictement négatif)
			NonNegative = 502 , //! is >=0 (positif)
			NonPositive = 503 , //! is <=0 (négatif)
			NonZero     = 504 , //! is !=0 (non nul)
			Zero        = 505   //! is ==0 (nul)

		};

#if HAVE_CXX11
#else
		};
#endif
	}

}

#endif // __LINBOX_linbox_tags_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
