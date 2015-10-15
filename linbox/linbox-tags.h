/* Copyright (C) 2010,2011,2012 LinBox
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

/*! @file linbox-tags.h
 * @ingroup algorithms
 * @brief Provides tags for various algorithms/solutions, à la \c FFLAS.
 * This is a subset of Fflas* enums.
 */

#ifndef __LINBOX_linbox_tags_H
#define __LINBOX_linbox_tags_H

#include <fflas-ffpack/fflas/fflas_enum.h>
#include "linbox/linbox-config.h"

#define LINBOX_enum(name) name

#if HAVE_CXX11
#define SCOPED(name) struct name : int32_t
#else
#define SCOPED(name) name 
#endif

namespace LinBox
{

	/*! Structure for tags.
	 * Tags are simple enums that set a choice in a routine.
	 * For instance, if the user wants a <i>right</i> nullspace,
	 * she will use a \c Tag::Right parameter.
	 *
	 * There it total compatiblity with \c FFLAS tags (cross link)
	 * For instance, in LinBox, it is similar to use \c Tag::Upper and
	 * <code>(Tag::Shape) FFLAS::FflasUpper</code>.
	 *
	 * @note Tags are not Methods.
	 */
	namespace Tag {
		//! Left/Right Tag
		enum SCOPED(Side)
			{
				Left  = FFLAS::FflasLeft, //!< Left
				Right = FFLAS::FflasRight  //!< Right
			};

					   //! (No)Transpose Tag
		enum SCOPED(Transpose)
		{

			NoTrans = FFLAS::FflasNoTrans,
			Trans   = FFLAS::FflasTrans
		};

		//! (Upp/Low)er Tag
		enum SCOPED(UpLo)
		{
			Upper = FFLAS::FflasUpper,
			Lower = FFLAS::FflasLower
		} ;

		using Shape = UpLo;

		//! (Non)Unit Diagonal Tag
		enum SCOPED(Diag)
		{
			NonUnit = FFLAS::FflasNonUnit,
			Unit    = FFLAS::FflasUnit
		} ;

		//! Dense format (table) output Tag
		enum SCOPED(FileFormat)
		{
			Plain = 0,
			Maple = 1,
			HTML  = 2,
			LaTeX = 3,
			Detect,
			Guillaume,
			Turner,
			Matlab,
			Pretty,
			MagmaCpt,
			OneBased,
			MatrixMarket

		} ;



		enum SCOPED(Direction)
		{
			Row = 10 ,
			Col = 11
		} ;

		enum SCOPED(Sign)
		{
			Positive    = 500 , //! is  >0 (strictement positif)
			Negative    = 501 , //! is  <0 (strictement négatif)
			NonNegative = 502 , //! is >=0 (positif)
			NonPositive = 503 , //! is <=0 (négatif)
			NonZero     = 504 , //! is !=0 (non nul)
			Zero        = 505   //! is ==0 (nul)

		};

	} // namespace Tag

} // namespace LinBox
#undef SCOPE
#endif // __LINBOX_linbox_tags_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
