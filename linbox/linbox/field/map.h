/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2008 LinBox
 * Written by JGD
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


#ifndef __LINBOX_map_H
#define __LINBOX_map_H

#include "linbox/field/hom.h"

namespace LinBox
{

	template< class Source, class Target >
	struct Map {

		typedef Hom<Source,Target> Homomorphism;

		Map(const Homomorphism& H) :
			_hom(H)
		{}
		Map(const Source& S, const Target& T) :
			_hom(S,T)
		{}


		template< typename ContainerTargetElement,
		typename ContainerSourceElement >
		ContainerTargetElement& operator() (
						    ContainerTargetElement& tgt,
						    const ContainerSourceElement& src)
		{

			tgt.resize(src.size());
			typename ContainerTargetElement::iterator tgt_it(tgt.begin());
			typename ContainerSourceElement::const_iterator src_it(src.begin());

			for( ; src_it != src.end(); ++src_it, ++tgt_it)
				this->_hom.image (*tgt_it, *src_it);

			return tgt;
		}


	private:
		Homomorphism _hom;
	};


	template< class Source, class Target >
	struct PreMap {

		typedef Hom<Source,Target> Homomorphism;

		PreMap(const Homomorphism& H) :
			_hom(H)
		{}
		PreMap(const Source& S, const Target& T) :
			_hom(S,T)
		{}

		template< typename ContainerSourceElement,
		typename ContainerTargetElement >
		ContainerSourceElement& operator() (
						    ContainerSourceElement& src,
						    ContainerTargetElement tgt)
		{

			src.resize(tgt.size());
			typename ContainerSourceElement::iterator src_it(src.begin());
			typename ContainerTargetElement::const_iterator tgt_it(tgt.begin());

			for( ; tgt_it != tgt.end(); ++src_it, ++tgt_it) {
				this->_hom.preimage (*src_it, *tgt_it);
			}

			return src;
		}


	private:
		Homomorphism _hom;
	};


}

#endif // __LINBOX_map_H

