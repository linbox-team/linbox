/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/blackbox/blackbox-interface.h
 * Copyright (C) 2002 LinBox
 * Written by David Saunders
 *
 *  ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#ifndef __LINBOX_blackbox_interface_H
#define __LINBOX_blackbox_interface_H

#include "linbox/element/archetype.h"

namespace LinBox
{
	// LinBox Blackbox Interface
	/*
	 * The LinBox @link BlackboxInterface@endlink common object @link Interfaces interface@endlink.
	 * The blackbox interface includes the public members defined in the archetype.
	 */
	/**
	 * \brief This blackbox base class exists solely to aid documentation organization.

	 *  For the general blackbox member function documentation consult the @link blackbox Archetype@endlink. For specific properties of individual representations consult the specific blackbox classes.
	 */
	class BlackboxInterface {
		/*
		   public:
		// this just demo's that some declarations could be here.
		typedef ElementArchetype Element;
		virtual Element& mul(Element& c, const Element& a, const Element& b) const = 0;
		*/
	};// empty class so doc++ makes a nice hierarchy.

} // namespace LinBox

#endif //  __LINBOX_blackbox_interface_H

