/* linbox/blackbox/blackbox-interface.h
 * Copyright (C) 2002 LinBox
 * Written by David Saunders
 *
 * For licensing information see COPYING
 */

#ifndef __LINBOX_blackbox_interface_H
#define __LINBOX_blackbox_interface_H

#include "linbox/element/archetype.h"

namespace LinBox
{
// LinBox Blackbox Interface
/*
 * The LinBox {@link BlackboxInterface} common object {@link Interfaces interface}.
 * The blackbox interface includes the public members defined in the archetype.
*/
/** 
 * \brief This blackbox base class exists solely to aid documentation organization.

 *  For the general blackbox member function documentation consult the {@link Blackbox Archetype}. For specific properties of individual representations consult the specific blackbox classes.
 */
class BlackboxInterface 
{
/*
    public:
	// this just demo's that some declarations could be here.
	typedef ElementArchetype Element; 
	virtual Element& mul(Element& c, const Element& a, const Element& b) const = 0;
*/
};// empty class so doc++ makes a nice hierarchy.

} // namespace LinBox

#endif //  __LINBOX_blackbox_interface_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
