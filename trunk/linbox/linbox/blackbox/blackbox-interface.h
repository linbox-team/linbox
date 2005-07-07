/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/blackbox-interface.h
 * Copyright (C) 2002 David Saunders
 *
 * For licensing information see COPYING
 */

#ifndef __BLACKBOX_INTERFACE_H
#define __BLACKBOX_INTERFACE_H
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

#endif // __BLACKBOX_INTERFACE_H
