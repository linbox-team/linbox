/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/util/timer.h
 * Copyright (C) 1994-1997 Givaro Team
 *
 * Written by T. Gautier
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added _start_t member to BaseTimer, so that stop () does not clobber the
 * class' memory of its start time. This allows it to be called repeatedly to
 * get elapsed times.
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 */

/*! @file util/timer.h
 * @ingroup util
 * LinBox timer is Givaro's.
 */

#ifndef __LINBOX_timer_H
#define __LINBOX_timer_H

#include <givaro/givtimer.h>

namespace LinBox {
typedef ::Givaro::Timer Timer  ;
typedef ::Givaro::BaseTimer BaseTimer ;
typedef ::Givaro::UserTimer UserTimer ;
typedef ::Givaro::SysTimer SysTimer ;
}

#endif  //__LINBOX_timer_H
