/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/modular-balanced.h
 * Copyright (C) 2010 LinBox
 * Written by Brice Boyer <brice.boyer@imag.fr>
 *
 * See COPYING for license information.
 */

/** @file field/modular-balanced.h
 * @ingroup field
 * @brief Common header for any modular-balanced field.
 * Use <code>\#include<modular-balanced></code> to get access to any modular
 * balanced representation.
 */

#ifndef __LINBOX_modular_balanced_H
#define __LINBOX_modular_balanced_H

#include "linbox/field/FFPACK/modular-balanced-float.h"
#include "linbox/field/FFPACK/modular-balanced-double.h"
#include "linbox/field/FFPACK/modular-balanced-int32.h"
#ifdef __LINBOX_HAVE_INT64
#include "linbox/field/FFPACK/modular-balanced-int64.h"
#endif

#endif
