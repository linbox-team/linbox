/* -*- mode: c; style: linux -*- */

/* unit.h
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef __MATH_UNIT_H
#define __MATH_UNIT_H

#include <gnome.h>

#include "math-object.h"

BEGIN_GNOME_DECLS

#define MATH_UNIT(obj)          GTK_CHECK_CAST (obj, math_unit_get_type (), MathUnit)
#define MATH_UNIT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, math_unit_get_type (), MathUnitClass)
#define IS_MATH_UNIT(obj)       GTK_CHECK_TYPE (obj, math_unit_get_type ())

typedef struct _MathUnit MathUnit;
typedef struct _MathUnitClass MathUnitClass;
typedef struct _MathUnitPrivate MathUnitPrivate;

struct _MathUnit 
{
	MathObject parent;

	MathUnitPrivate *p;
};

struct _MathUnitClass 
{
	MathObjectClass math_object_class;
};

guint       math_unit_get_type        (void);

MathObject *math_unit_get_superscript (MathUnit *math_unit);
MathObject *math_unit_get_subscript   (MathUnit *math_unit);
void        math_unit_set_superscript (MathUnit *math_unit,
				       MathObject *math_object);
void        math_unit_set_subscript   (MathUnit *math_unit,
				       MathObject *math_object);

END_GNOME_DECLS

#endif /* __MATH_UNIT_H */
