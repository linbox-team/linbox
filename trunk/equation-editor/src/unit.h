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

#ifndef __UNIT_H
#define __UNIT_H

#include <gnome.h>

#include "math-object.h"

BEGIN_GNOME_DECLS

#define UNIT(obj)          GTK_CHECK_CAST (obj, unit_get_type (), Unit)
#define UNIT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, unit_get_type (), UnitClass)
#define IS_UNIT(obj)       GTK_CHECK_TYPE (obj, unit_get_type ())

typedef struct _Unit Unit;
typedef struct _UnitClass UnitClass;
typedef struct _UnitPrivate UnitPrivate;

struct _Unit 
{
	MathObject parent;

	UnitPrivate *p;
};

struct _UnitClass 
{
	MathObjectClass math_object_class;
};

guint       unit_get_type        (void);

MathObject *unit_get_superscript (Unit *unit);
MathObject *unit_get_subscript   (Unit *unit);
void        unit_set_superscript (Unit *unit, MathObject *math_object);
void        unit_set_subscript   (Unit *unit, MathObject *math_object);

END_GNOME_DECLS

#endif /* __UNIT_H */
