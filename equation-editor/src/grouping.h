/* -*- mode: c; style: linux -*- */

/* grouping.h
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

#ifndef __GROUPING_H
#define __GROUPING_H

#include <gnome.h>

#include "unit.h"

BEGIN_GNOME_DECLS

#define GROUPING(obj)          GTK_CHECK_CAST (obj, grouping_get_type (), Grouping)
#define GROUPING_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, grouping_get_type (), GroupingClass)
#define IS_GROUPING(obj)       GTK_CHECK_TYPE (obj, grouping_get_type ())

typedef struct _Grouping Grouping;
typedef struct _GroupingClass GroupingClass;
typedef struct _GroupingPrivate GroupingPrivate;

struct _Grouping 
{
	Unit parent;

	GroupingPrivate *p;
};

struct _GroupingClass 
{
	UnitClass unit_class;
};

guint grouping_get_type         (void);

GtkObject *grouping_new         (void);

END_GNOME_DECLS

#endif /* __GROUPING_H */
