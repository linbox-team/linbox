/* -*- mode: c; style: linux -*- */

/* unit-layout.h
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

#ifndef __UNIT_LAYOUT_H
#define __UNIT_LAYOUT_H

#include <gnome.h>

#include "layout.h"

BEGIN_GNOME_DECLS

#define UNIT_LAYOUT(obj)          GTK_CHECK_CAST (obj, unit_layout_get_type (), UnitLayout)
#define UNIT_LAYOUT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, unit_layout_get_type (), UnitLayoutClass)
#define IS_UNIT_LAYOUT(obj)       GTK_CHECK_TYPE (obj, unit_layout_get_type ())

typedef struct _UnitLayout UnitLayout;
typedef struct _UnitLayoutClass UnitLayoutClass;
typedef struct _UnitLayoutPrivate UnitLayoutPrivate;

struct _UnitLayout 
{
	Layout parent;

	UnitLayoutPrivate *p;
};

struct _UnitLayoutClass 
{
	LayoutClass layout_class;
};

guint unit_layout_get_type         (void);

GtkObject *unit_layout_new         (void);

END_GNOME_DECLS

#endif /* __UNIT_LAYOUT_H */
