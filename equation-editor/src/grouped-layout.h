/* -*- mode: c; style: linux -*- */

/* grouped-layout.h
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

#ifndef __GROUPED_LAYOUT_H
#define __GROUPED_LAYOUT_H

#include <gnome.h>

#include "unit-layout.h"

BEGIN_GNOME_DECLS

#define GROUPED_LAYOUT(obj)          GTK_CHECK_CAST (obj, grouped_layout_get_type (), GroupedLayout)
#define GROUPED_LAYOUT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, grouped_layout_get_type (), GroupedLayoutClass)
#define IS_GROUPED_LAYOUT(obj)       GTK_CHECK_TYPE (obj, grouped_layout_get_type ())

typedef struct _GroupedLayout GroupedLayout;
typedef struct _GroupedLayoutClass GroupedLayoutClass;
typedef struct _GroupedLayoutPrivate GroupedLayoutPrivate;

struct _GroupedLayout 
{
	UnitLayout parent;

	GroupedLayoutPrivate *p;
};

struct _GroupedLayoutClass 
{
	UnitLayoutClass unit_layout_class;
};

guint grouped_layout_get_type         (void);

GtkObject *grouped_layout_new         (void);

END_GNOME_DECLS

#endif /* __GROUPED_LAYOUT_H */
