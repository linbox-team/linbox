/* -*- mode: c; style: linux -*- */

/* widget-class-name.h
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

#ifndef __WIDGET_CLASS_NAME_H
#define __WIDGET_CLASS_NAME_H

#include <gnome.h>

#include "parent-class-name.h"

BEGIN_GNOME_DECLS

#define WIDGET_CLASS_NAME(obj)          GTK_CHECK_CAST (obj, widget_class_name_get_type (), WidgetClassName)
#define WIDGET_CLASS_NAME_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, widget_class_name_get_type (), WidgetClassNameClass)
#define IS_WIDGET_CLASS_NAME(obj)       GTK_CHECK_TYPE (obj, widget_class_name_get_type ())

typedef struct _WidgetClassName WidgetClassName;
typedef struct _WidgetClassNameClass WidgetClassNameClass;
typedef struct _WidgetClassNamePrivate WidgetClassNamePrivate;

struct _WidgetClassName 
{
	ParentClassName parent;

	WidgetClassNamePrivate *p;
};

struct _WidgetClassNameClass 
{
	ParentClassNameClass parent_class_name_class;
};

guint widget_class_name_get_type         (void);

GtkObject *widget_class_name_new         (void);

END_GNOME_DECLS

#endif /* __WIDGET_CLASS_NAME_H */
