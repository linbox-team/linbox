/* -*- mode: c; style: linux -*- */

/* math-object.h
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

#ifndef __MATH_OBJECT_H
#define __MATH_OBJECT_H

#include <gnome.h>


BEGIN_GNOME_DECLS

#define MATH_OBJECT(obj)          GTK_CHECK_CAST (obj, math_object_get_type (), MathObject)
#define MATH_OBJECT_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, math_object_get_type (), MathObjectClass)
#define IS_MATH_OBJECT(obj)       GTK_CHECK_TYPE (obj, math_object_get_type ())

typedef struct _MathObject MathObject;
typedef struct _MathObjectClass MathObjectClass;

struct _MathObject 
{
	GtkObject parent;
};

struct _MathObjectClass 
{
	GtkObjectClass gtk_object_class;
};

guint math_object_get_type         (void);

GtkObject *math_object_new         (void);

END_GNOME_DECLS

#endif /* __MATH_OBJECT_H */
