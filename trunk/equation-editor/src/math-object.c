/* -*- mode: c; style: linux -*- */

/* math-object.c
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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "math-object.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static GtkObjectClass *parent_class;

static void math_object_init        (MathObject *math_object);
static void math_object_class_init  (MathObjectClass *class);

static void math_object_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void math_object_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
math_object_get_type (void)
{
	static guint math_object_type = 0;

	if (!math_object_type) {
		GtkTypeInfo math_object_info = {
			"MathObject",
			sizeof (MathObject),
			sizeof (MathObjectClass),
			(GtkClassInitFunc) math_object_class_init,
			(GtkObjectInitFunc) math_object_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_object_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &math_object_info);
	}

	return math_object_type;
}

static void
math_object_init (MathObject *math_object)
{
}

static void
math_object_class_init (MathObjectClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MathObject::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = math_object_set_arg;
	object_class->get_arg = math_object_get_arg;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
math_object_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathObject *math_object;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (object));

	math_object = MATH_OBJECT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
math_object_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathObject *math_object;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (object));

	math_object = MATH_OBJECT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
math_object_new (void) 
{
	return gtk_object_new (math_object_get_type (),
			       NULL);
}
