/* -*- mode: c; style: linux -*- */

/* identifier.c
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

#include "identifier.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static MathObjectClass *math_object_class;

static void identifier_init        (Identifier *identifier);
static void identifier_class_init  (IdentifierClass *class);

static void identifier_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void identifier_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
identifier_get_type (void)
{
	static guint identifier_type = 0;

	if (!identifier_type) {
		GtkTypeInfo identifier_info = {
			"Identifier",
			sizeof (Identifier),
			sizeof (IdentifierClass),
			(GtkClassInitFunc) identifier_class_init,
			(GtkObjectInitFunc) identifier_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		identifier_type = 
			gtk_type_unique (math_object_get_type (), 
					 &identifier_info);
	}

	return identifier_type;
}

static void
identifier_init (Identifier *identifier)
{
}

static void
identifier_class_init (IdentifierClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Identifier::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = identifier_set_arg;
	object_class->get_arg = identifier_get_arg;

	parent_class = MATH_OBJECT_CLASS
		(gtk_type_class (math_object_get_type ()));
}

static void
identifier_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Identifier *identifier;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_IDENTIFIER (object));

	identifier = IDENTIFIER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
identifier_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Identifier *identifier;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_IDENTIFIER (object));

	identifier = IDENTIFIER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
identifier_new (void) 
{
	return gtk_object_new (identifier_get_type (),
			       NULL);
}
