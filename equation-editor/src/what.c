/* -*- mode: c; style: linux -*- */

/* what.c
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

#include "what.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _WhatPrivate 
{
	/* Private data members */
};

static Class *parent_class;

static void what_init        (What *what);
static void what_class_init  (WhatClass *class);

static void what_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void what_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void what_finalize    (GtkObject *object);

guint
what_get_type (void)
{
	static guint what_type = 0;

	if (!what_type) {
		GtkTypeInfo what_info = {
			"What",
			sizeof (What),
			sizeof (WhatClass),
			(GtkClassInitFunc) what_class_init,
			(GtkObjectInitFunc) what_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		what_type = 
			gtk_type_unique (_get_type (), 
					 &what_info);
	}

	return what_type;
}

static void
what_init (What *what)
{
	what->p = g_new0 (WhatPrivate, 1);
}

static void
what_class_init (WhatClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("What::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = what_finalize;
	object_class->set_arg = what_set_arg;
	object_class->get_arg = what_get_arg;

	parent_class = _CLASS
		(gtk_type_class (_get_type ()));
}

static void
what_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	What *what;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WHAT (object));

	what = WHAT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
what_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	What *what;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WHAT (object));

	what = WHAT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
what_finalize (GtkObject *object) 
{
	What *what;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WHAT (object));

	what = WHAT (object);

	g_free (what->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
what_new (void) 
{
	return gtk_object_new (what_get_type (),
			       NULL);
}
