/* -*- mode: c; style: linux -*- */

/* layout.c
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

#include "layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _LayoutPrivate 
{
	/* Private data members */
};

static GtkObjectClass *parent_class;

static void layout_init        (Layout *layout);
static void layout_class_init  (LayoutClass *class);

static void layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void layout_finalize    (GtkObject *object);

guint
layout_get_type (void)
{
	static guint layout_type = 0;

	if (!layout_type) {
		GtkTypeInfo layout_info = {
			"Layout",
			sizeof (Layout),
			sizeof (LayoutClass),
			(GtkClassInitFunc) layout_class_init,
			(GtkObjectInitFunc) layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		layout_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &layout_info);
	}

	return layout_type;
}

static void
layout_init (Layout *layout)
{
	layout->p = g_new0 (LayoutPrivate, 1);
}

static void
layout_class_init (LayoutClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Layout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = layout_finalize;
	object_class->set_arg = layout_set_arg;
	object_class->get_arg = layout_get_arg;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Layout *layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_LAYOUT (object));

	layout = LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Layout *layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_LAYOUT (object));

	layout = LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
layout_finalize (GtkObject *object) 
{
	Layout *layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_LAYOUT (object));

	layout = LAYOUT (object);

	g_free (layout->p);
}

GtkObject *
layout_new (void) 
{
	return gtk_object_new (layout_get_type (),
			       NULL);
}
