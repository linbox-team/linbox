/* -*- mode: c; style: linux -*- */

/* cursor.c
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

#include "cursor.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _CursorPrivate 
{
	/* Private data members */
};

static GtkObjectClass *parent_class;

static void cursor_init        (Cursor *cursor);
static void cursor_class_init  (CursorClass *class);

static void cursor_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void cursor_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void cursor_finalize    (GtkObject *object);

guint
cursor_get_type (void)
{
	static guint cursor_type = 0;

	if (!cursor_type) {
		GtkTypeInfo cursor_info = {
			"Cursor",
			sizeof (Cursor),
			sizeof (CursorClass),
			(GtkClassInitFunc) cursor_class_init,
			(GtkObjectInitFunc) cursor_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		cursor_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &cursor_info);
	}

	return cursor_type;
}

static void
cursor_init (Cursor *cursor)
{
	cursor->p = g_new0 (CursorPrivate, 1);
}

static void
cursor_class_init (CursorClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Cursor::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = cursor_finalize;
	object_class->set_arg = cursor_set_arg;
	object_class->get_arg = cursor_get_arg;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
cursor_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
cursor_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
cursor_finalize (GtkObject *object) 
{
	Cursor *cursor;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CURSOR (object));

	cursor = CURSOR (object);

	g_free (cursor->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
cursor_new (void) 
{
	return gtk_object_new (cursor_get_type (),
			       NULL);
}
