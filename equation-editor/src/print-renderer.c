/* -*- mode: c; style: linux -*- */

/* print-renderer.c
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

#include "print-renderer.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static RendererClass *renderer_class;

static void print_renderer_init        (PrintRenderer *print_renderer);
static void print_renderer_class_init  (PrintRendererClass *class);

static void print_renderer_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void print_renderer_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
print_renderer_get_type (void)
{
	static guint print_renderer_type = 0;

	if (!print_renderer_type) {
		GtkTypeInfo print_renderer_info = {
			"PrintRenderer",
			sizeof (PrintRenderer),
			sizeof (PrintRendererClass),
			(GtkClassInitFunc) print_renderer_class_init,
			(GtkObjectInitFunc) print_renderer_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		print_renderer_type = 
			gtk_type_unique (renderer_get_type (), 
					 &print_renderer_info);
	}

	return print_renderer_type;
}

static void
print_renderer_init (PrintRenderer *print_renderer)
{
}

static void
print_renderer_class_init (PrintRendererClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("PrintRenderer::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = print_renderer_set_arg;
	object_class->get_arg = print_renderer_get_arg;

	parent_class = RENDERER_CLASS
		(gtk_type_class (renderer_get_type ()));
}

static void
print_renderer_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	PrintRenderer *print_renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_PRINT_RENDERER (object));

	print_renderer = PRINT_RENDERER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
print_renderer_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	PrintRenderer *print_renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_PRINT_RENDERER (object));

	print_renderer = PRINT_RENDERER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
print_renderer_new (void) 
{
	return gtk_object_new (print_renderer_get_type (),
			       NULL);
}
