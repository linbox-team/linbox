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

static void layout_init              (Layout *layout);
static void layout_class_init        (LayoutClass *class);

static void layout_set_arg           (GtkObject *object, 
				      GtkArg *arg, 
				      guint arg_id);
static void layout_get_arg           (GtkObject *object, 
				      GtkArg *arg, 
				      guint arg_id);

static void layout_finalize          (GtkObject *object);

static void layout_real_render       (Layout *layout, 
				      MathObject *math_object,
				      Renderer *renderer,
				      GdkRectangle *full_area,
				      GdkRectangle *clip_area);

static void layout_real_size_request (Layout *layout,
				      Renderer *renderer,
				      MathObject *math_object,
				      gdouble *width,
				      gdouble *height,
				      gdouble *ascent,
				      gdouble *descent);

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

	class->render = layout_real_render;
	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = layout_finalize;
	object_class->set_arg = layout_set_arg;
	object_class->get_arg = layout_get_arg;

	class->render = layout_real_render;
	class->size_request = layout_real_size_request;

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

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
layout_new (void) 
{
	return gtk_object_new (layout_get_type (),
			       NULL);
}

void
layout_render (Layout *layout, MathObject *math_object, Renderer *renderer,
	       GdkRectangle *full_area, GdkRectangle *clip_area) 
{
	g_return_if_fail (layout != NULL);
	g_return_if_fail (IS_LAYOUT (layout));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));
	g_return_if_fail (full_area != NULL);

	LAYOUT_CLASS (GTK_OBJECT (layout)->klass)->render
		(layout, math_object, renderer, full_area, clip_area);
}

void
layout_size_request (Layout *layout, Renderer *renderer,
		     MathObject *math_object,
		     gdouble *width, gdouble *height,
		     gdouble *ascent, gdouble *descent)
{
	g_return_if_fail (layout != NULL);
	g_return_if_fail (IS_LAYOUT (layout));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	LAYOUT_CLASS (GTK_OBJECT (layout)->klass)->size_request
		(layout, renderer, math_object, 
		 width, height, ascent, descent);
}

static void
layout_real_render (Layout *layout, MathObject *math_object, 
		    Renderer *renderer,
		    GdkRectangle *full_area, GdkRectangle *clip_area)
{
	g_warning("Pure virtual method Layout::render invoked");
}

static void
layout_real_size_request (Layout *layout, Renderer *renderer, 
			  MathObject *math_object,
			  gdouble *width, gdouble *height,
			  gdouble *ascent, gdouble *descent)
{
	g_warning("Pure virtual method Layout::size_request invoked");
}

