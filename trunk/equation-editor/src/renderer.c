/* -*- mode: c; style: linux -*- */

/* renderer.c
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

#include "renderer.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _RendererPrivate 
{
	/* Private data members */
};

static GtkObjectClass *parent_class;

static void renderer_init                 (Renderer *renderer);
static void renderer_class_init           (RendererClass *class);

static void renderer_set_arg              (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void renderer_get_arg              (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void renderer_finalize             (GtkObject *object);

static void renderer_real_render_line     (Renderer *renderer,
					   gdouble x1, gdouble y1, 
					   gdouble x2, gdouble y2,
					   gdouble thickness);
static void renderer_real_render_glyph    (Renderer *renderer,
					   gint code, gdouble x, gdouble y,
					   gdouble scale);
static void renderer_real_render_number   (Renderer *renderer,
					   gdouble value, gdouble x, gdouble y,
					   gdouble scale, gdouble pres);
static void renderer_real_render_string   (Renderer *renderer,
					   const gchar *string, 
					   gdouble x, gdouble y,
					   gdouble scale);

static void renderer_real_get_glyph_geom  (Renderer *renderer, gint code,
					   gdouble *width, gdouble *height,
					   gdouble *ascent, gdouble *descent);
static void renderer_real_get_number_geom (Renderer *renderer, gdouble value,
					   gdouble *width, gdouble *height,
					   gdouble *ascent, gdouble *descent);
static void renderer_real_get_string_geom (Renderer *renderer, gchar *string,
					   gdouble *width, gdouble *height,
					   gdouble *ascent, gdouble *descent);

guint
renderer_get_type (void)
{
	static guint renderer_type = 0;

	if (!renderer_type) {
		GtkTypeInfo renderer_info = {
			"Renderer",
			sizeof (Renderer),
			sizeof (RendererClass),
			(GtkClassInitFunc) renderer_class_init,
			(GtkObjectInitFunc) renderer_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		renderer_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &renderer_info);
	}

	return renderer_type;
}

static void
renderer_init (Renderer *renderer)
{
	renderer->p = g_new0 (RendererPrivate, 1);
}

static void
renderer_class_init (RendererClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Renderer::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = renderer_finalize;
	object_class->set_arg = renderer_set_arg;
	object_class->get_arg = renderer_get_arg;

	class->render_line = renderer_real_render_line;
	class->render_glyph = renderer_real_render_glyph;
	class->render_number = renderer_real_render_number;
	class->render_string = renderer_real_render_string;

	class->get_glyph_geom = renderer_real_get_glyph_geom;
	class->get_number_geom = renderer_real_get_number_geom;
	class->get_string_geom = renderer_real_get_string_geom;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
renderer_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Renderer *renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_RENDERER (object));

	renderer = RENDERER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
renderer_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Renderer *renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_RENDERER (object));

	renderer = RENDERER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
renderer_finalize (GtkObject *object) 
{
	Renderer *renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_RENDERER (object));

	renderer = RENDERER (object);

	g_free (renderer->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

/**
 * renderer_render_line:
 * @renderer: 
 * @x1: 
 * @y1: 
 * @x2: 
 * @y2: 
 * @thickness: 
 * 
 * Render a line with the specified endpoints and thickness onto the specified 
 * canvas
 **/

void
renderer_render_line (Renderer *renderer,
		      gdouble x1, gdouble y1, 
		      gdouble x2, gdouble y2,
		      gdouble thickness)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->render_line
		(renderer, x1, y1, x2, y2, thickness);
}

/**
 * renderer_render_glyph:
 * @renderer: 
 * @code: 
 * @x: 
 * @y: 
 * @scale: 
 * 
 * Render a glyph on the canvas, presented in Unicode
 **/

void
renderer_render_glyph (Renderer *renderer,
		       gint code, gdouble x, gdouble y,
		       gdouble scale)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->render_glyph
		(renderer, code, x, y, scale);
}

/**
 * renderer_render_number:
 * @renderer: 
 * @value: 
 * @x: 
 * @y: 
 * @scale: 
 * @pres: 
 * 
 * Render a number to the canvas
 **/

void
renderer_render_number (Renderer *renderer,
			gdouble value, gdouble x, gdouble y,
			gdouble scale, gdouble pres)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->render_number
		(renderer, value, x, y, scale, pres);
}

/**
 * renderer_render_string:
 * @renderer: 
 * @string: 
 * @x: 
 * @y: 
 * @scale: 
 * 
 * Render a string to the canvas
 **/

void
renderer_render_string (Renderer *renderer,
			const gchar *string, gdouble x, gdouble y,
			gdouble scale)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->render_string
		(renderer, string, x, y, scale);
}

void
renderer_get_glyph_geom (Renderer *renderer, gint code,
			 gdouble *width, gdouble *height,
			 gdouble *ascent, gdouble *descent)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->get_glyph_geom
		(renderer, code, width, height, ascent, descent);
}

void
renderer_get_number_geom (Renderer *renderer, gdouble value,
			  gdouble *width, gdouble *height,
			  gdouble *ascent, gdouble *descent)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->get_number_geom
		(renderer, value, width, height, ascent, descent);
}

void
renderer_get_string_geom (Renderer *renderer, gchar *string,
			  gdouble *width, gdouble *height,
			  gdouble *ascent, gdouble *descent)
{
	g_return_if_fail (renderer != NULL);
	g_return_if_fail (IS_RENDERER (renderer));

	RENDERER_CLASS (GTK_OBJECT (renderer)->klass)->get_string_geom
		(renderer, string, width, height, ascent, descent);
}

static void
renderer_real_render_line (Renderer *renderer,
			   gdouble x1, gdouble y1, 
			   gdouble x2, gdouble y2,
			   gdouble thickness)
{
	g_warning ("Pure virtual method Renderer::render_line called");
}

static void
renderer_real_render_glyph (Renderer *renderer,
			    gint code, gdouble x, gdouble y,
			    gdouble scale)
{
	g_warning ("Pure virtual method Renderer::render_glyph called");
}

static void
renderer_real_render_number (Renderer *renderer,
			     gdouble value, gdouble x, gdouble y,
			     gdouble scale, gdouble pres)
{
	g_warning ("Pure virtual method Renderer::render_number called");
}

static void
renderer_real_render_string (Renderer *renderer,
			     const gchar *string, gdouble x, gdouble y,
			     gdouble scale)
{
	g_warning ("Pure virtual method Renderer::render_string called");
}

static void
renderer_real_get_glyph_geom (Renderer *renderer, gint code,
			      gdouble *width, gdouble *height,
			      gdouble *ascent, gdouble *descent)
{
	g_warning ("Pure virtual method Renderer::get_glyph_geom called");
}

static void
renderer_real_get_number_geom (Renderer *renderer, gdouble value,
			       gdouble *width, gdouble *height,
			       gdouble *ascent, gdouble *descent)
{
	g_warning ("Pure virtual method Renderer::get_number_geom called");
}

static void
renderer_real_get_string_geom (Renderer *renderer, gchar *string,
			       gdouble *width, gdouble *height,
			       gdouble *ascent, gdouble *descent)
{
	g_warning ("Pure virtual method Renderer::get_string_geom called");
}

