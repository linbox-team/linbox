/* -*- mode: c; style: linux -*- */

/* canvas-renderer.c
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

#include <ctype.h>

#include <freetype/freetype.h>
#include <gdk-pixbuf/gdk-pixbuf.h>

#include "canvas-renderer.h"

enum {
	ARG_0,
	ARG_CANVAS
};

struct _CanvasRendererPrivate 
{
	GtkWidget *canvas;
	GdkFont *font;
};

static RendererClass *parent_class;
static FT_Library library;
static FT_Face normal_face, italic_face, blex_face, blsy_face, rblmi_face;

static void canvas_renderer_init            (CanvasRenderer *canvas_renderer);
static void canvas_renderer_class_init      (CanvasRendererClass *class);

static void canvas_renderer_set_arg         (GtkObject *object, 
					     GtkArg *arg, 
					     guint arg_id);
static void canvas_renderer_get_arg         (GtkObject *object, 
					     GtkArg *arg, 
					     guint arg_id);

static void canvas_renderer_finalize        (GtkObject *object);

static void canvas_renderer_render_line     (Renderer *renderer,
					     gdouble x1, gdouble y1, 
					     gdouble x2, gdouble y2,
					     gdouble thickness);
static void canvas_renderer_render_box      (Renderer *renderer,
					     gdouble x1, gdouble y1, 
					     gdouble x2, gdouble y2,
					     gdouble thickness);
static void canvas_renderer_render_string   (Renderer *renderer,
					     const gchar *string, 
					     gdouble x, gdouble y,
					     gdouble ascent, gdouble descent);

static void canvas_renderer_get_string_geom (Renderer *renderer, gchar *string,
					     gdouble *width, gdouble *height,
					     gdouble *ascent,
					     gdouble *descent);

static FT_Face load_glyph                   (gint code, gint ascent,
					     gint descent);
static void render_glyph                    (GtkWidget *canvas, 
					     gint code, gint *x, gint y,
					     gint ascent, gint descent);
static void render_ft_bitmap                (FT_Bitmap *bitmap,
					     GdkWindow *window, GdkGC *gc, 
					     int x, int y);
static gint utf8_to_unicode                 (const gchar **string);

guint
canvas_renderer_get_type (void)
{
	static guint canvas_renderer_type = 0;

	if (!canvas_renderer_type) {
		GtkTypeInfo canvas_renderer_info = {
			"CanvasRenderer",
			sizeof (CanvasRenderer),
			sizeof (CanvasRendererClass),
			(GtkClassInitFunc) canvas_renderer_class_init,
			(GtkObjectInitFunc) canvas_renderer_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		canvas_renderer_type = 
			gtk_type_unique (renderer_get_type (), 
					 &canvas_renderer_info);
	}

	return canvas_renderer_type;
}

static void
canvas_renderer_init (CanvasRenderer *canvas_renderer)
{
	canvas_renderer->p = g_new0 (CanvasRendererPrivate, 1);
	canvas_renderer->p->font =
		gdk_font_load
		("-adobe-helvetica-medium-r-normal--24-*-*-*-*-*-iso8859-1");
}

static void
canvas_renderer_class_init (CanvasRendererClass *class) 
{
	GtkObjectClass *object_class;
	RendererClass *renderer_class;

	gtk_object_add_arg_type ("CanvasRenderer::canvas",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_CANVAS);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = canvas_renderer_finalize;
	object_class->set_arg = canvas_renderer_set_arg;
	object_class->get_arg = canvas_renderer_get_arg;

	renderer_class = RENDERER_CLASS (class);
	renderer_class->render_line = canvas_renderer_render_line;
	renderer_class->render_box = canvas_renderer_render_box;
	renderer_class->render_string = canvas_renderer_render_string;

	renderer_class->get_string_geom = canvas_renderer_get_string_geom;

	parent_class = RENDERER_CLASS
		(gtk_type_class (renderer_get_type ()));
}

static void
canvas_renderer_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	CanvasRenderer *canvas_renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CANVAS_RENDERER (object));

	canvas_renderer = CANVAS_RENDERER (object);

	switch (arg_id) {
	case ARG_CANVAS:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  GTK_IS_WIDGET (GTK_VALUE_POINTER (*arg)));

		if (canvas_renderer->p->canvas != NULL)
			gtk_object_unref
				(GTK_OBJECT (canvas_renderer->p->canvas));

		canvas_renderer->p->canvas = GTK_VALUE_POINTER (*arg);

		if (canvas_renderer->p->canvas != NULL)
			gtk_object_ref
				(GTK_OBJECT (canvas_renderer->p->canvas));
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
canvas_renderer_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	CanvasRenderer *canvas_renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CANVAS_RENDERER (object));

	canvas_renderer = CANVAS_RENDERER (object);

	switch (arg_id) {
	case ARG_CANVAS:
		GTK_VALUE_POINTER (*arg) = canvas_renderer->p->canvas;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
canvas_renderer_finalize (GtkObject *object) 
{
	CanvasRenderer *canvas_renderer;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CANVAS_RENDERER (object));

	canvas_renderer = CANVAS_RENDERER (object);

	g_free (canvas_renderer->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
canvas_renderer_new (GtkWidget *canvas) 
{
	return gtk_object_new (canvas_renderer_get_type (),
			       "canvas", canvas,
			       NULL);
}

static void
canvas_renderer_render_line (Renderer *renderer,
			     gdouble x1, gdouble y1, 
			     gdouble x2, gdouble y2,
			     gdouble thickness)
{
	CanvasRenderer *canvas_renderer;

	canvas_renderer = CANVAS_RENDERER (renderer);

	if (!GTK_WIDGET_REALIZED (canvas_renderer->p->canvas)) return;

	gdk_draw_line (canvas_renderer->p->canvas->window,
		       canvas_renderer->p->canvas->style->black_gc,
		       x1, y1, x2, y2);
}

static void
canvas_renderer_render_box (Renderer *renderer,
			    gdouble x1, gdouble y1, 
			    gdouble x2, gdouble y2,
			    gdouble thickness)
{
	CanvasRenderer *canvas_renderer;

	canvas_renderer = CANVAS_RENDERER (renderer);

	if (!GTK_WIDGET_REALIZED (canvas_renderer->p->canvas)) return;

	gdk_draw_rectangle (canvas_renderer->p->canvas->window,
			    canvas_renderer->p->canvas->style->black_gc,
			    0, x1, y1, x2, y2);
}

static void
canvas_renderer_render_string (Renderer *renderer,
			       const gchar *string, gdouble x, gdouble y,
			       gdouble ascent, gdouble descent)
{
	CanvasRenderer *canvas_renderer;
	gint current_x;

	canvas_renderer = CANVAS_RENDERER (renderer);

	if (!GTK_WIDGET_REALIZED (canvas_renderer->p->canvas)) return;

	current_x = x;

	/* N.B.: This will fail *badly* with malformed strings. It may walk
	 * off the end of the array
	 */

	while (*string)
		render_glyph (canvas_renderer->p->canvas,
			      utf8_to_unicode (&string),
			      &current_x, y, ascent, descent);
}

static void
canvas_renderer_get_string_geom (Renderer *renderer, gchar *string,
				 gdouble *width, gdouble *height,
				 gdouble *ascent, gdouble *descent)
{
	CanvasRenderer *canvas_renderer;
	gint lbearing, rbearing, i_width = 0, i_ascent = 0, i_descent = 0;
	FT_Face face;

	canvas_renderer = CANVAS_RENDERER (renderer);

	/* N.B.: See canvas_renderer_render_string above */

	while (*string) {
		face = load_glyph (utf8_to_unicode (&string), 20, 0);

		i_width += face->glyph->metrics.horiAdvance >> 6;
		i_ascent = MAX (i_ascent,
				face->glyph->metrics.horiBearingY >> 6);
		i_descent = MAX (i_descent, (face->glyph->metrics.height -
				 face->glyph->metrics.horiBearingY) >> 6);
	}

	if (width != NULL)
		*width = i_width;
	if (height != NULL)
		*height = i_ascent + i_descent;
	if (ascent != NULL)
		*ascent = i_ascent;
	if (descent != NULL)
		*descent = i_descent;
}

/**
 * load_glyph:
 * @code: 
 * @ascent: 
 * @descent: 
 * 
 * Load a glyph into the face structure
 * 
 * Return value: The font face used, or NULL on error
 **/

static FT_Face
load_glyph (gint code, gint ascent, gint descent) 
{
	gint error = 0;
	FT_UInt glyph_index;
	FT_Face face;

	if (!library)
		FT_Init_FreeType (&library);

	if (isalpha (code)) {
		if (!italic_face) {
			error = FT_New_Face (library, "timesi.ttf", 0, 
					     &italic_face);
			FT_Select_Charmap (italic_face, ft_encoding_unicode);
		}

		face = italic_face;
	}
	else if (code < 0x0080) {
		if (!normal_face) {
			error = FT_New_Face (library, "times.ttf", 0, 
					     &normal_face);
			FT_Select_Charmap (normal_face, ft_encoding_unicode);
		}

		face = normal_face;
	}
	else if (code >= 0x2200 && code <= 0x22ff) {
		if (!blex_face) {
			error = FT_New_Face (library, "blex.ttf", 0,
					     &blex_face);
			FT_Select_Charmap (blex_face, ft_encoding_unicode);
		}

		face = blex_face;
	}
	else if (code >= 0x0370 && code <= 0x03ff) {
		if (!rblmi_face) {
			error = FT_New_Face (library, "rblmi.ttf", 0,
					     &rblmi_face);
			FT_Select_Charmap (rblmi_face, ft_encoding_unicode);
		}

		face = rblmi_face;
	}

	if (error != 0) return NULL;

	FT_Set_Pixel_Sizes (face, 0, ascent + descent);
	glyph_index = FT_Get_Char_Index (face, code);
	if (FT_Load_Glyph (face, glyph_index, FT_LOAD_DEFAULT)) return NULL;

	return face;
}

/**
 * render_glyph:
 * @code: 
 * @x: 
 * @y: 
 * @ascent: 
 * @descent: 
 * 
 * Render a (Unicode) glyph on the screen and update the x-coordonite
 **/

static void
render_glyph (GtkWidget *canvas, gint code, gint *x, gint y,
	      gint ascent, gint descent) 
{
	FT_UInt error;
	FT_Face face;

	face = load_glyph (code, ascent, descent);
	if (face == NULL) return;
	error = FT_Render_Glyph (face->glyph, ft_render_mode_normal);
	if (error || !face->glyph->bitmap.buffer) return;

	render_ft_bitmap (&face->glyph->bitmap, canvas->window,
			  canvas->style->black_gc, *x, y);
	(*x) += face->glyph->advance.x >> 6;
}

/**
 * render_ft_bitmap:
 * @pixmap: 
 * @window: 
 * @gc: 
 * @x: 
 * @y: 
 * 
 * Render a Freetype-produced bitmap onto the given drawable
 **/

static void
render_ft_bitmap (FT_Bitmap *bitmap, GdkWindow *window, GdkGC *gc, 
		  int x, int y) 
{
	GdkPixbuf *pixbuf;
	GdkPixmap *gdk_pixmap, *gdk_mask;
	GdkColormap *colormap;
	GdkColor fg_color, bg_color;
	int i, j, k;
	char *buffer;

	buffer = g_new (char, bitmap->rows * bitmap->width * 3);

	for (i = 0; i < bitmap->rows * bitmap->width * 3; i++)
		buffer[i] = 0xff;

	colormap = gdk_window_get_colormap (window);
	gdk_color_black (colormap, &fg_color);
	gdk_color_white (colormap, &bg_color);

	for (i = 0; i < bitmap->rows; i++)
		for (j = 0; j < bitmap->width; j++)
			for (k = 0; k < 3; k++)
				buffer[(i * bitmap->width + j) * 3 + k] =
					(char) (255 - (int) bitmap->buffer[i * bitmap->pitch + j]);

	pixbuf = gdk_pixbuf_new_from_data
		(buffer, GDK_COLORSPACE_RGB, FALSE, 8,
		 bitmap->width, bitmap->rows, bitmap->width * 3,
		 (GdkPixbufDestroyNotify) g_free, NULL);

	gdk_pixbuf_render_pixmap_and_mask
		(pixbuf, &gdk_pixmap, &gdk_mask, 1);

	gdk_draw_pixmap (window, gc, gdk_pixmap, 0, 0, x, y, 
			 bitmap->width, bitmap->rows);

	gdk_pixbuf_unref (pixbuf);
	gdk_pixmap_unref (gdk_pixmap);
}

/**
 * utf8_to_unicode:
 * @string: 
 * 
 * Converts a UTF-8 encoded string character into unicode and updates the
 * string pointer 
 * 
 * Return value: The Unicode character
 **/

static gint
utf8_to_unicode (const gchar **string) 
{
	gint retval = 0;

	if (*(guchar *)*string < 0x80) {
		retval = **string;
		(*string)++;
	}
	else if (*(gushort *)*string < 0x800) {
		retval = ((**string & 0x1f) << 6) | (*(*string + 1) & 0x3f);
		*string += 2;
	}
	else if (*(guint *)*string < 0x10000) {
		retval = ((**string & 0x0f) << 12) |
			((*(*string + 1) & 0x3f) << 6) |
			(*(*string + 2) & 0x3f);
		*string += 3;
	} else {
		retval = ((**string & 0x07) << 18) ||
			((*(*string + 1) & 0x3f) << 12) ||
			((*(*string + 2) & 0x3f) << 6) ||
			(*(*string + 3) & 0x3f);
		*string += 4;
	}

	return retval;
}
