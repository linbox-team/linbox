/* -*- mode: c; style: linux -*- */

/* fraction-block-layout.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *            Anthony Asher, Rob Wehde, Matt Spilich
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

#include "fraction-block-layout.h"
#include "fraction-block.h"
#include "cursor.h"

enum {
	ARG_0,
	ARG_CURSOR
};

struct _FractionBlockLayoutPrivate 
{
	Cursor        *cursor;
};

static BlockLayoutClass *parent_class;

static void fraction_block_layout_init        (FractionBlockLayout *fraction_block_layout);
static void fraction_block_layout_class_init  (FractionBlockLayoutClass *class);

static void fraction_block_layout_set_arg (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void fraction_block_layout_get_arg (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void fraction_block_layout_finalize (GtkObject *object);

static void fraction_block_layout_render   (Layout *layout,
					    MathObject *object,
					    Renderer *renderer,
					    GdkRectangle *full_area,
					    GdkRectangle *clip_area);

static void fraction_block_layout_size_request (Layout *layout,
						Renderer *renderer,
						MathObject *object,
						gdouble *width,
						gdouble *height,
						gdouble *ascent,
						gdouble *descent);

guint
fraction_block_layout_get_type (void)
{
	static guint fraction_block_layout_type = 0;

	if (!fraction_block_layout_type) {
		GtkTypeInfo fraction_block_layout_info = {
			"FractionBlockLayout",
			sizeof (FractionBlockLayout),
			sizeof (FractionBlockLayoutClass),
			(GtkClassInitFunc) fraction_block_layout_class_init,
			(GtkObjectInitFunc) fraction_block_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		fraction_block_layout_type = 
			gtk_type_unique (block_layout_get_type (), 
					 &fraction_block_layout_info);
	}

	return fraction_block_layout_type;
}

static void
fraction_block_layout_init (FractionBlockLayout *fraction_block_layout)
{
	fraction_block_layout->p = g_new0 (FractionBlockLayoutPrivate, 1);
}

static void
fraction_block_layout_class_init (FractionBlockLayoutClass *class) 
{
	GtkObjectClass *object_class;
	LayoutClass *layout_class;

	gtk_object_add_arg_type ("FractionBlockLayout::cursor",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_CURSOR);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = fraction_block_layout_finalize;
	object_class->set_arg = fraction_block_layout_set_arg;
	object_class->get_arg = fraction_block_layout_get_arg;

	layout_class = LAYOUT_CLASS (class);
	layout_class->render = fraction_block_layout_render;
	layout_class->size_request = fraction_block_layout_size_request;

	parent_class = BLOCK_LAYOUT_CLASS
		(gtk_type_class (block_layout_get_type ()));
}

static void
fraction_block_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	FractionBlockLayout *fraction_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK_LAYOUT (object));

	fraction_block_layout = FRACTION_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_CURSOR:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) != NULL);
		g_return_if_fail (IS_CURSOR (GTK_VALUE_POINTER (*arg)));

		if (fraction_block_layout->p->cursor != NULL)
			gtk_object_unref (GTK_OBJECT
					  (fraction_block_layout->p->cursor));

		fraction_block_layout->p->cursor = GTK_VALUE_POINTER (*arg);

		if (fraction_block_layout->p->cursor != NULL)
			gtk_object_ref (GTK_OBJECT
					(fraction_block_layout->p->cursor));
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
fraction_block_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	FractionBlockLayout *fraction_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK_LAYOUT (object));

	fraction_block_layout = FRACTION_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_CURSOR:
		GTK_VALUE_POINTER (*arg) = fraction_block_layout->p->cursor;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
fraction_block_layout_finalize (GtkObject *object) 
{
	FractionBlockLayout *fraction_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK_LAYOUT (object));

	fraction_block_layout = FRACTION_BLOCK_LAYOUT (object);

	g_free (fraction_block_layout->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
fraction_block_layout_new (void) 
{
	return gtk_object_new (fraction_block_layout_get_type (),
			       NULL);
}

/* The following adapted from work by
 * Chris Lahey <clahey@helixcode.com>
 */

static void
fraction_block_layout_render (Layout *layout,
			      MathObject *object,
			      Renderer *renderer,
			      GdkRectangle *full_area,
			      GdkRectangle *clip_area)
{
	FractionBlockLayout *fraction_block_layout;
	Layout *obj_layout_N;
	Layout *obj_layout_D;
	MathObject *object_N;
	MathObject *object_D;
	GdkRectangle num_full_area;
	GdkRectangle den_full_area;
	GdkRectangle num_clip_area;
	GdkRectangle den_clip_area;

	gdouble num_w;
	gdouble den_w;
	gdouble num_h;
	gdouble den_h;

	g_return_if_fail (IS_FRACTION_BLOCK (object));

	fraction_block_layout = FRACTION_BLOCK_LAYOUT (layout);

	object_N = fraction_block_get_numerator (FRACTION_BLOCK (object));
	object_D = fraction_block_get_denominator (FRACTION_BLOCK (object));

	obj_layout_N = math_object_get_layout (object_N);
	obj_layout_D = math_object_get_layout (object_D);

	gtk_object_set (GTK_OBJECT (obj_layout_N), "cursor",
			fraction_block_layout->p->cursor, NULL);
	gtk_object_set (GTK_OBJECT (obj_layout_D), "cursor",
			fraction_block_layout->p->cursor, NULL);
	
	layout_size_request (LAYOUT (obj_layout_N), renderer, object_N,
			     &num_w, &num_h, NULL, NULL);

	layout_size_request (LAYOUT (obj_layout_D), renderer, object_D,
			     &den_w, &den_h, NULL, NULL);

	num_full_area.x = full_area->x;
	num_full_area.y = full_area->y;
	num_full_area.width = full_area->width;
	num_full_area.height = num_h;

	den_full_area.x = full_area->x;
	den_full_area.y = full_area->y + num_h + 10;
	den_full_area.width = full_area->width;
	den_full_area.height = den_h;

	renderer_render_line (renderer, full_area->x, 
			      full_area->y + num_h + 5, 
			      full_area->x + full_area->width,
			      full_area->y + num_h + 5, 1);

	layout_render (LAYOUT (obj_layout_N), object_N, renderer,
		       &num_full_area, &num_clip_area);

	layout_render (LAYOUT (obj_layout_D), object_D, renderer,
		       &den_full_area, &den_clip_area);
	
}

static void
fraction_block_layout_size_request (Layout *layout, Renderer *renderer,
				    MathObject *object,
				    gdouble *width, gdouble *height,
				    gdouble *ascent, gdouble *descent)
{
        FractionBlockLayout *fraction_block_layout;
	MathObject *num, *den;
	Layout *obj_layout;
	gdouble n_width, n_height, d_width, d_height;

        g_return_if_fail (IS_FRACTION_BLOCK (object));

        fraction_block_layout = FRACTION_BLOCK_LAYOUT (layout);

	num = fraction_block_get_numerator (FRACTION_BLOCK (object));
	den = fraction_block_get_denominator (FRACTION_BLOCK (object));

	obj_layout = math_object_get_layout (num);
	layout_size_request (obj_layout, renderer, num,
			     &n_width, &n_height, NULL, NULL);
	gtk_object_unref (GTK_OBJECT (obj_layout));

	obj_layout = math_object_get_layout (den);
	layout_size_request (obj_layout, renderer, den,
			     &d_width, &d_height, NULL, NULL);
	gtk_object_unref (GTK_OBJECT (obj_layout));

	if (width != NULL)
		*width = MAX (n_width, d_width) + 10;
	if (height != NULL)
		*height = n_height + d_height + 10;
	if (ascent != NULL)
		*ascent = n_height + d_height + 10;
	if (descent != NULL)
		*descent = 0;
}
