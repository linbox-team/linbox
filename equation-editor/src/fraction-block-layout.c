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

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _FractionBlockLayoutPrivate 
{
	/* Private data members */
	Renderer      *current_renderer;
	gdouble        current_x;
	GdkRectangle  *full_area;
	GdkRectangle  *clip_area;

	gdouble       *current_width;
	gdouble       *current_height;
	gdouble       *current_ascent;
	gdouble       *current_descent;
};

static BlockLayoutClass *parent_class;

static void fraction_block_layout_init        (FractionBlockLayout *fraction_block_layout);
static void fraction_block_layout_class_init  (FractionBlockLayoutClass *class);

static void fraction_block_layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void fraction_block_layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void fraction_block_layout_finalize    (GtkObject *object);

static void fraction_block_layout_render       (Layout *layout,
					   MathObject *object,
					   Renderer *renderer,
					   GdkRectangle *full_area,
					   GdkRectangle *clip_area);

static int render_cb                      (FractionBlock *block,
					   MathObject *object, 
					   FractionBlockLayout *layout);

static void fraction_block_layout_size_request (Layout *layout,
					   Renderer *renderer,
					   MathObject *object,
					   gdouble *width,
					   gdouble *height,
					   gdouble *ascent,
					   gdouble *descent);

static int size_request_cb                (FractionBlock *block,
					   MathObject *object, 
					   FractionBlockLayout *layout);


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

	gtk_object_add_arg_type ("FractionBlockLayout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

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
	case ARG_SAMPLE:
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
	case ARG_SAMPLE:
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
fraction_block_layout_render (Layout *layout, MathObject *object,
			 Renderer *renderer, GdkRectangle *full_area,
			 GdkRectangle *clip_area)
{
	FractionBlockLayout *fraction_block_layout;
	Layout *obj_layout_N;
	Layout *obj_layout_D;
	MathObject *object_D;
	MathObject *object_N;
	GdkRectangle num_area;
	GdkRectangle den_area;

	gdouble num_w;
	gdouble den_w;
	gdouble num_h;
	gdouble den_h;
	gdouble ascent;
	gdouble descent;

	g_return_if_fail (IS_FRACTION_BLOCK (object));


	fraction_block_layout = FRACTION_BLOCK_LAYOUT (layout);
	fraction_block_layout->p->current_renderer = renderer;
	fraction_block_layout->p->current_x = full_area->x;
	fraction_block_layout->p->full_area = full_area;
	fraction_block_layout->p->clip_area = clip_area;
	
	printf("call to get numerator\n");
	object_N = fraction_block_get_numerator( FRACTION_BLOCK (object));
	printf("done\n");

	printf("call to get denominator\n");
	object_D = fraction_block_get_denominator( FRACTION_BLOCK (object));
	printf("done\n");

	obj_layout_N = math_object_get_layout (object_N);
	obj_layout_D = math_object_get_layout (object_D);
	
	layout_size_request ( LAYOUT(obj_layout_N), renderer, object_N,
			&num_w, &num_h, &ascent, &descent);

	layout_size_request ( LAYOUT(obj_layout_D), renderer, object_D,
			&den_w, &den_h, &ascent, &descent);

	g_print ("rendering line\n");

	renderer_render_line( renderer, clip_area->x, 
	    clip_area->y + num_h + 25, clip_area->x + MAX(num_w, den_w),
	    clip_area->y + num_h + 25, 20);

	printf("line rendered\n");
	printf("rendering numerator\n");

	layout_render (LAYOUT(obj_layout_N),object_N,renderer,
		       &num_area, &clip_area);

	printf("numerator rendered\n");
	printf("rendering denominator\n");

	layout_render (LAYOUT(obj_layout_D),object_D,renderer,
		       &den_area, &clip_area);
	
	printf("denominator rendered\n");
}

static void
fraction_block_layout_size_request (Layout *layout, Renderer *renderer,
			       MathObject *object,
			       gdouble *width, gdouble *height,
			       gdouble *ascent, gdouble *descent)
{
	FractionBlockLayout *fraction_block_layout;

	g_return_if_fail (IS_FRACTION_BLOCK (object));

	if (ascent != NULL) *ascent = 0;
	if (descent != NULL) *descent = 0;
	if (width != NULL) *width = 0;
	if (height != NULL) *height = 0;

	fraction_block_layout = FRACTION_BLOCK_LAYOUT (layout);
}
