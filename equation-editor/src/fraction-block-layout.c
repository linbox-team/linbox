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

static void fraction_block_layout_render(Layout *layout,
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
fraction_block_layout_render(Layout *layout,
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

	num_full_area.x = full_area->x;
	num_full_area.y = full_area->y;
	num_full_area.width = num_w;
	num_full_area.height = num_h;

	num_clip_area.x = full_area->x;
	num_clip_area.y = full_area->y;
	num_clip_area.width = num_w;
	num_clip_area.height = num_h;

	den_full_area.x = full_area->x;
	den_full_area.y = full_area->y + 40;
	den_full_area.width = den_w;
	den_full_area.height = den_h;

	den_clip_area.x = full_area->x;
	den_clip_area.y = full_area->y + 20;
	den_clip_area.width = den_w;
	den_clip_area.height = den_h;


	object_N = fraction_block_get_numerator( FRACTION_BLOCK (object));

	object_D = fraction_block_get_denominator( FRACTION_BLOCK (object));

	obj_layout_N = math_object_get_layout (object_N);
	obj_layout_D = math_object_get_layout (object_D);
	
	layout_size_request ( LAYOUT(obj_layout_N), renderer, object_N,
			&num_w, &num_h, NULL, NULL);

	layout_size_request ( LAYOUT(obj_layout_D), renderer, object_D,
			&den_w, &den_h, NULL, NULL);


	renderer_render_line( renderer, full_area->x+55, 
	    full_area->y + num_h + 40, full_area->x + MAX(num_w, den_w)+40,
	    full_area->y + num_h + 40, 20);

	layout_render (LAYOUT(obj_layout_N),object_N,renderer,
		       &num_full_area, &num_full_area);

	layout_render (LAYOUT(obj_layout_D),object_D,renderer,
		       &den_full_area, &den_full_area);
	
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
        fraction_block_layout->p->current_width = width;
        fraction_block_layout->p->current_height = height + 10;
        /* the plus 10 constant is the approximate height of the line */
        fraction_block_layout->p->current_ascent = ascent;
        fraction_block_layout->p->current_descent = descent;
        fraction_block_layout->p->current_renderer = renderer;

        block_foreach (BLOCK (object), (BlockIteratorCB) size_request_cb,
                       fraction_block_layout);
}

static int 
size_request_cb ( FractionBlock *block, MathObject *object,
                  FractionBlockLayout *layout)
{

        Layout *obj_layout;
        double obj_width, obj_height, obj_ascent, obj_descent;

        obj_layout = math_object_get_layout (object);

        layout_size_request (obj_layout, layout->p->current_renderer, object,
                             &obj_width, &obj_height,
                             &obj_ascent, &obj_descent);

        if (layout->p->current_width != NULL)
                *layout->p->current_width = MAX (obj_width,
                                *layout->p->current_width) + 5;

        if (layout->p->current_height != NULL &&
            obj_ascent > *layout->p->current_height)
                *layout->p->current_height += obj_height;
        if (layout->p->current_ascent != NULL &&
            obj_ascent > *layout->p->current_ascent)
                *layout->p->current_ascent = obj_ascent;
        if (layout->p->current_descent != NULL &&
            obj_descent > *layout->p->current_descent)
                *layout->p->current_descent = obj_descent;

        return 0;
 
}

