/* -*- mode: c; style: linux -*- */

/* row-block-layout.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *  	      Rob Wede, Anthony Asher, Matt Spilich
 * Based on work by Chris Lahey <clahey@helixcode.com>
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

#include "row-block-layout.h"
#include "row-block.h"
#include "cursor.h"

enum {
	ARG_0,
	ARG_CURSOR
};

struct _RowBlockLayoutPrivate 
{
	Renderer      *current_renderer;
	gdouble        current_x;
	gdouble	       current_y;
	GdkRectangle  *full_area;
	GdkRectangle  *clip_area;

	gdouble       *current_width;
	gdouble       *current_height;
	gdouble       *current_ascent;
	gdouble       *current_descent;

	Cursor        *cursor;
};

static BlockLayoutClass *parent_class;

static void row_block_layout_init         (RowBlockLayout *row_block_layout);
static void row_block_layout_class_init   (RowBlockLayoutClass *class);

static void row_block_layout_set_arg      (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void row_block_layout_get_arg      (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void row_block_layout_finalize     (GtkObject *object);

static void row_block_layout_render       (Layout *layout,
					   MathObject *object,
					   Renderer *renderer,
					   GdkRectangle *full_area,
					   GdkRectangle *clip_area);

static int render_cb                      (RowBlock *block,
					   MathObject *object, 
					   RowBlockLayout *layout);

static void row_block_layout_size_request (Layout *layout,
					   Renderer *renderer,
					   MathObject *object,
					   gdouble *width,
					   gdouble *height,
					   gdouble *ascent,
					   gdouble *descent);

static int size_request_cb                (RowBlock *block,
					   MathObject *object, 
					   RowBlockLayout *layout);

guint
row_block_layout_get_type (void)
{
	static guint row_block_layout_type = 0;

	if (!row_block_layout_type) {
		GtkTypeInfo row_block_layout_info = {
			"RowBlockLayout",
			sizeof (RowBlockLayout),
			sizeof (RowBlockLayoutClass),
			(GtkClassInitFunc) row_block_layout_class_init,
			(GtkObjectInitFunc) row_block_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		row_block_layout_type = 
			gtk_type_unique (block_layout_get_type (), 
					 &row_block_layout_info);
	}

	return row_block_layout_type;
}

static void
row_block_layout_init (RowBlockLayout *row_block_layout)
{
	row_block_layout->p = g_new0 (RowBlockLayoutPrivate, 1);
}

static void
row_block_layout_class_init (RowBlockLayoutClass *class) 
{
	GtkObjectClass *object_class;
	LayoutClass *layout_class;

	gtk_object_add_arg_type ("RowBlockLayout::cursor",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_CURSOR);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = row_block_layout_finalize;
	object_class->set_arg = row_block_layout_set_arg;
	object_class->get_arg = row_block_layout_get_arg;

	layout_class = LAYOUT_CLASS (class);
	layout_class->render = row_block_layout_render;
	layout_class->size_request = row_block_layout_size_request;

	parent_class = BLOCK_LAYOUT_CLASS
		(gtk_type_class (block_layout_get_type ()));
}

static void
row_block_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	RowBlockLayout *row_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK_LAYOUT (object));

	row_block_layout = ROW_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_CURSOR:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) != NULL);
		g_return_if_fail (IS_CURSOR (GTK_VALUE_POINTER (*arg)));

		if (row_block_layout->p->cursor != NULL)
			gtk_object_unref (GTK_OBJECT
					  (row_block_layout->p->cursor));

		row_block_layout->p->cursor = GTK_VALUE_POINTER (*arg);

		if (row_block_layout->p->cursor != NULL)
			gtk_object_ref (GTK_OBJECT
					(row_block_layout->p->cursor));
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
row_block_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	RowBlockLayout *row_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK_LAYOUT (object));

	row_block_layout = ROW_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_CURSOR:
		GTK_VALUE_POINTER (*arg) = row_block_layout->p->cursor;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
row_block_layout_finalize (GtkObject *object) 
{
	RowBlockLayout *row_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK_LAYOUT (object));

	row_block_layout = ROW_BLOCK_LAYOUT (object);

	g_free (row_block_layout->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
row_block_layout_new (void) 
{
	return gtk_object_new (row_block_layout_get_type (),
			       NULL);
}

/* The following four methods adapted from work by
 * Chris Lahey <clahey@helixcode.com>
 */

static void
row_block_layout_render (Layout *layout, MathObject *object,
			 Renderer *renderer, GdkRectangle *full_area,
			 GdkRectangle *clip_area)
{
	RowBlockLayout *row_block_layout;
	gdouble width, height, ascent, descent;

	g_return_if_fail (IS_ROW_BLOCK (object));

	row_block_layout_size_request (layout, renderer, object,
				       &width, &height, &ascent, &descent);

	full_area->x = full_area->x + (full_area->width - width) / 2;
	full_area->y = full_area->y + (full_area->height - height) / 2;
	full_area->width = width;
	full_area->height = height;

	row_block_layout = ROW_BLOCK_LAYOUT (layout);
	row_block_layout->p->current_renderer = renderer;
	row_block_layout->p->current_x = full_area->x;
	row_block_layout->p->full_area = full_area;
	row_block_layout->p->clip_area = clip_area;

	if (ROW_BLOCK_LAYOUT (layout)->p->cursor != NULL &&
	    cursor_currently_in (ROW_BLOCK_LAYOUT (layout)->p->cursor, object))
		renderer_render_box (renderer, 
				     full_area->x, full_area->y,
				     full_area->width, full_area->height,
				     1);

	block_foreach (BLOCK (object), (BlockIteratorCB) render_cb,
		       row_block_layout);
}

static int
render_cb (RowBlock *block, MathObject *object, RowBlockLayout *layout) 
{
	Layout *obj_layout;
	gdouble width, height;
	GdkRectangle object_full_area;
	GdkRectangle object_clip_area;

	g_return_val_if_fail (object != NULL, 1);
	g_return_val_if_fail (IS_MATH_OBJECT (object), 1);

	obj_layout = math_object_get_layout (object);
	gtk_object_set (GTK_OBJECT (obj_layout), 
			"cursor", layout->p->cursor, NULL);
	layout_size_request (obj_layout, layout->p->current_renderer, object,
			     &width, &height, NULL, NULL);

	object_full_area.x = layout->p->current_x;
	object_full_area.y = layout->p->full_area->y;
	object_full_area.width = width;
	object_full_area.height = height;

	object_clip_area.x =
		MAX (layout->p->current_x, layout->p->clip_area->x);
	object_clip_area.y = layout->p->clip_area->y;
	object_clip_area.width = MIN (width, 
				      layout->p->clip_area->width - 
				      (layout->p->current_x - 
				       layout->p->full_area->x));
	object_clip_area.height = MIN (height, layout->p->clip_area->height);

	layout_render (obj_layout, object, layout->p->current_renderer,
		       &object_full_area, &object_clip_area);

	layout->p->current_x += width + 2;

	gtk_object_unref (GTK_OBJECT (obj_layout));

	return 0;
}

static void
row_block_layout_size_request (Layout *layout, Renderer *renderer,
			       MathObject *object,
			       gdouble *width, gdouble *height,
			       gdouble *ascent, gdouble *descent)
{
	RowBlockLayout *row_block_layout;

	g_return_if_fail (IS_ROW_BLOCK (object));

	if (ascent != NULL) *ascent = 0;
	if (descent != NULL) *descent = 0;
	if (width != NULL) *width = 0;
	if (height != NULL) *height = 0;

	row_block_layout = ROW_BLOCK_LAYOUT (layout);
	row_block_layout->p->current_width = width;
	row_block_layout->p->current_height = height;
	row_block_layout->p->current_ascent = ascent;
	row_block_layout->p->current_descent = descent;
	row_block_layout->p->current_renderer = renderer;

	block_foreach (BLOCK (object), (BlockIteratorCB) size_request_cb,
		       row_block_layout);
}

static int
size_request_cb (RowBlock *block, MathObject *object, RowBlockLayout *layout) 
{
	Layout *obj_layout;
	double obj_width, obj_height, obj_ascent, obj_descent;

	obj_layout = math_object_get_layout (object);

	layout_size_request (obj_layout, layout->p->current_renderer, object,
			     &obj_width, &obj_height,
			     &obj_ascent, &obj_descent);

	if (layout->p->current_width != NULL)
		*layout->p->current_width += obj_width + 2;
	if (layout->p->current_height != NULL && 
	    obj_ascent > *layout->p->current_height)
		*layout->p->current_height = obj_height;
	if (layout->p->current_ascent != NULL && 
	    obj_ascent > *layout->p->current_ascent)
		*layout->p->current_ascent = obj_ascent;
	if (layout->p->current_descent != NULL && 
	    obj_descent > *layout->p->current_descent)
		*layout->p->current_descent = obj_descent;

	return 0;
}
