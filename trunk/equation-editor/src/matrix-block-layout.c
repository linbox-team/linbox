/* -*- mode: c; style: linux -*- */

/* matrix-block-layout.c
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

#include "matrix-block-layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _MatrixBlockLayoutPrivate 
{
	/* Private data members */
};

static BlockLayoutClass *parent_class;

static void matrix_block_layout_init        (MatrixBlockLayout *matrix_block_layout);
static void matrix_block_layout_class_init  (MatrixBlockLayoutClass *class);

static void matrix_block_layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void matrix_block_layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void matrix_block_layout_finalize    (GtkObject *object);

guint
matrix_block_layout_get_type (void)
{
	static guint matrix_block_layout_type = 0;

	if (!matrix_block_layout_type) {
		GtkTypeInfo matrix_block_layout_info = {
			"MatrixBlockLayout",
			sizeof (MatrixBlockLayout),
			sizeof (MatrixBlockLayoutClass),
			(GtkClassInitFunc) matrix_block_layout_class_init,
			(GtkObjectInitFunc) matrix_block_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		matrix_block_layout_type = 
			gtk_type_unique (block_layout_get_type (), 
					 &matrix_block_layout_info);
	}

	return matrix_block_layout_type;
}

static void
matrix_block_layout_init (MatrixBlockLayout *matrix_block_layout)
{
	matrix_block_layout->p = g_new0 (MatrixBlockLayoutPrivate, 1);
}

static void
matrix_block_layout_class_init (MatrixBlockLayoutClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MatrixBlockLayout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = matrix_block_layout_finalize;
	object_class->set_arg = matrix_block_layout_set_arg;
	object_class->get_arg = matrix_block_layout_get_arg;

	parent_class = BLOCK_LAYOUT_CLASS
		(gtk_type_class (block_layout_get_type ()));
}

static void
matrix_block_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlockLayout *matrix_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK_LAYOUT (object));

	matrix_block_layout = MATRIX_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
matrix_block_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlockLayout *matrix_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK_LAYOUT (object));

	matrix_block_layout = MATRIX_BLOCK_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
matrix_block_layout_finalize (GtkObject *object) 
{
	MatrixBlockLayout *matrix_block_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK_LAYOUT (object));

	matrix_block_layout = MATRIX_BLOCK_LAYOUT (object);

	g_free (matrix_block_layout->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
matrix_block_layout_new (void) 
{
	return gtk_object_new (matrix_block_layout_get_type (),
			       NULL);
}
