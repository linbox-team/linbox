/* -*- mode: c; style: linux -*- */

/* matrix-block.c
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

#include "matrix-block.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static BlockClass *block_class;

static void matrix_block_init        (MatrixBlock *matrix_block);
static void matrix_block_class_init  (MatrixBlockClass *class);

static void matrix_block_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void matrix_block_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
matrix_block_get_type (void)
{
	static guint matrix_block_type = 0;

	if (!matrix_block_type) {
		GtkTypeInfo matrix_block_info = {
			"MatrixBlock",
			sizeof (MatrixBlock),
			sizeof (MatrixBlockClass),
			(GtkClassInitFunc) matrix_block_class_init,
			(GtkObjectInitFunc) matrix_block_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		matrix_block_type = 
			gtk_type_unique (block_get_type (), 
					 &matrix_block_info);
	}

	return matrix_block_type;
}

static void
matrix_block_init (MatrixBlock *matrix_block)
{
}

static void
matrix_block_class_init (MatrixBlockClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MatrixBlock::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = matrix_block_set_arg;
	object_class->get_arg = matrix_block_get_arg;

	parent_class = BLOCK_CLASS
		(gtk_type_class (block_get_type ()));
}

static void
matrix_block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlock *matrix_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (object));

	matrix_block = MATRIX_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
matrix_block_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MatrixBlock *matrix_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATRIX_BLOCK (object));

	matrix_block = MATRIX_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
matrix_block_new (void) 
{
	return gtk_object_new (matrix_block_get_type (),
			       NULL);
}
