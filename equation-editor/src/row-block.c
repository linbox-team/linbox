/* -*- mode: c; style: linux -*- */

/* row-block.c
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

#include "row-block.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _RowBlockPrivate 
{
	/* Private data members */
};

static BlockClass *parent_class;

static void row_block_init        (RowBlock *row_block);
static void row_block_class_init  (RowBlockClass *class);

static void row_block_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void row_block_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void row_block_finalize    (GtkObject *object);

guint
row_block_get_type (void)
{
	static guint row_block_type = 0;

	if (!row_block_type) {
		GtkTypeInfo row_block_info = {
			"RowBlock",
			sizeof (RowBlock),
			sizeof (RowBlockClass),
			(GtkClassInitFunc) row_block_class_init,
			(GtkObjectInitFunc) row_block_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		row_block_type = 
			gtk_type_unique (block_get_type (), 
					 &row_block_info);
	}

	return row_block_type;
}

static void
row_block_init (RowBlock *row_block)
{
	row_block->p = g_new0 (RowBlockPrivate, 1);
}

static void
row_block_class_init (RowBlockClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("RowBlock::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = row_block_finalize;
	object_class->set_arg = row_block_set_arg;
	object_class->get_arg = row_block_get_arg;

	parent_class = BLOCK_CLASS
		(gtk_type_class (block_get_type ()));
}

static void
row_block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	RowBlock *row_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK (object));

	row_block = ROW_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
row_block_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	RowBlock *row_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK (object));

	row_block = ROW_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
row_block_finalize (GtkObject *object) 
{
	RowBlock *row_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_ROW_BLOCK (object));

	row_block = ROW_BLOCK (object);

	g_free (row_block->p);
}

GtkObject *
row_block_new (void) 
{
	return gtk_object_new (row_block_get_type (),
			       NULL);
}
