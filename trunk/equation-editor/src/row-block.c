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
#include "row-block-layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _RowBlockPrivate 
{
	GList *objects;
};

static BlockClass *parent_class;

static RowBlockLayout *layout;

static void row_block_init        (RowBlock *row_block);
static void row_block_class_init  (RowBlockClass *class);

static void row_block_set_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);
static void row_block_get_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);

static void row_block_finalize    (GtkObject *object);

static const Layout *row_block_get_layout (MathObject *math_object);

static void row_block_foreach     (Block *block,
				   BlockIteratorCB callback,
				   gpointer data);

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
	MathObjectClass *math_object_class;
	BlockClass *block_class;

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

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = row_block_get_layout;

	block_class = BLOCK_CLASS (class);
	block_class->foreach = row_block_foreach;

	if (layout == NULL)
		layout = ROW_BLOCK_LAYOUT (row_block_layout_new ());
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

/**
 * row_block_new:
 * @void: 
 * 
 * Factory method
 * 
 * Return value: 
 **/

GtkObject *
row_block_new (void) 
{
	return gtk_object_new (row_block_get_type (),
			       NULL);
}

/**
 * row_block_insert:
 * @row_block: 
 * @math_object: 
 * @before: The object before which to insert, NULL to insert at the end
 * 
 * Insert the given math object before the object `before'; append to the end
 * if the `before' is not located in the row block
 **/

void 
row_block_insert (RowBlock *row_block, MathObject *math_object,
		  MathObject *before)
{
	gint position;

	g_return_if_fail (row_block != NULL);
	g_return_if_fail (IS_ROW_BLOCK (row_block));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));
	g_return_if_fail (before == NULL || IS_MATH_OBJECT (before));

	position = g_list_index (row_block->p->objects, before);
	row_block->p->objects = g_list_insert (row_block->p->objects,
					       math_object, position);
}

/**
 * row_block_insert_at:
 * @row_block: 
 * @math_object: 
 * @position: 
 * 
 * Insert a math object at the given position
 **/

void 
row_block_insert_at (RowBlock *row_block, MathObject *math_object,
		     gint position)
{
	g_return_if_fail (row_block != NULL);
	g_return_if_fail (IS_ROW_BLOCK (row_block));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	row_block->p->objects = g_list_insert (row_block->p->objects,
					       math_object, position);
}

/**
 * row_block_get_object_at:
 * @row_block: 
 * @position: Position at which to get object (0 <= position < length)
 * 
 * Get the object at the specified position
 * 
 * Return value: 
 **/

MathObject *
row_block_get_object_at (RowBlock *row_block, gint position)
{
	g_return_val_if_fail (row_block != NULL, NULL);
	g_return_val_if_fail (IS_ROW_BLOCK (row_block), NULL);
	g_return_val_if_fail (position > 
			      g_list_length (row_block->p->objects), NULL);

	return g_list_nth_data (row_block->p->objects, position);
}

/**
 * row_block_get_length:
 * @row_block: 
 * 
 * Get the number of objects in the row block
 * 
 * Return value: 
 **/

guint
row_block_get_length (RowBlock *row_block) 
{
	g_return_val_if_fail (row_block != NULL, 0);
	g_return_val_if_fail (IS_ROW_BLOCK (row_block), 0);

	return g_list_length (row_block->p->objects);
}

static const Layout *
row_block_get_layout (MathObject *math_object) 
{
	return LAYOUT (layout);
}

static void
row_block_foreach (Block *block, BlockIteratorCB callback, gpointer data) 
{
	RowBlock *row_block;
	GList *node;

	row_block = ROW_BLOCK (block);

	for (node = row_block->p->objects; node; node = node->next)
		if (callback (block, MATH_OBJECT (node->data), data)) break;
}
