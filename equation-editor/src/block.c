/* -*- mode: c; style: linux -*- */

/* block.c
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
 *
 * Block class: Generic multi-object containers
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "block.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _BlockPrivate 
{
};

static MathObjectClass *parent_class;

static void block_init         (Block *block);
static void block_class_init   (BlockClass *class);

static void block_set_arg      (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);
static void block_get_arg      (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);

static void block_finalize     (GtkObject *object);

static void block_real_foreach (Block *block,
				BlockIteratorCB callback,
				gpointer data);

/**
 * block_get_type
 *
 * Return type identifier and register if necessary; see Gtk+ docs for details
 */

guint
block_get_type (void)
{
	static guint block_type = 0;

	if (!block_type) {
		GtkTypeInfo block_info = {
			"Block",
			sizeof (Block),
			sizeof (BlockClass),
			(GtkClassInitFunc) block_class_init,
			(GtkObjectInitFunc) block_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		block_type = 
			gtk_type_unique (math_object_get_type (), 
					 &block_info);
	}

	return block_type;
}

/**
 * block_init
 *
 * Instance initialization function; see Gtk+ docs for details
 */

static void
block_init (Block *block)
{
	block->p = g_new0 (BlockPrivate, 1);
}

/**
 * block_class_init
 *
 * Class initialization function; see Gtk+ docs for details
 */

static void
block_class_init (BlockClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Block::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = block_finalize;
	object_class->set_arg = block_set_arg;
	object_class->get_arg = block_get_arg;

	parent_class = MATH_OBJECT_CLASS
		(gtk_type_class (math_object_get_type ()));

	class->foreach = block_real_foreach;
}

/**
 * block_set_arg
 *
 * Argument set function; see Gtk+ docs for details
 */

static void
block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Block *block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_BLOCK (object));

	block = BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

/**
 * block_get_arg
 *
 * Argument get function; see Gtk+ docs for details
 */

static void
block_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Block *block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_BLOCK (object));

	block = BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

/**
 * block_finalize
 *
 * Implementation of gtk_object_finalize
 */

static void
block_finalize (GtkObject *object) 
{
	Block *block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_BLOCK (object));

	block = BLOCK (object);

	g_free (block->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

/**
 * block_foreach:
 * @block: object
 * @callback: Callback to invoke on each object
 * @data: Data to pass to the callback on each invocation
 * 
 * Call the function @callback for each math object in the block, passing the
 * block, the math object, and any arbitrary data (given by the data
 * parameter) to the object. If the return value of the callback is nonzero,
 * return immediately; otherwise continue with the remaining objects.
 **/

void
block_foreach (Block *block, BlockIteratorCB callback, gpointer data) 
{
	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_BLOCK (block));
	g_return_if_fail (callback != NULL);

	BLOCK_CLASS (GTK_OBJECT (block)->klass)->foreach
		(block, callback, data);
}

/**
 * block_real_foreach:
 *
 * Default implementation of block_foreach
 **/

static void
block_real_foreach (Block *block, BlockIteratorCB callback, gpointer data) 
{
	g_warning ("Invoked pure virtual method Block::foreach");
}

