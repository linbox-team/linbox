/* -*- mode: c; style: linux -*- */

/* fraction-block.c
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

#include "fraction-block.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static BlockClass *parent_class;

static void fraction_block_init        (FractionBlock *fraction_block);
static void fraction_block_class_init  (FractionBlockClass *class);

static void fraction_block_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void fraction_block_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
fraction_block_get_type (void)
{
	static guint fraction_block_type = 0;

	if (!fraction_block_type) {
		GtkTypeInfo fraction_block_info = {
			"FractionBlock",
			sizeof (FractionBlock),
			sizeof (FractionBlockClass),
			(GtkClassInitFunc) fraction_block_class_init,
			(GtkObjectInitFunc) fraction_block_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		fraction_block_type = 
			gtk_type_unique (block_get_type (), 
					 &fraction_block_info);
	}

	return fraction_block_type;
}

static void
fraction_block_init (FractionBlock *fraction_block)
{
}

static void
fraction_block_class_init (FractionBlockClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("FractionBlock::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = fraction_block_set_arg;
	object_class->get_arg = fraction_block_get_arg;

	parent_class = BLOCK_CLASS
		(gtk_type_class (block_get_type ()));
}

static void
fraction_block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	FractionBlock *fraction_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (object));

	fraction_block = FRACTION_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
fraction_block_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	FractionBlock *fraction_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (object));

	fraction_block = FRACTION_BLOCK (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
fraction_block_new (void) 
{
	return gtk_object_new (fraction_block_get_type (),
			       NULL);
}
