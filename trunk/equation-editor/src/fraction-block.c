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
#include "fraction-block-layout.h"

enum {
	ARG_0,
	ARG_NUMERATOR,
	ARG_DENOMINATOR
};

struct _FractionBlockPrivate 
{
	MathObject *numerator;
	MathObject *denominator;
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

static void fraction_block_finalize    (GtkObject *object);

static Layout *fraction_block_get_layout (MathObject *math_object);

static void fraction_block_foreach     (Block *block,
					BlockIteratorCB callback,
					gpointer data);

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
	fraction_block->p = g_new0 (FractionBlockPrivate, 1);
}

static void
fraction_block_class_init (FractionBlockClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;
	BlockClass *block_class;

	gtk_object_add_arg_type ("FractionBlock::numerator",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_NUMERATOR);

	gtk_object_add_arg_type ("FractionBlock::denominator",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_DENOMINATOR);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = fraction_block_finalize;
	object_class->set_arg = fraction_block_set_arg;
	object_class->get_arg = fraction_block_get_arg;

	parent_class = BLOCK_CLASS
		(gtk_type_class (block_get_type ()));

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = fraction_block_get_layout;

	block_class = BLOCK_CLASS (class);
	block_class->foreach = fraction_block_foreach;
}

static void
fraction_block_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	FractionBlock *fraction_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (object));

	fraction_block = FRACTION_BLOCK (object);

	switch (arg_id) {
	case ARG_NUMERATOR:
		g_return_if_fail ((GTK_VALUE_POINTER (*arg) == NULL) ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (fraction_block->p->numerator != NULL)
			gtk_object_unref
				(GTK_OBJECT (fraction_block->p->numerator));

		fraction_block->p->numerator = GTK_VALUE_POINTER (*arg);

		if (fraction_block->p->numerator != NULL)
			gtk_object_ref
				(GTK_OBJECT (fraction_block->p->numerator));

		gtk_signal_emit_by_name (GTK_OBJECT (fraction_block),
					 "changed", NULL);
		break;

	case ARG_DENOMINATOR:
		g_return_if_fail ((GTK_VALUE_POINTER (*arg) == NULL) ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (fraction_block->p->denominator != NULL)
			gtk_object_unref
				(GTK_OBJECT (fraction_block->p->denominator));

		fraction_block->p->denominator = GTK_VALUE_POINTER (*arg);

		if (fraction_block->p->denominator != NULL)
			gtk_object_ref
				(GTK_OBJECT (fraction_block->p->denominator));

		gtk_signal_emit_by_name (GTK_OBJECT (fraction_block),
					 "changed", NULL);
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
	case ARG_NUMERATOR:
		GTK_VALUE_POINTER (*arg) = fraction_block->p->numerator;
		break;

	case ARG_DENOMINATOR:
		GTK_VALUE_POINTER (*arg) = fraction_block->p->denominator;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
fraction_block_finalize (GtkObject *object) 
{
	FractionBlock *fraction_block;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (object));

	fraction_block = FRACTION_BLOCK (object);

	g_free (fraction_block->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

/**
 * fraction_block_new:
 * @numerator: 
 * @denominator: 
 * 
 * Factory method
 * 
 * Return value: 
 **/

GtkObject *
fraction_block_new (MathObject *numerator, MathObject *denominator) 
{
	return gtk_object_new (fraction_block_get_type (),
			       "numerator", numerator,
			       "denominator", denominator,
			       NULL);
}

/**
 * fraction_block_get_numerator:
 * @block: 
 * 
 * Get the numerator for the fraction
 * 
 * Return value: 
 **/

MathObject *
fraction_block_get_numerator (FractionBlock *block)
{
	g_return_val_if_fail (block != NULL, NULL);
	g_return_val_if_fail (IS_FRACTION_BLOCK (block), NULL);

	return block->p->numerator;
}

/**
 * fraction_block_get_denominator:
 * @block: 
 * 
 * Get the denominator for the fraction
 * 
 * Return value: 
 **/

MathObject *
fraction_block_get_denominator (FractionBlock *block)
{
	g_return_val_if_fail (block != NULL, NULL);
	g_return_val_if_fail (IS_FRACTION_BLOCK (block), NULL);

	return block->p->denominator;
}

/**
 * fraction_block_set_numerator:
 * @block: 
 * @numerator: 
 * 
 * Set the numerator for the fraction; wrapper for gtk_object_set
 **/

void
fraction_block_set_numerator (FractionBlock *block, MathObject *numerator)
{
	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (block));
	g_return_if_fail (numerator == NULL || IS_MATH_OBJECT (numerator));

	gtk_object_set (GTK_OBJECT (block), "numerator", numerator, NULL);
}

/**
 * fraction_block_set_denominator:
 * @block: 
 * @denominator: 
 * 
 * Set the denominator for the fraction; wrapper for gtk_object_set
 **/

void
fraction_block_set_denominator (FractionBlock *block, MathObject *denominator)
{
	g_return_if_fail (block != NULL);
	g_return_if_fail (IS_FRACTION_BLOCK (block));
	g_return_if_fail (denominator == NULL || IS_MATH_OBJECT (denominator));

	gtk_object_set (GTK_OBJECT (block), "denominator", denominator, NULL);
}

static Layout *
fraction_block_get_layout (MathObject *math_object)
{
	return LAYOUT (fraction_block_layout_new ());
}

static void
fraction_block_foreach (Block *block, BlockIteratorCB callback, gpointer data)
{
	FractionBlock *fraction_block;

	fraction_block = FRACTION_BLOCK (block);

	if (callback (block, fraction_block->p->numerator, data)) return;
	callback (block, fraction_block->p->denominator, data);
}
