/* -*- mode: c; style: linux -*- */

/* math-unit.c
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

#include "math-unit.h"

enum {
	ARG_0,
	ARG_SUPERSCRIPT,
	ARG_SUBSCRIPT
};

struct _MathUnitPrivate 
{
	MathObject *superscript;
	MathObject *subscript;
};

static MathObjectClass *parent_class;

static void math_unit_init        (MathUnit *math_unit);
static void math_unit_class_init  (MathUnitClass *class);

static void math_unit_set_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);
static void math_unit_get_arg     (GtkObject *object, 
				   GtkArg *arg, 
				   guint arg_id);

static void math_unit_finalize    (GtkObject *object);

guint
math_unit_get_type (void)
{
	static guint math_unit_type = 0;

	if (!math_unit_type) {
		GtkTypeInfo math_unit_info = {
			"MathUnit",
			sizeof (MathUnit),
			sizeof (MathUnitClass),
			(GtkClassInitFunc) math_unit_class_init,
			(GtkObjectInitFunc) math_unit_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_unit_type = 
			gtk_type_unique (math_object_get_type (), 
					 &math_unit_info);
	}

	return math_unit_type;
}

static void
math_unit_init (MathUnit *math_unit)
{
	math_unit->p = g_new0 (MathUnitPrivate, 1);
}

static void
math_unit_class_init (MathUnitClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MathUnit::superscript",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SUPERSCRIPT);
	gtk_object_add_arg_type ("MathUnit::subscript",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SUBSCRIPT);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = math_unit_finalize;
	object_class->set_arg = math_unit_set_arg;
	object_class->get_arg = math_unit_get_arg;

	parent_class = MATH_OBJECT_CLASS
		(gtk_type_class (math_object_get_type ()));
}

static void
math_unit_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathUnit *math_unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_UNIT (object));

	math_unit = MATH_UNIT (object);

	switch (arg_id) {
	case ARG_SUPERSCRIPT:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (math_unit->p->superscript != NULL)
			gtk_object_unref
				(GTK_OBJECT (math_unit->p->superscript));

		math_unit->p->superscript = 
			MATH_OBJECT (GTK_VALUE_POINTER (*arg));

		if (math_unit->p->superscript != NULL)
			gtk_object_ref
				(GTK_OBJECT (math_unit->p->superscript));
		break;

	case ARG_SUBSCRIPT:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (math_unit->p->subscript != NULL)
			gtk_object_unref
				(GTK_OBJECT (math_unit->p->subscript));

		math_unit->p->subscript = 
			MATH_OBJECT (GTK_VALUE_POINTER (*arg));

		if (math_unit->p->subscript != NULL)
			gtk_object_ref (GTK_OBJECT (math_unit->p->subscript));
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
math_unit_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathUnit *math_unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_UNIT (object));

	math_unit = MATH_UNIT (object);

	switch (arg_id) {
	case ARG_SUPERSCRIPT:
		GTK_VALUE_POINTER (*arg) = math_unit->p->superscript;
		break;

	case ARG_SUBSCRIPT:
		GTK_VALUE_POINTER (*arg) = math_unit->p->subscript;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
math_unit_finalize (GtkObject *object) 
{
	MathUnit *math_unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_UNIT (object));

	math_unit = MATH_UNIT (object);

	if (math_unit->p->superscript != NULL)
		gtk_object_unref (GTK_OBJECT (math_unit->p->superscript));

	if (math_unit->p->subscript != NULL)
		gtk_object_unref (GTK_OBJECT (math_unit->p->subscript));

	g_free (math_unit->p);
}

/**
 * math_unit_get_superscript:
 * @math_unit: 
 * 
 * Get this math_unit's superscript
 * 
 * Return value: 
 **/

MathObject *
math_unit_get_superscript (MathUnit *math_unit)
{
	g_return_val_if_fail (math_unit != NULL, NULL);
	g_return_val_if_fail (IS_MATH_UNIT (math_unit), NULL);

	return math_unit->p->superscript;
}

/**
 * math_unit_get_subscript:
 * @math_unit: 
 * 
 * Get this math_unit's subscript
 * 
 * Return value: 
 **/

MathObject *
math_unit_get_subscript (MathUnit *math_unit)
{
	g_return_val_if_fail (math_unit != NULL, NULL);
	g_return_val_if_fail (IS_MATH_UNIT (math_unit), NULL);

	return math_unit->p->subscript;
}

/**
 * math_unit_set_superscript:
 * @math_unit: 
 * @math_object: 
 * 
 * Set this math_unit's superscript to a new object; wrapper for gtk_object_set
 **/

void
math_unit_set_superscript (MathUnit *math_unit, MathObject *math_object)
{
	g_return_if_fail (math_unit != NULL);
	g_return_if_fail (IS_MATH_UNIT (math_unit));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	gtk_object_set (GTK_OBJECT (math_unit), "superscript", math_object, NULL);
}

/**
 * math_unit_set_subscript:
 * @math_unit: 
 * @math_object: 
 * 
 * Set this math_unit's subscript to a new value, wrapper for gtk_object_set
 **/

void
math_unit_set_subscript (MathUnit *math_unit, MathObject *math_object)
{
	g_return_if_fail (math_unit != NULL);
	g_return_if_fail (IS_MATH_UNIT (math_unit));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	gtk_object_set (GTK_OBJECT (math_unit), "subscript", math_object, NULL);
}

