/* -*- mode: c; style: linux -*- */

/* unit.c
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

#include "unit.h"

enum {
	ARG_0,
	ARG_SUPERSCRIPT,
	ARG_SUBSCRIPT
};

struct _UnitPrivate 
{
	MathObject *superscript;
	MathObject *subscript;
};

static MathObjectClass *parent_class;

static void unit_init        (Unit *unit);
static void unit_class_init  (UnitClass *class);

static void unit_set_arg     (GtkObject *object, 
			      GtkArg *arg, 
			      guint arg_id);
static void unit_get_arg     (GtkObject *object, 
			      GtkArg *arg, 
			      guint arg_id);

static void unit_finalize    (GtkObject *object);

guint
unit_get_type (void)
{
	static guint unit_type = 0;

	if (!unit_type) {
		GtkTypeInfo unit_info = {
			"Unit",
			sizeof (Unit),
			sizeof (UnitClass),
			(GtkClassInitFunc) unit_class_init,
			(GtkObjectInitFunc) unit_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		unit_type = 
			gtk_type_unique (math_object_get_type (), 
					 &unit_info);
	}

	return unit_type;
}

static void
unit_init (Unit *unit)
{
	unit->p = g_new0 (UnitPrivate, 1);
}

static void
unit_class_init (UnitClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Unit::superscript",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SUPERSCRIPT);
	gtk_object_add_arg_type ("Unit::subscript",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SUBSCRIPT);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = unit_finalize;
	object_class->set_arg = unit_set_arg;
	object_class->get_arg = unit_get_arg;

	parent_class = MATH_OBJECT_CLASS
		(gtk_type_class (math_object_get_type ()));
}

static void
unit_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Unit *unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT (object));

	unit = UNIT (object);

	switch (arg_id) {
	case ARG_SUPERSCRIPT:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (unit->p->superscript != NULL)
			gtk_object_unref (GTK_OBJECT (unit->p->superscript));

		unit->p->superscript = MATH_OBJECT (GTK_VALUE_POINTER (*arg));

		if (unit->p->superscript != NULL)
			gtk_object_ref (GTK_OBJECT (unit->p->superscript));
		break;

	case ARG_SUBSCRIPT:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (unit->p->subscript != NULL)
			gtk_object_unref (GTK_OBJECT (unit->p->subscript));

		unit->p->subscript = MATH_OBJECT (GTK_VALUE_POINTER (*arg));

		if (unit->p->subscript != NULL)
			gtk_object_ref (GTK_OBJECT (unit->p->subscript));
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
unit_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Unit *unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT (object));

	unit = UNIT (object);

	switch (arg_id) {
	case ARG_SUPERSCRIPT:
		GTK_VALUE_POINTER (*arg) = unit->p->superscript;
		break;

	case ARG_SUBSCRIPT:
		GTK_VALUE_POINTER (*arg) = unit->p->subscript;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
unit_finalize (GtkObject *object) 
{
	Unit *unit;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT (object));

	unit = UNIT (object);

	if (unit->p->superscript != NULL)
		gtk_object_unref (unit->p->superscript);

	if (unit->p->subscript != NULL)
		gtk_object_unref (unit->p->subscript);

	g_free (unit->p);
}

/**
 * unit_get_superscript:
 * @unit: 
 * 
 * Get this unit's superscript
 * 
 * Return value: 
 **/

MathObject *
unit_get_superscript (Unit *unit)
{
	g_return_val_if_fail (unit != NULL, NULL);
	g_return_val_if_fail (IS_UNIT (unit), NULL);

	return unit->p->superscript;
}

/**
 * unit_get_subscript:
 * @unit: 
 * 
 * Get this unit's subscript
 * 
 * Return value: 
 **/

MathObject *
unit_get_subscript (Unit *unit)
{
	g_return_val_if_fail (unit != NULL, NULL);
	g_return_val_if_fail (IS_UNIT (unit), NULL);

	return unit->p->subscript;
}

/**
 * unit_set_superscript:
 * @unit: 
 * @math_object: 
 * 
 * Set this unit's superscript to a new object; wrapper for gtk_object_set
 **/

void
unit_set_superscript (Unit *unit, MathObject *math_object)
{
	g_return_if_fail (unit != NULL);
	g_return_if_fail (IS_UNIT (unit));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	gtk_object_set (GTK_OBJECT (unit), "superscript", math_object, NULL);
}

/**
 * unit_set_subscript:
 * @unit: 
 * @math_object: 
 * 
 * Set this unit's subscript to a new value, wrapper for gtk_object_set
 **/

void
unit_set_subscript (Unit *unit, MathObject *math_object)
{
	g_return_if_fail (unit != NULL);
	g_return_if_fail (IS_UNIT (unit));
	g_return_if_fail (math_object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (math_object));

	gtk_object_set (GTK_OBJECT (unit), "subscript", math_object, NULL);
}

