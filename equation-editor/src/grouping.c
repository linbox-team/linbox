/* -*- mode: c; style: linux -*- */

/* grouping.c
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

#include "grouping.h"
/*  #include "grouping-layout.h" */

enum {
	ARG_0,
	ARG_TYPE,
	ARG_MATH_OBJECT
};

struct _GroupingPrivate 
{
	GroupingType type;
	MathObject *math_object;
};

static MathUnitClass *parent_class;

static void grouping_init        (Grouping *grouping);
static void grouping_class_init  (GroupingClass *class);

static void grouping_set_arg     (GtkObject *object, 
				  GtkArg *arg, 
				  guint arg_id);
static void grouping_get_arg     (GtkObject *object, 
				  GtkArg *arg, 
				  guint arg_id);

static void grouping_finalize    (GtkObject *object);

/*  static Layout *grouping_get_layout (MathObject *math_object); */

guint
grouping_get_type (void)
{
	static guint grouping_type = 0;

	if (!grouping_type) {
		GtkTypeInfo grouping_info = {
			"Grouping",
			sizeof (Grouping),
			sizeof (GroupingClass),
			(GtkClassInitFunc) grouping_class_init,
			(GtkObjectInitFunc) grouping_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		grouping_type = 
			gtk_type_unique (math_unit_get_type (), 
					 &grouping_info);
	}

	return grouping_type;
}

static void
grouping_init (Grouping *grouping)
{
	grouping->p = g_new0 (GroupingPrivate, 1);
}

static void
grouping_class_init (GroupingClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;

	gtk_object_add_arg_type ("Grouping::type",
				 GTK_TYPE_INT,
				 GTK_ARG_READWRITE,
				 ARG_TYPE);
	gtk_object_add_arg_type ("Grouping::math-object",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_MATH_OBJECT);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = grouping_finalize;
	object_class->set_arg = grouping_set_arg;
	object_class->get_arg = grouping_get_arg;

	math_object_class = MATH_OBJECT_CLASS (class);
/*  	math_object_class->get_layout = grouping_get_layout; */

	parent_class = MATH_UNIT_CLASS
		(gtk_type_class (math_unit_get_type ()));
}

static void
grouping_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Grouping *grouping;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPING (object));

	grouping = GROUPING (object);

	switch (arg_id) {
	case ARG_TYPE:
		grouping->p->type = GTK_VALUE_INT (*arg);
		gtk_signal_emit_by_name (GTK_OBJECT (grouping),
					 "changed", NULL);
		break;

	case ARG_MATH_OBJECT:
		g_return_if_fail ((GTK_VALUE_POINTER (*arg) == NULL) ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (grouping->p->math_object != NULL)
			gtk_object_unref
				(GTK_OBJECT (grouping->p->math_object));

		grouping->p->math_object = GTK_VALUE_POINTER (*arg);

		if (grouping->p->math_object != NULL)
			gtk_object_ref
				(GTK_OBJECT (grouping->p->math_object));

		gtk_signal_emit_by_name (GTK_OBJECT (grouping),
					 "changed", NULL);
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
grouping_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Grouping *grouping;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPING (object));

	grouping = GROUPING (object);

	switch (arg_id) {
	case ARG_TYPE:
		GTK_VALUE_INT (*arg) = grouping->p->type;
		break;

	case ARG_MATH_OBJECT:
		GTK_VALUE_POINTER (*arg) = grouping->p->math_object;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
grouping_finalize (GtkObject *object) 
{
	Grouping *grouping;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPING (object));

	grouping = GROUPING (object);

	g_free (grouping->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
grouping_new (GroupingType type, MathObject *math_object) 
{
	return gtk_object_new (grouping_get_type (),
			       "type", type,
			       "math-object", math_object,
			       NULL);
}

void 
grouping_set_math_object (Grouping *grouping, MathObject *math_object)
{
	g_return_if_fail (grouping != NULL);
	g_return_if_fail (IS_GROUPING (grouping));
	g_return_if_fail (math_object == NULL || IS_MATH_OBJECT (math_object));

	gtk_object_set (GTK_OBJECT (grouping),
			"math-object", math_object,
			NULL);
}

MathObject *
grouping_get_math_object (Grouping *grouping)
{
	g_return_val_if_fail (grouping != NULL, NULL);
	g_return_val_if_fail (IS_GROUPING (grouping), NULL);

	return grouping->p->math_object;
}

void
grouping_set_grouping_type (Grouping *grouping, GroupingType type)
{
	g_return_if_fail (grouping != NULL);
	g_return_if_fail (IS_GROUPING (grouping));

	gtk_object_set (GTK_OBJECT (grouping),
			"type", type,
			NULL);
}

GroupingType
grouping_get_grouping_type (Grouping *grouping)
{
	g_return_val_if_fail (grouping != NULL, GROUPING_PAREN);
	g_return_val_if_fail (IS_GROUPING (grouping), GROUPING_PAREN);

	return grouping->p->type;
}

/*  static Layout * */
/*  grouping_get_layout (MathObject *math_object) */
/*  { */
/*  	return LAYOUT (grouping_layout_new ()); */
/*  } */


