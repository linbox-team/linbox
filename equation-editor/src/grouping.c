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

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _GroupingPrivate 
{
	/* Private data members */
};

static UnitClass *parent_class;

static void grouping_init        (Grouping *grouping);
static void grouping_class_init  (GroupingClass *class);

static void grouping_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void grouping_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void grouping_finalize    (GtkObject *object);

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
			gtk_type_unique (unit_get_type (), 
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

	gtk_object_add_arg_type ("Grouping::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = grouping_finalize;
	object_class->set_arg = grouping_set_arg;
	object_class->get_arg = grouping_get_arg;

	parent_class = UNIT_CLASS
		(gtk_type_class (unit_get_type ()));
}

static void
grouping_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Grouping *grouping;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPING (object));

	grouping = GROUPING (object);

	switch (arg_id) {
	case ARG_SAMPLE:
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
	case ARG_SAMPLE:
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
}

GtkObject *
grouping_new (void) 
{
	return gtk_object_new (grouping_get_type (),
			       NULL);
}
