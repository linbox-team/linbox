/* -*- mode: c; style: linux -*- */

/* unit-layout.c
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

#include "unit-layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _UnitLayoutPrivate 
{
	/* Private data members */
};

static LayoutClass *parent_class;

static void unit_layout_init        (UnitLayout *unit_layout);
static void unit_layout_class_init  (UnitLayoutClass *class);

static void unit_layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void unit_layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void unit_layout_finalize    (GtkObject *object);

guint
unit_layout_get_type (void)
{
	static guint unit_layout_type = 0;

	if (!unit_layout_type) {
		GtkTypeInfo unit_layout_info = {
			"UnitLayout",
			sizeof (UnitLayout),
			sizeof (UnitLayoutClass),
			(GtkClassInitFunc) unit_layout_class_init,
			(GtkObjectInitFunc) unit_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		unit_layout_type = 
			gtk_type_unique (layout_get_type (), 
					 &unit_layout_info);
	}

	return unit_layout_type;
}

static void
unit_layout_init (UnitLayout *unit_layout)
{
	unit_layout->p = g_new0 (UnitLayoutPrivate, 1);
}

static void
unit_layout_class_init (UnitLayoutClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("UnitLayout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = unit_layout_finalize;
	object_class->set_arg = unit_layout_set_arg;
	object_class->get_arg = unit_layout_get_arg;

	parent_class = LAYOUT_CLASS
		(gtk_type_class (layout_get_type ()));
}

static void
unit_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	UnitLayout *unit_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT_LAYOUT (object));

	unit_layout = UNIT_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
unit_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	UnitLayout *unit_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT_LAYOUT (object));

	unit_layout = UNIT_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
unit_layout_finalize (GtkObject *object) 
{
	UnitLayout *unit_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_UNIT_LAYOUT (object));

	unit_layout = UNIT_LAYOUT (object);

	g_free (unit_layout->p);
}

GtkObject *
unit_layout_new (void) 
{
	return gtk_object_new (unit_layout_get_type (),
			       NULL);
}
