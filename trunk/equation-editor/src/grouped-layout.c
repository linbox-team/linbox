/* -*- mode: c; style: linux -*- */

/* grouped-layout.c
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

#include "grouped-layout.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _GroupedLayoutPrivate 
{
	/* Private data members */
};

static UnitLayoutClass *parent_class;

static void grouped_layout_init        (GroupedLayout *grouped_layout);
static void grouped_layout_class_init  (GroupedLayoutClass *class);

static void grouped_layout_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void grouped_layout_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void grouped_layout_finalize    (GtkObject *object);

guint
grouped_layout_get_type (void)
{
	static guint grouped_layout_type = 0;

	if (!grouped_layout_type) {
		GtkTypeInfo grouped_layout_info = {
			"GroupedLayout",
			sizeof (GroupedLayout),
			sizeof (GroupedLayoutClass),
			(GtkClassInitFunc) grouped_layout_class_init,
			(GtkObjectInitFunc) grouped_layout_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		grouped_layout_type = 
			gtk_type_unique (unit_layout_get_type (), 
					 &grouped_layout_info);
	}

	return grouped_layout_type;
}

static void
grouped_layout_init (GroupedLayout *grouped_layout)
{
	grouped_layout->p = g_new0 (GroupedLayoutPrivate, 1);
}

static void
grouped_layout_class_init (GroupedLayoutClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("GroupedLayout::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = grouped_layout_finalize;
	object_class->set_arg = grouped_layout_set_arg;
	object_class->get_arg = grouped_layout_get_arg;

	parent_class = UNIT_LAYOUT_CLASS
		(gtk_type_class (unit_layout_get_type ()));
}

static void
grouped_layout_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	GroupedLayout *grouped_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPED_LAYOUT (object));

	grouped_layout = GROUPED_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
grouped_layout_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	GroupedLayout *grouped_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPED_LAYOUT (object));

	grouped_layout = GROUPED_LAYOUT (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
grouped_layout_finalize (GtkObject *object) 
{
	GroupedLayout *grouped_layout;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_GROUPED_LAYOUT (object));

	grouped_layout = GROUPED_LAYOUT (object);

	g_free (grouped_layout->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
grouped_layout_new (void) 
{
	return gtk_object_new (grouped_layout_get_type (),
			       NULL);
}
