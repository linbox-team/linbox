/* -*- mode: c; style: linux -*- */

/* widget-class-name.c
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

#include "widget-class-name.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _WidgetClassNamePrivate 
{
	/* Private data members */
};

static ParentClassNameClass *parent_class;

static void widget_class_name_init        (WidgetClassName *widget_class_name);
static void widget_class_name_class_init  (WidgetClassNameClass *class);

static void widget_class_name_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void widget_class_name_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void widget_class_name_finalize    (GtkObject *object);

guint
widget_class_name_get_type (void)
{
	static guint widget_class_name_type = 0;

	if (!widget_class_name_type) {
		GtkTypeInfo widget_class_name_info = {
			"WidgetClassName",
			sizeof (WidgetClassName),
			sizeof (WidgetClassNameClass),
			(GtkClassInitFunc) widget_class_name_class_init,
			(GtkObjectInitFunc) widget_class_name_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		widget_class_name_type = 
			gtk_type_unique (parent_class_name_get_type (), 
					 &widget_class_name_info);
	}

	return widget_class_name_type;
}

static void
widget_class_name_init (WidgetClassName *widget_class_name)
{
	widget_class_name->p = g_new0 (WidgetClassNamePrivate, 1);
}

static void
widget_class_name_class_init (WidgetClassNameClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("WidgetClassName::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = widget_class_name_finalize;
	object_class->set_arg = widget_class_name_set_arg;
	object_class->get_arg = widget_class_name_get_arg;

	parent_class = PARENT_CLASS_NAME_CLASS
		(gtk_type_class (parent_class_name_get_type ()));
}

static void
widget_class_name_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	WidgetClassName *widget_class_name;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WIDGET_CLASS_NAME (object));

	widget_class_name = WIDGET_CLASS_NAME (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
widget_class_name_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	WidgetClassName *widget_class_name;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WIDGET_CLASS_NAME (object));

	widget_class_name = WIDGET_CLASS_NAME (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
widget_class_name_finalize (GtkObject *object) 
{
	WidgetClassName *widget_class_name;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_WIDGET_CLASS_NAME (object));

	widget_class_name = WIDGET_CLASS_NAME (object);

	g_free (widget_class_name->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
widget_class_name_new (void) 
{
	return gtk_object_new (widget_class_name_get_type (),
			       NULL);
}
