/* -*- mode: c; style: linux -*- */

/* controller.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 *            Rob Wehde, Matt Spilich, Anthony Asher
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

#include "controller.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _ControllerPrivate 
{
     /* Cursor current_pos; */
};

static GtkObjectClass *parent_class;

static void controller_init        (Controller *controller);
static void controller_class_init  (ControllerClass *class);

static void controller_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void controller_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void controller_finalize    (GtkObject *object);



guint
controller_get_type (void)
{
	static guint controller_type = 0;

	if (!controller_type) {
		GtkTypeInfo controller_info = {
			"Controller",
			sizeof (Controller),
			sizeof (ControllerClass),
			(GtkClassInitFunc) controller_class_init,
			(GtkObjectInitFunc) controller_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		controller_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &controller_info);
	}

	return controller_type;
}

static void
controller_init (Controller *controller)
{
	controller->p = g_new0 (ControllerPrivate, 1);
}

static void
controller_class_init (ControllerClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Controller::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = controller_finalize;
	object_class->set_arg = controller_set_arg;
	object_class->get_arg = controller_get_arg;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
controller_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Controller *controller;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CONTROLLER (object));

	controller = CONTROLLER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
controller_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Controller *controller;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CONTROLLER (object));

	controller = CONTROLLER (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
controller_finalize (GtkObject *object) 
{
	Controller *controller;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_CONTROLLER (object));

	controller = CONTROLLER (object);

	g_free (controller->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
controller_new (void) 
{
	return gtk_object_new (controller_get_type (),
			       NULL);
}



/**Notethis only works for the key #1 being pressed.  The rest are easy to add
once we get the 1 working.**/

static void controller_insert(GdkEventKey *event, MathObject *obj, int pos) 
{

	gchar *this_keypressed;

	g_return_if_fail (obj != NULL);
	g_return_if_fail (IS_MATH_OBJECT (obj));

	this_keypressed = event->string;

	g_warning(*this_keypressed);
	
}
