/* -*- mode: c; style: linux -*- */

/* math-object.c
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
 *
 * Math object class: Generic mathematical objects
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "math-object.h"
#include "layout.h"

enum {
	ARG_0,
	ARG_PARENT
};

enum {
	CHANGED_SIGNAL,
	LAST_SIGNAL
};

struct _MathObjectPrivate 
{
	MathObject *parent;
};

static GtkObjectClass *parent_class;

static gint math_object_signals[LAST_SIGNAL] = { 0 };

static void math_object_init        (MathObject *math_object);
static void math_object_class_init  (MathObjectClass *class);

static void math_object_set_arg     (GtkObject *object, 
				     GtkArg *arg, 
				     guint arg_id);
static void math_object_get_arg     (GtkObject *object, 
				     GtkArg *arg, 
				     guint arg_id);

static void math_object_finalize    (GtkObject *object);

static Layout *math_object_real_get_layout (MathObject *math_object);

/**
 * math_object_get_type
 *
 * Return type identifier and register if necessary; see Gtk+ docs for details
 */

guint
math_object_get_type (void)
{
	static guint math_object_type = 0;

	if (!math_object_type) {
		GtkTypeInfo math_object_info = {
			"MathObject",
			sizeof (MathObject),
			sizeof (MathObjectClass),
			(GtkClassInitFunc) math_object_class_init,
			(GtkObjectInitFunc) math_object_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_object_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &math_object_info);
	}

	return math_object_type;
}

/**
 * math_object_init
 *
 * Instance initialization function; see Gtk+ docs for details
 */

static void
math_object_init (MathObject *math_object)
{
	math_object->p = g_new0 (MathObjectPrivate, 1);
}

/**
 * math_object_class_init
 *
 * Class initialization function; see Gtk+ docs for details
 */

static void
math_object_class_init (MathObjectClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MathObject::parent",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_PARENT);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = math_object_finalize;
	object_class->set_arg = math_object_set_arg;
	object_class->get_arg = math_object_get_arg;

	math_object_signals[CHANGED_SIGNAL] =
		gtk_signal_new ("changed", GTK_RUN_FIRST,
				object_class->type,
				GTK_SIGNAL_OFFSET (MathObjectClass, changed),
				gtk_signal_default_marshaller,
				GTK_TYPE_NONE, 0);

	gtk_object_class_add_signals (object_class, math_object_signals,
				      LAST_SIGNAL);

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));

	class->get_layout = math_object_real_get_layout;
	class->changed = NULL;
}

/**
 * math_object_set_arg
 *
 * Argument set function; see Gtk+ docs for details
 */

static void
math_object_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathObject *math_object;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (object));

	math_object = MATH_OBJECT (object);

	switch (arg_id) {
	case ARG_PARENT:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (math_object->p->parent != NULL)
			gtk_object_unref (GTK_OBJECT (math_object->p->parent));

		math_object->p->parent = GTK_VALUE_POINTER (*arg);

		if (math_object->p->parent != NULL)
			gtk_object_ref (GTK_OBJECT (math_object->p->parent));

		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

/**
 * math_object_get_arg
 *
 * Argument get function; see Gtk+ docs for details
 */

static void
math_object_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathObject *math_object;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (object));

	math_object = MATH_OBJECT (object);

	switch (arg_id) {
	case ARG_PARENT:
		GTK_VALUE_POINTER (*arg) = math_object->p->parent;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

/**
 * math_object_finalize
 *
 * Implementation of gtk_object_finalize
 */

static void
math_object_finalize (GtkObject *object) 
{
	MathObject *math_object;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_OBJECT (object));

	math_object = MATH_OBJECT (object);

	g_free (math_object->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

/**
 * math_object_get_layout:
 * @math_object: object
 * 
 * Return a layout object used to render the math object
 * 
 * Return value: Layout object; should be unrefed when done
 **/

Layout *
math_object_get_layout (MathObject *math_object) 
{
	g_return_val_if_fail (math_object != NULL, NULL);
	g_return_val_if_fail (IS_MATH_OBJECT (math_object), NULL);

	return MATH_OBJECT_CLASS (GTK_OBJECT (math_object)->klass)->
		get_layout (math_object);
}

/**
 * math_object_get_parent:
 * @math_object: 
 * 
 * Returns a pointer to the math object's parent object
 * 
 * Return value: Parent object; should be unrefed
 **/

MathObject *
math_object_get_parent (MathObject *math_object)
{
	g_return_val_if_fail (math_object != NULL, NULL);
	g_return_val_if_fail (IS_MATH_OBJECT (math_object), NULL);

	if (math_object->p->parent != NULL)
		gtk_object_ref (GTK_OBJECT (math_object->p->parent));
	return math_object->p->parent;
}

/**
 * math_object_real_get_layout:
 *
 * Default implementation of math_object_get_layout
 **/

static Layout *
math_object_real_get_layout (MathObject *math_object) 
{
	g_warning ("Invoked pure virtual method MathObject::get_layout");
	return NULL;
}
