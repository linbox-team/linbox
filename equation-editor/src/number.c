/* -*- mode: c; style: linux -*- */

/* number.c
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

#include "number.h"
#include "glyph-layout.h"

enum {
	ARG_0,
	ARG_VALUE
};

struct _NumberPrivate 
{
	gdouble value;
};

static MathUnitClass *parent_class;

static GlyphLayout *layout;

static void number_init        (Number *number);
static void number_class_init  (NumberClass *class);

static void number_set_arg     (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);
static void number_get_arg     (GtkObject *object, 
				GtkArg *arg, 
				guint arg_id);

static void number_finalize    (GtkObject *object);

static const Layout *number_get_layout (MathObject *math_object);

guint
number_get_type (void)
{
	static guint number_type = 0;

	if (!number_type) {
		GtkTypeInfo number_info = {
			"Number",
			sizeof (Number),
			sizeof (NumberClass),
			(GtkClassInitFunc) number_class_init,
			(GtkObjectInitFunc) number_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		number_type = 
			gtk_type_unique (math_unit_get_type (), 
					 &number_info);
	}

	return number_type;
}

static void
number_init (Number *number)
{
	number->p = g_new0 (NumberPrivate, 1);
	number->p->value = 0.0;
}

static void
number_class_init (NumberClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;

	gtk_object_add_arg_type ("Number::value",
				 GTK_TYPE_FLOAT,
				 GTK_ARG_READWRITE,
				 ARG_VALUE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = number_finalize;
	object_class->set_arg = number_set_arg;
	object_class->get_arg = number_get_arg;

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = number_get_layout;

	parent_class = MATH_UNIT_CLASS
		(gtk_type_class (math_unit_get_type ()));

	layout = GLYPH_LAYOUT (glyph_layout_new ());
}

static void
number_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Number *number;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_NUMBER (object));

	number = NUMBER (object);

	switch (arg_id) {
	case ARG_VALUE:
		number->p->value = GTK_VALUE_FLOAT (*arg);
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
number_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Number *number;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_NUMBER (object));

	number = NUMBER (object);

	switch (arg_id) {
	case ARG_VALUE:
		GTK_VALUE_FLOAT (*arg) = number->p->value;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
number_finalize (GtkObject *object) 
{
	Number *number;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_NUMBER (object));

	number = NUMBER (object);

	g_free (number->p);
}

GtkObject *
number_new (gdouble value) 
{
	return gtk_object_new (number_get_type (),
			       "value", value,
			       NULL);
}

gdouble
number_get_value (Number *number)
{
	g_return_val_if_fail (number != NULL, 0.0);
	g_return_val_if_fail (IS_NUMBER (number), 0.0);

	return number->p->value;
}

void
number_set_value (Number *number, gdouble value)
{
	g_return_if_fail (number != NULL);
	g_return_if_fail (IS_NUMBER (number));

	number->p->value = value;
}

static const Layout *
number_get_layout (MathObject *math_object)
{
	return LAYOUT (glyph_layout_new ());
}

