/* -*- mode: c; style: linux -*- */

/* identifier.c
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

#include "identifier.h"
#include "glyph-layout.h"

enum {
	ARG_0,
	ARG_STRING
};

struct _IdentifierPrivate 
{
	char *string;
};

static MathUnitClass *parent_class;

static void identifier_init        (Identifier *identifier);
static void identifier_class_init  (IdentifierClass *class);

static void identifier_set_arg     (GtkObject *object, 
				    GtkArg *arg, 
				    guint arg_id);
static void identifier_get_arg     (GtkObject *object, 
				    GtkArg *arg, 
				    guint arg_id);

static void identifier_finalize    (GtkObject *object);

static Layout *identifier_get_layout (MathObject *math_object);

guint
identifier_get_type (void)
{
	static guint identifier_type = 0;

	if (!identifier_type) {
		GtkTypeInfo identifier_info = {
			"Identifier",
			sizeof (Identifier),
			sizeof (IdentifierClass),
			(GtkClassInitFunc) identifier_class_init,
			(GtkObjectInitFunc) identifier_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		identifier_type = 
			gtk_type_unique (math_unit_get_type (), 
					 &identifier_info);
	}

	return identifier_type;
}

static void
identifier_init (Identifier *identifier)
{
	identifier->p = g_new0 (IdentifierPrivate, 1);
}

static void
identifier_class_init (IdentifierClass *class) 
{
	GtkObjectClass *object_class;
	MathObjectClass *math_object_class;

	gtk_object_add_arg_type ("Identifier::string",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_STRING);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = identifier_finalize;
	object_class->set_arg = identifier_set_arg;
	object_class->get_arg = identifier_get_arg;

	math_object_class = MATH_OBJECT_CLASS (class);
	math_object_class->get_layout = identifier_get_layout;

	parent_class = MATH_UNIT_CLASS
		(gtk_type_class (math_unit_get_type ()));
}

static void
identifier_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Identifier *identifier;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_IDENTIFIER (object));

	identifier = IDENTIFIER (object);

	switch (arg_id) {
	case ARG_STRING:
		if (identifier->p->string != NULL)
			g_free (identifier->p->string);

		identifier->p->string = g_strdup (GTK_VALUE_POINTER (*arg));

		gtk_signal_emit_by_name (object, "changed", NULL);
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
identifier_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	Identifier *identifier;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_IDENTIFIER (object));

	identifier = IDENTIFIER (object);

	switch (arg_id) {
	case ARG_STRING:
		GTK_VALUE_POINTER (*arg) = identifier->p->string;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
identifier_finalize (GtkObject *object) 
{
	Identifier *identifier;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_IDENTIFIER (object));

	identifier = IDENTIFIER (object);

	g_free (identifier->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkObject *
identifier_new (gchar *string) 
{
	return gtk_object_new (identifier_get_type (),
			       "string", string,
			       NULL);
}

void
identifier_set_string (Identifier *identifier, const gchar *string)
{
	g_return_if_fail (identifier != NULL);
	g_return_if_fail (IS_IDENTIFIER (identifier));

	gtk_object_set (GTK_OBJECT (identifier), "string", string, NULL);
}

const gchar *
identifier_get_string (Identifier *identifier)
{
	g_return_val_if_fail (identifier != NULL, NULL);
	g_return_val_if_fail (IS_IDENTIFIER (identifier), NULL);

	return identifier->p->string;
}

static Layout *
identifier_get_layout (MathObject *math_object)
{
	return LAYOUT (glyph_layout_new ());
}

