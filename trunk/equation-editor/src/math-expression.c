/* -*- mode: c; style: linux -*- */

/* math-expression.c
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

#include "math-expression.h"

enum {
	ARG_0,
	ARG_TOPLEVEL
};

struct _MathExpressionPrivate 
{
	MathObject *toplevel;
};

static GtkObjectClass *parent_class;

static void math_expression_init        (MathExpression *math_expression);
static void math_expression_class_init  (MathExpressionClass *class);

static void math_expression_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void math_expression_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

static void math_expression_finalize    (GtkObject *object);

guint
math_expression_get_type (void)
{
	static guint math_expression_type = 0;

	if (!math_expression_type) {
		GtkTypeInfo math_expression_info = {
			"MathExpression",
			sizeof (MathExpression),
			sizeof (MathExpressionClass),
			(GtkClassInitFunc) math_expression_class_init,
			(GtkObjectInitFunc) math_expression_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_expression_type = 
			gtk_type_unique (gtk_object_get_type (), 
					 &math_expression_info);
	}

	return math_expression_type;
}

static void
math_expression_init (MathExpression *math_expression)
{
	math_expression->p = g_new0 (MathExpressionPrivate, 1);
}

static void
math_expression_class_init (MathExpressionClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MathExpression::toplevel",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_TOPLEVEL);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = math_expression_finalize;
	object_class->set_arg = math_expression_set_arg;
	object_class->get_arg = math_expression_get_arg;

	parent_class = GTK_OBJECT_CLASS
		(gtk_type_class (gtk_object_get_type ()));
}

static void
math_expression_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathExpression *math_expression;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION (object));

	math_expression = MATH_EXPRESSION (object);

	switch (arg_id) {
	case ARG_TOPLEVEL:
		g_return_if_fail ((GTK_VALUE_POINTER (*arg) == NULL) ||
				  IS_MATH_OBJECT (GTK_VALUE_POINTER (*arg)));

		if (math_expression->p->toplevel != NULL)
			gtk_object_unref 
				(GTK_OBJECT (math_expression->p->toplevel));

		math_expression->p->toplevel = 
			MATH_OBJECT (GTK_VALUE_POINTER (*arg));
		gtk_object_ref (GTK_OBJECT (math_expression->p->toplevel));

		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
math_expression_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathExpression *math_expression;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION (object));

	math_expression = MATH_EXPRESSION (object);

	switch (arg_id) {
	case ARG_TOPLEVEL:
		GTK_VALUE_POINTER (*arg) = math_expression->p->toplevel;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
math_expression_finalize (GtkObject *object) 
{
	MathExpression *math_expression;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION (object));

	math_expression = MATH_EXPRESSION (object);

	g_free (math_expression->p);
}

/**
 * math_expression_new:
 * @toplevel: 
 * 
 * Factory method
 * 
 * Return value: 
 **/

GtkObject *
math_expression_new (MathObject *toplevel) 
{
	return gtk_object_new (math_expression_get_type (),
			       "toplevel", toplevel,
			       NULL);
}

/**
 * math_expression_get_toplevel:
 * @expression: 
 * 
 * Return the top level math object
 * 
 * Return value: Pointer to top level math object; should be refed if needed
 * for long-term use
 **/

MathObject *
math_expression_get_toplevel (MathExpression *expression)
{
	g_return_if_fail (expression != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION (expression));

	return expression->p->toplevel;
}
