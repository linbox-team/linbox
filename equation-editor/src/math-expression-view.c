/* -*- mode: c; style: linux -*- */

/* math-expression-view.c
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

#include "math-expression-view.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

static GtkWidgetClass *gtk_widget_class;

static void math_expression_view_init        (MathExpressionView *math_expression_view);
static void math_expression_view_class_init  (MathExpressionViewClass *class);

static void math_expression_view_set_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);
static void math_expression_view_get_arg     (GtkObject *object, 
					   GtkArg *arg, 
					   guint arg_id);

guint
math_expression_view_get_type (void)
{
	static guint math_expression_view_type = 0;

	if (!math_expression_view_type) {
		GtkTypeInfo math_expression_view_info = {
			"MathExpressionView",
			sizeof (MathExpressionView),
			sizeof (MathExpressionViewClass),
			(GtkClassInitFunc) math_expression_view_class_init,
			(GtkObjectInitFunc) math_expression_view_init,
			(GtkArgSetFunc) NULL,
			(GtkArgGetFunc) NULL
		};

		math_expression_view_type = 
			gtk_type_unique (gtk_widget_get_type (), 
					 &math_expression_view_info);
	}

	return math_expression_view_type;
}

static void
math_expression_view_init (MathExpressionView *math_expression_view)
{
}

static void
math_expression_view_class_init (MathExpressionViewClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("MathExpressionView::sample",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_SAMPLE);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->set_arg = math_expression_view_set_arg;
	object_class->get_arg = math_expression_view_get_arg;

	parent_class = GTK_WIDGET_CLASS
		(gtk_type_class (gtk_widget_get_type ()));
}

static void
math_expression_view_set_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (object));

	math_expression_view = MATH_EXPRESSION_VIEW (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument set");
		break;
	}
}

static void
math_expression_view_get_arg (GtkObject *object, GtkArg *arg, guint arg_id) 
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (object));

	math_expression_view = MATH_EXPRESSION_VIEW (object);

	switch (arg_id) {
	case ARG_SAMPLE:
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

GtkObject *
math_expression_view_new (void) 
{
	return gtk_object_new (math_expression_view_get_type (),
			       NULL);
}
