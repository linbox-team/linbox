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
#include "renderer.h"

enum {
	ARG_0,
	ARG_EXPRESSION
};

struct _MathExpressionViewPrivate 
{
	MathExpression *expression;
	Renderer *renderer;
};

static GtkWidgetClass *parent_class;

static void math_expression_view_init        (MathExpressionView *math_expression_view);
static void math_expression_view_class_init  (MathExpressionViewClass *class);

static void math_expression_view_set_arg     (GtkObject *object, 
					      GtkArg *arg, 
					      guint arg_id);
static void math_expression_view_get_arg     (GtkObject *object, 
					      GtkArg *arg, 
					      guint arg_id);

static void math_expression_view_finalize    (GtkObject *object);

static void math_expression_view_realize     (GtkWidget *widget);
static void math_expression_view_expose      (GtkWidget *widget,
					      GdkEventExpose *event);
static void math_expression_view_key_press   (GtkWidget *widget,
					      GdkEventKey *event);
static void math_expression_view_button_press (GtkWidget *widget,
					       GdkEventButton *event);
static void math_expression_view_button_release (GtkWidget *widget,
						 GdkEventButton *event);
static void math_expression_view_selection_get (GtkWidget *widget,
						GtkSelectionData *data,
						guint info,
						guint time);
static void math_expression_view_selection_received (GtkWidget *widget,
						     GtkSelectionData *data,
						     guint time);

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
	math_expression_view->p = g_new0 (MathExpressionViewPrivate, 1);
	math_expression_view->p->renderer = 
		canvas_renderer_new (GTK_WIDGET (math_expression_view));
}

static void
math_expression_view_class_init (MathExpressionViewClass *class) 
{
	GtkObjectClass *object_class;
	GtkWidgetClass *widget_class;

	gtk_object_add_arg_type ("MathExpressionView::expression",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_EXPRESSION);

	object_class = GTK_OBJECT_CLASS (class);
	object_class->finalize = math_expression_view_finalize;
	object_class->set_arg = math_expression_view_set_arg;
	object_class->get_arg = math_expression_view_get_arg;

	widget_class = GTK_WIDGET_CLASS (class);
	widget_class->realize = math_expression_view_realize;
	widget_class->expose_event = math_expression_view_expose;
	widget_class->key_press_event = math_expression_view_key_press;
	widget_class->button_press_event = math_expression_view_button_press;
	widget_class->button_release_event =
		math_expression_view_button_release;
	widget_class->selection_get = math_expression_view_selection_get;
	widget_class->selection_received = 
		math_expression_view_selection_received;

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
	case ARG_EXPRESSION:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) == NULL ||
				  IS_MATH_EXPRESSION
				  (GTK_VALUE_POINTER (*arg)));

		if (math_expression_view->p->expression != NULL)
			gtk_object_unref
				(GTK_OBJECT
				 (math_expression_view->p->expression));

		math_expression_view->p->expression =
			MATH_EXPRESSION (GTK_VALUE_POINTER (*arg));

		if (math_expression_view->p->expression != NULL)
			gtk_object_ref
				(GTK_OBJECT
				 (math_expression_view->p->expression));
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
	case ARG_EXPRESSION:
		GTK_VALUE_POINTER (*arg) = math_expression_view->p->expression;
		break;

	default:
		g_warning ("Bad argument get");
		break;
	}
}

static void
math_expression_view_finalize (GtkObject *object) 
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (object != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (object));

	math_expression_view = MATH_EXPRESSION_VIEW (object);

	g_free (math_expression_view->p);

	GTK_OBJECT_CLASS (parent_class)->finalize (object);
}

GtkWidget *
math_expression_view_new (MathExpression *expr) 
{
	return gtk_widget_new (math_expression_view_get_type (),
			       "expression", expr,
			       NULL);
}

/**
 * math_expression_view_render:
 * @view: 
 * @area: The box to render; NULL to render the whole view
 * 
 * Render a portion of this view; render the full view if area is NULL
 **/

void
math_expression_view_render (MathExpressionView *view, GdkRectangle *area) 
{
	MathObject *toplevel;
	Layout *main_layout;
	GdkRectangle full_area;

	g_return_if_fail (view != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (view));

	if (!GTK_WIDGET_REALIZED (GTK_WIDGET (view))) return;

	full_area.x = full_area.y = 0;
	full_area.width = GTK_WIDGET (view)->allocation.width;
	full_area.height = GTK_WIDGET (view)->allocation.height;

	toplevel = math_expression_get_toplevel (view->p->expression);
	main_layout = math_object_get_layout (toplevel);
	layout_render (main_layout, toplevel, view->p->renderer, 
		       &full_area, area);
}

static void
math_expression_view_realize (GtkWidget *widget)
{
	MathExpressionView *math_expression_view;
	GdkWindowAttr attributes;
	gint attributes_mask;
	GdkColor color;
	GdkColormap *colormap;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	GTK_WIDGET_SET_FLAGS (widget, GTK_REALIZED);
	GTK_WIDGET_SET_FLAGS (widget, GTK_CAN_FOCUS);

	attributes.window_type = GDK_WINDOW_CHILD;
	attributes.x = widget->allocation.x;
	attributes.y = widget->allocation.y;
	attributes.width = MAX (widget->allocation.width, 400);
	attributes.height = MAX (widget->allocation.height, 300);
	attributes.wclass = GDK_INPUT_OUTPUT;
	attributes.visual = gtk_widget_get_visual (widget);
	attributes.colormap = gtk_widget_get_colormap (widget);
	attributes.event_mask = gtk_widget_get_events (widget) |
		GDK_ALL_EVENTS_MASK;

	attributes_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_VISUAL |
		GDK_WA_COLORMAP;

	widget->window = gdk_window_new (gtk_widget_get_parent_window (widget),
					 &attributes, attributes_mask);
	gdk_window_set_user_data (widget->window, widget);

	widget->style = gtk_style_attach (widget->style, widget->window);
	gtk_style_set_background (widget->style, widget->window,
				  GTK_STATE_NORMAL);

	colormap = gtk_widget_get_colormap (widget);
	gdk_color_white (colormap, &color);
	gdk_window_set_background (widget->window, &color);
}

static void
math_expression_view_expose (GtkWidget *widget, GdkEventExpose *event)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);
	math_expression_view_render (widget, &event->area);
}

static void
math_expression_view_key_press (GtkWidget *widget, GdkEventKey *event)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);


	g_warning("Key Pressed");
	/* FIXME */
}

static void
math_expression_view_button_press (GtkWidget *widget, GdkEventButton *button)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	/* FIXME */
}

static void
math_expression_view_button_release (GtkWidget *widget, GdkEventButton *button)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	/* FIXME */
}

static void
math_expression_view_selection_get (GtkWidget *widget, GtkSelectionData *data,
				    guint info, guint time)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	/* FIXME */
}

static void
math_expression_view_selection_received (GtkWidget *widget, 
					 GtkSelectionData *data, guint time)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	/* FIXME */
}
