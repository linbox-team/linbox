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
#include "symbol.h"
#include "number.h"
#include "row-block.h"

enum {
	ARG_0,
	ARG_EXPRESSION
};

struct _MathExpressionViewPrivate 
{
	MathExpression *expression;
	Renderer *renderer;
	/* Controller *controlObj; */
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
static void math_expression_view_enter       (GtkWidget *widget,
					      GdkEventCrossing *event);

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
	widget_class->enter_notify_event = math_expression_view_enter;

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
		GDK_KEY_PRESS_MASK | GDK_EXPOSURE_MASK | GDK_ENTER_NOTIFY_MASK;

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
	MathObject *toplevel;
	MathObject *this_key;
	Symbol *symbol;
	Number *number;
	MathExpression *expr;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	math_expression_view = MATH_EXPRESSION_VIEW (widget);

	expr = math_expression_view->p->expression;

	toplevel = math_expression_get_toplevel(expr);

	if ((event->state & GDK_SHIFT_MASK) == GDK_SHIFT_MASK) {
		if (event->keyval == GDK_L)
			g_warning("Shift L entered");
	}	

	
	switch(*(event->string))
	{
	case '+':g_warning("+ entered");
	  	symbol = SYMBOL( symbol_new('+'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '-':g_warning("- entered");
	  	symbol = SYMBOL( symbol_new('-'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '*':g_warning("* entered");
	  	symbol = SYMBOL( symbol_new('*'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '/':g_warning("/ entered");
	  	symbol = SYMBOL( symbol_new('/'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '^':g_warning("^ entered");
	  	symbol = SYMBOL( symbol_new('^'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '(':g_warning("( entered");
	  	symbol = SYMBOL( symbol_new('('));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case ')':g_warning(") entered");
	  	symbol = SYMBOL( symbol_new(')'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '<':g_warning("< entered");
	  	symbol = SYMBOL( symbol_new('<'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '>':g_warning("> entered");
	  	symbol = SYMBOL( symbol_new('>'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '&':g_warning("& entered");
	  	symbol = SYMBOL( symbol_new('&'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '|':g_warning("| entered");
	  	symbol = SYMBOL( symbol_new('|'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '=':g_warning("= entered");
	  	symbol = SYMBOL( symbol_new('='));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'a':g_warning("a entered");
	  	symbol = SYMBOL( symbol_new('a'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'b':g_warning("b entered");
	  	symbol = SYMBOL( symbol_new('b'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'c':g_warning("c entered");
	  	symbol = SYMBOL( symbol_new('c'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'd':g_warning("d entered");
	  	symbol = SYMBOL( symbol_new('d'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'e':g_warning("e entered");
	  	symbol = SYMBOL( symbol_new('e'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'f':g_warning("f entered");
	  	symbol = SYMBOL( symbol_new('f'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case 'g':g_warning("g entered");
	  	symbol = SYMBOL( symbol_new('g'));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(symbol),NULL); break;
	case '0':g_warning("0 entered");
	  	number = NUMBER( number_new(0));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '1':g_warning("1 entered");
	  	number = NUMBER( number_new(1));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '2':g_warning("2 entered");
	  	number = NUMBER( number_new(2));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '3':g_warning("3 entered");
	  	number = NUMBER( number_new(3));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '4':g_warning("4 entered");
	  	number = NUMBER( number_new(4));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '5':g_warning("5 entered");
	  	number = NUMBER( number_new(5));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '6':g_warning("6 entered");
	  	number = NUMBER( number_new(6));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '7':g_warning("7 entered");
	  	number = NUMBER( number_new(7));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '8':g_warning("8 entered");
	  	number = NUMBER( number_new(8));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	case '9':g_warning("9 entered");
	  	number = NUMBER( number_new(9));
	  	row_block_insert(ROW_BLOCK(toplevel),
		MATH_OBJECT(number),NULL); break;
	default: break;
	}

/**	if ( *(event->string) == '0') {
		g_warning("0 entered");
		num = NUMBER( number_new(0));
		row_block_insert(ROW_BLOCK(toplevel), MATH_OBJECT(num), NULL); 
	}

	if ( *(event->string) == 'a')
		g_warning("a entered");

	if ( *(event->string) == 'b')
		g_warning("b entered");

**/	/* FIXME */
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

static void 
math_expression_view_enter (GtkWidget *widget, GdkEventCrossing *event)
{
	MathExpressionView *math_expression_view;

	g_return_if_fail (widget != NULL);
	g_return_if_fail (IS_MATH_EXPRESSION_VIEW (widget));

	gtk_widget_grab_focus (widget);
}
