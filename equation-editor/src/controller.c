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

#include <ctype.h>

#include "controller.h"
#include "cursor.h"
#include "math-expression-view.h"
#include "row-block.h"
#include "fraction-block.h"
#include "math-atom.h"

typedef enum _ControllerState ControllerState;

enum {
	ARG_0,
	ARG_EXPR,
	ARG_VIEW
};

enum _ControllerState {
	STATE_NORMAL, STATE_TEXT_ENTRY
};

struct _ControllerPrivate 
{
	MathExpression *expr;
	MathExpressionView *view;
	Cursor *cursor;
	ControllerState state;
};

static GtkObjectClass *parent_class;

static void     controller_init        (Controller *controller);
static void     controller_class_init  (ControllerClass *class);

static void     controller_set_arg     (GtkObject *object, 
					GtkArg *arg, 
					guint arg_id);
static void     controller_get_arg     (GtkObject *object, 
					GtkArg *arg, 
					guint arg_id);

static void     controller_finalize    (GtkObject *object);

static MathObject *make_fraction       (Controller *controller);
static void        key_press           (Controller *controller,
					GdkEventKey *event);
static gboolean    do_insert_character (Controller *controller,
					MathObject *math_object,
				        gint keycode);

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
	controller->p->state = STATE_NORMAL;
}

static void
controller_class_init (ControllerClass *class) 
{
	GtkObjectClass *object_class;

	gtk_object_add_arg_type ("Controller::expr",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_EXPR);
	gtk_object_add_arg_type ("Controller::view",
				 GTK_TYPE_POINTER,
				 GTK_ARG_READWRITE,
				 ARG_VIEW);

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
	case ARG_EXPR:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) != NULL);
		g_return_if_fail (IS_MATH_EXPRESSION 
				  (GTK_VALUE_POINTER (*arg)));

		if (controller->p->expr != NULL) {
			gtk_object_unref (GTK_OBJECT (controller->p->expr));
			gtk_object_unref (GTK_OBJECT (controller->p->cursor));
		}

		controller->p->expr = GTK_VALUE_POINTER (*arg);

		if (controller->p->expr != NULL) {
			gtk_object_ref (GTK_OBJECT (controller->p->expr));
			controller->p->cursor = 
				CURSOR (cursor_new (controller->p->expr));

			if (controller->p->view != NULL) {
				gtk_signal_connect_object
					(GTK_OBJECT (controller->p->cursor),
					 "moved",
					 GTK_SIGNAL_FUNC (math_expression_view_render_by_object),
					 GTK_OBJECT (controller->p->view));
			}
		}

		break;

	case ARG_VIEW:
		g_return_if_fail (GTK_VALUE_POINTER (*arg) != NULL);
		g_return_if_fail (IS_MATH_EXPRESSION_VIEW
				  (GTK_VALUE_POINTER (*arg)));

		if (controller->p->view != NULL)
			gtk_object_unref (GTK_OBJECT (controller->p->view));

		controller->p->view = GTK_VALUE_POINTER (*arg);

		if (controller->p->view != NULL) {
			gtk_object_ref (GTK_OBJECT (controller->p->view));

			/* Connect view signals for the key press event */
			gtk_signal_connect_object
				(GTK_OBJECT (controller->p->view),
				 "key-press-event",
				 GTK_SIGNAL_FUNC (key_press),
				 GTK_OBJECT (controller));

			if (controller->p->cursor != NULL) {
				gtk_signal_connect_object
					(GTK_OBJECT (controller->p->cursor),
					 "moved",
					 GTK_SIGNAL_FUNC (math_expression_view_render_by_object),
					 GTK_OBJECT (controller->p->view));
			}

			gtk_object_set (GTK_OBJECT (controller->p->view),
					"controller", controller,
					NULL);
		}

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
	case ARG_EXPR:
		GTK_VALUE_POINTER (*arg) = controller->p->expr;
		break;

	case ARG_VIEW:
		GTK_VALUE_POINTER (*arg) = controller->p->view;
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
controller_new (MathExpression *expr, MathExpressionView *view) 
{
	return gtk_object_new (controller_get_type (),
			       "expr", expr,
			       "view", view,
			       NULL);
}

/**
 * controller_get_cursor:
 * @controller: 
 * 
 * Get the cursor for this controller
 * 
 * Return value: The cursor; should be unrefed
 **/

Cursor *
controller_get_cursor (Controller *controller)
{
	g_return_val_if_fail (controller != NULL, NULL);
	g_return_val_if_fail (IS_CONTROLLER (controller), NULL);

	if (controller->p->cursor != NULL)
		gtk_object_ref (GTK_OBJECT (controller->p->cursor));
	return controller->p->cursor;
}

/**
 * make_fraction:
 * @controller: 
 * 
 * Create a fraction with empty row blocks for the numerator and denominator
 * and return it
 **/

static MathObject *
make_fraction (Controller *controller) 
{
	MathObject *new_object = NULL, *num, *den;

	num = MATH_OBJECT (row_block_new ());
	gtk_signal_connect_object
		(GTK_OBJECT (num), "changed",
		 GTK_SIGNAL_FUNC (math_expression_view_render_by_object),
		 GTK_OBJECT (controller->p->view));
	den = MATH_OBJECT (row_block_new ());
	gtk_signal_connect_object
		(GTK_OBJECT (den), "changed",
		 GTK_SIGNAL_FUNC (math_expression_view_render_by_object),
		 GTK_OBJECT (controller->p->view));
	new_object = MATH_OBJECT (fraction_block_new (num, den));

	return new_object;
}

/**
 * key_press:
 * @controller: object
 * @event: The key press event
 * 
 * Key press event handler; interprets the keypress and dispatches the
 * appropriate event
 **/

static void
key_press (Controller *controller, GdkEventKey *event)
{
	MathObject *new_object = NULL, *current;

	g_return_if_fail (controller != NULL);
	g_return_if_fail (IS_CONTROLLER (controller));

	current = cursor_get_current_object (controller->p->cursor);

	if (event->state & GDK_CONTROL_MASK) {
		switch (event->keyval) {
		case GDK_r: /* Insert fraction */
			new_object = make_fraction (controller);
			break;

		case GDK_a:
			cursor_move_to_beginning (controller->p->cursor);
			break;

		case GDK_e:
			cursor_move_to_end (controller->p->cursor);
			break;

		case GDK_f:
			cursor_move_right (controller->p->cursor);
			break;

		case GDK_b:
			cursor_move_left (controller->p->cursor);
			break;

		default: /* Invalid keypress */
			break;
		}
	}
	else if (event->state & GDK_MOD1_MASK) {
		switch (event->keyval) {
		case GDK_t:  /* Enter "text edit" mode */
			new_object = MATH_OBJECT 
				(math_atom_new (MATH_ATOM_DIVSTRING));
			controller->p->state = STATE_TEXT_ENTRY;

		default: /* Invalid keypress */
			break;
		}
	} else {
		if (do_insert_character (controller, current, event->keyval)) {
			math_atom_insert
				(MATH_ATOM (current),
				 cursor_get_insertion_point
				 (controller->p->cursor),
				 event->keyval);
			cursor_move_right (controller->p->cursor);
		} else {
			switch (event->keyval) {
			case GDK_Left:
				cursor_move_left (controller->p->cursor);
				break;
			case GDK_Right:
				cursor_move_right (controller->p->cursor);
				break;
			case GDK_BackSpace:
				cursor_move_left (controller->p->cursor);
/*  				eat_next_object (controller); */
				break;
			default:
				if (event->keyval < 128) {
					new_object = MATH_OBJECT
						(math_atom_new
						 (MATH_ATOM_DIVSTRING));
					math_atom_append
						(MATH_ATOM (new_object),
						 event->keyval);
					break;
				}
			}
		}
	}

	if (new_object != NULL) {
		g_assert (IS_ROW_BLOCK (current));

		gtk_signal_connect_object
			(GTK_OBJECT (new_object), "changed",
			 GTK_SIGNAL_FUNC
			 (math_expression_view_render_by_object),
			 GTK_OBJECT (controller->p->view));

		row_block_insert_at (ROW_BLOCK (current), new_object,
				     cursor_get_insertion_point
				     (controller->p->cursor));
		cursor_move_right (controller->p->cursor);
	}

	gtk_object_unref (GTK_OBJECT (current));
}

/**
 * do_insert_character:
 * @math_object: 
 * @keycode: The keycode as passed by GDK
 * 
 * Check if the key indicated should be appended onto the math object
 * 
 * Return value: TRUE if it should be appended, FALSE otherwise
 **/

static gboolean
do_insert_character (Controller *controller, MathObject *math_object,
		     gint keycode) 
{
	if (IS_MATH_ATOM (math_object) &&
	    (isdigit (keycode) ||
	     (isalpha (keycode) &&
	      controller->p->state == STATE_TEXT_ENTRY)))
	    return TRUE;

	controller->p->state = STATE_NORMAL;
	return FALSE;
}
