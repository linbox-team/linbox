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
#include "row-block.h"
#include "fraction-block.h"
#include "symbol.h"
#include "number.h"
#include "identifier.h"

enum {
	ARG_0,
	ARG_SAMPLE
};

struct _ControllerPrivate 
{
     /* Cursor current_pos; */
     MathObject *toplevel;
     MathObject *current_obj;
     int pos; 
     MathObject *previous_obj;
     MathObject *next_obj;
     MathObject *parent_obj;
	/* if current obj is a rowblock, pos gives what element in
	the rowblock is actually the current object */
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
	controller->p->current_obj = NULL;
	controller->p->next_obj = NULL;
	controller->p->previous_obj = NULL;
	controller->p->pos = 0;
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

void controller_initialize( Controller *controller, MathObject *toplevel)
{
	g_return_if_fail(IS_MATH_OBJECT(toplevel));
	g_return_if_fail(IS_CONTROLLER(controller));
	controller->p->current_obj = toplevel;
	controller->p->pos = 0;
}

void controller_insert (Controller *controller, GdkEventKey *event) 
{
	Symbol *symbol;
	gchar *this_keypressed;
	g_warning("In controller_insert");
	g_return_if_fail (controller->p->current_obj != NULL);
	g_return_if_fail (IS_MATH_OBJECT (controller->p->current_obj));
	this_keypressed = (event->string);
	g_return_if_fail (event != NULL);

/**	if (event->keyval == GDK_Tab) 
	{ controller_movenext( controller, controller->p->current_obj); }
	else 
**/
	if ( IS_ROW_BLOCK (controller->p->current_obj) ) {
		if(*(event->string)){
        	if ((event->state & GDK_CONTROL_MASK) == GDK_CONTROL_MASK) {
                if (event->keyval == GDK_f)
                        g_warning("Insert fraction call");  }
	else {
           symbol = SYMBOL( symbol_new(*(event->string)));
           row_block_insert_at(ROW_BLOCK(controller->p->current_obj),
           MATH_OBJECT(symbol), controller->p->pos); 
	   (controller->p->pos)++;  }	

	}
}}

#if 0

static void controller_movenext(Controller *controller, MathObject
*obj) {
int row_cnt, posit;
MathObject *thisobj;
posit = controller->p->pos - 1;
if (posit < 0) {posit = 0;}

if ( IS_ROW_BLOCK (obj) ) {
	row_cnt = row_block_get_length (obj);
	g_warning("pos: %d", posit); 
	g_warning("row_len: %d", row_block_get_length(obj));
	thisobj = row_block_get_object_at ( obj, posit );
	if ( IS_SYMBOL(thisobj) || IS_NUMBER(thisobj) ||
	     IS_IDENTIFIER(thisobj) ) {
		/* then we don't have to navigate down to the child
		   object - just advance the position */
		if ( controller->p->pos < row_cnt )
			(controller->p->pos)++;
		else { controller->p->pos = 0; }

	}
	if ( IS_FRACTION_BLOCK(thisobj) ) {
		/* then we need to set the current object to be the
		   fraction block */
	}
	else {
	   if (obj == controller->p->toplevel) 
		{ controller->p->pos = 0; }
	}
/   else {
		if parent_obj != NULL
	thisobj = row_block_get_object_at (obj, pos);
	if ( IS_FRACTION_BLOCK (thisobj) ) {
        	obj = fraction_block_get_numerator(thisobj);
        	controller->p->pos = 0;
        }
/
}

#endif

 

