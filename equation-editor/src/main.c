/* -*- mode: c; style: linux -*- */

/* main.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen (hovinen@helixcode.com)
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
#  include <config.h>
#endif

#include <gnome.h>
#include <bonobo.h>
#include <liboaf/liboaf.h>

#include "math-expression.h"
#include "math-expression-view.h"
#include "row-block.h"

#include "number.h"
#include "symbol.h"

static void
view_activate_cb (BonoboView *view, gboolean activate, void *closure) 
{
	bonobo_view_activate_notify (view, activate);
}

static BonoboView *
equation_view_factory (BonoboEmbeddable *embeddable,
		       const Bonobo_ViewFrame view_frame,
		       void *closure) 
{
	MathExpression *expr;
	MathExpressionView *view;
	BonoboView *bonobo_view;

	expr = MATH_EXPRESSION (closure);
	view = MATH_EXPRESSION_VIEW (math_expression_view_new (expr));
	bonobo_view = bonobo_view_new (GTK_WIDGET (view));

        gtk_signal_connect (GTK_OBJECT (bonobo_view), "activate",
                            GTK_SIGNAL_FUNC (view_activate_cb), NULL);

	return bonobo_view;
}

static BonoboObject *
equation_factory (BonoboGenericFactory *factory, gpointer data) 
{
	BonoboEmbeddable *embeddable;
	MathExpression *expr;
	RowBlock *toplevel;

	Number *num1, *num2;
	Symbol *add_op;

	num1 = NUMBER (number_new (1));
	add_op = SYMBOL (symbol_new ('+'));
	num2 = NUMBER (number_new (2));

	toplevel = ROW_BLOCK (row_block_new ());
	row_block_insert (toplevel, MATH_OBJECT (num1), NULL);
	row_block_insert (toplevel, MATH_OBJECT (add_op), NULL);
	row_block_insert (toplevel, MATH_OBJECT (num2), NULL);

	expr = MATH_EXPRESSION (math_expression_new (MATH_OBJECT (toplevel)));

	embeddable = bonobo_embeddable_new (equation_view_factory, expr);

	return BONOBO_OBJECT (embeddable);
}

static void
bonobo_equation_factory_init (void) 
{
	static BonoboGenericFactory *equation_factory_obj = NULL;

	if (equation_factory_obj != NULL)
		return;

	equation_factory_obj =
		bonobo_generic_factory_new
		("OAFIID:embeddable-factory:bonobo-equation:cc87cc19-f9b1-449d-862d-3744186dc937",
		 equation_factory, NULL);

	if (equation_factory_obj == NULL)
		g_error ("Could not register factory");
}

int
main (int argc, char **argv) 
{
	CORBA_Environment ev;
	CORBA_ORB orb;

        bindtextdomain (PACKAGE, GNOMELOCALEDIR);
        textdomain (PACKAGE);

	CORBA_exception_init (&ev);

	gnome_init_with_popt_table ("equation-editor", VERSION, argc, argv,
				    oaf_popt_options, 0, NULL);

	orb = oaf_init (argc, argv);

	if (bonobo_init (orb, NULL, NULL) == FALSE)
		g_error ("Could not initialize Bonobo");

	bonobo_equation_factory_init ();
	bonobo_main ();

	return 0;
}
