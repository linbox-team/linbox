/* -*- mode: c; style: linux -*- */

/* test.c
 * Copyright (C) 2000 Helix Code, Inc.
 *
 * Written by Bradford Hovinen (hovinen@helixcode.com)
 *            Anthony Asher, Rob Wehde, Matt Spilich
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

#include "math-expression.h"
#include "number.h"
#include "symbol.h"
#include "row-block.h"
#include "fraction-block.h"

static void about_cb (GtkWidget *widget);

static GnomeUIInfo file_menu[] = {
	GNOMEUIINFO_MENU_NEW_ITEM (N_("New formula"), 
				   N_("Create a new empty formula"),
				   NULL, NULL),
	GNOMEUIINFO_MENU_OPEN_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_SAVE_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_SAVE_AS_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_REVERT_ITEM (NULL, NULL),
	GNOMEUIINFO_SEPARATOR,
	GNOMEUIINFO_MENU_PRINT_ITEM (NULL, NULL),
	GNOMEUIINFO_SEPARATOR,
	GNOMEUIINFO_MENU_CLOSE_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_EXIT_ITEM (gtk_main_quit, NULL),
	GNOMEUIINFO_END
};

static GnomeUIInfo edit_menu[] = {
	GNOMEUIINFO_MENU_UNDO_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_REDO_ITEM (NULL, NULL),
	GNOMEUIINFO_SEPARATOR,
	GNOMEUIINFO_MENU_CUT_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_COPY_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_PASTE_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_SELECT_ALL_ITEM (NULL, NULL),
	GNOMEUIINFO_MENU_CLEAR_ITEM (NULL, NULL),
	GNOMEUIINFO_END
};

static GnomeUIInfo settings_menu[] = {
	GNOMEUIINFO_MENU_PROPERTIES_ITEM (NULL, NULL),
	GNOMEUIINFO_END
};

static GnomeUIInfo help_menu[] = {
	GNOMEUIINFO_MENU_ABOUT_ITEM (about_cb, NULL),
	GNOMEUIINFO_END
};

static GnomeUIInfo menu_bar[] = {
	GNOMEUIINFO_MENU_FILE_TREE (file_menu),
	GNOMEUIINFO_MENU_EDIT_TREE (edit_menu),
	GNOMEUIINFO_MENU_SETTINGS_TREE (settings_menu),
	GNOMEUIINFO_MENU_HELP_TREE (help_menu),
	GNOMEUIINFO_END
};

/**
 * setup_app_window:
 * @void: 
 * 
 * Create a basic shell window with suitable menus
 * 
 * Return value: 
 **/

static GtkWidget *
setup_app_window (MathExpression *expr) 
{
	GtkWidget *app;

	app = gnome_app_new ("test-program", "Test Program");
	gnome_app_create_menus (GNOME_APP (app), menu_bar);
	gnome_app_set_contents (app, math_expression_view_new (expr));

	gtk_window_set_default_size (GTK_WINDOW (app), 400, 300);

	gtk_signal_connect (GTK_OBJECT (app), "destroy", gtk_main_quit, NULL);

	gtk_widget_set_events( app, GDK_KEY_PRESS);

	gtk_widget_show_all (app);
	return app;
}
/**
*****************************************************
int
main (int argc, char **argv) 
{
	MathExpression *expr;
	Number *num1, *num2;
	Symbol *add_op;
	RowBlock *toplevel;

        bindtextdomain (PACKAGE, GNOMELOCALEDIR);
        textdomain (PACKAGE);

	gnome_init ("test-program", VERSION, argc, argv);

	num1 = NUMBER (number_new (1));
	add_op = SYMBOL (symbol_new ('+'));
	num2 = NUMBER (number_new (2));

	toplevel = ROW_BLOCK (row_block_new ());
	row_block_insert (toplevel, MATH_OBJECT (num1), NULL);
	row_block_insert (toplevel, MATH_OBJECT (add_op), NULL);
	row_block_insert (toplevel, MATH_OBJECT (num2), NULL);

	expr = math_expression_new (toplevel);

	setup_app_window (expr);

	gtk_main ();

	return 0;
}
**********************************************************
**/
int
main (int argc, char **argv) 
{
	MathExpression *expr;
	Number *num1, *num2;
	Symbol *add_op;
	FractionBlock *toplevel;
	RowBlock *N;
	RowBlock *D;

        bindtextdomain (PACKAGE, GNOMELOCALEDIR);
        textdomain (PACKAGE);

	gnome_init ("test-program", VERSION, argc, argv);

	num1 = NUMBER (number_new (1));
	/**	add_op = SYMBOL (symbol_new ('+')); **/
	num2 = NUMBER (number_new (2));

	N = ROW_BLOCK (row_block_new ());
	D = ROW_BLOCK (row_block_new ());

	row_block_insert (N, MATH_OBJECT (num1), NULL);
	row_block_insert (D, MATH_OBJECT (num2), NULL);

	toplevel = FRACTION_BLOCK (fraction_block_new(MATH_OBJECT(N), MATH_OBJECT(D)));
	/*fraction_block_set_numerator (toplevel, MATH_OBJECT (num1));
	row_block_insert (toplevel, MATH_OBJECT (add_op), NULL);
	fraction_block_set_denominator (toplevel, MATH_OBJECT (num2)); */

	expr = math_expression_new (toplevel);

	setup_app_window (expr);

	gtk_main ();

	return 0;
}

static void about_cb (GtkWidget *widget)
{
	static GtkWidget *about_dialog = NULL;
	static gchar *authors[] = {
		"Bradford Hovinen <hovinen@helixcode.com>",
		"Rob Wehde <robw@udel.edu>",
		"Matt Spilich <mspilich@udel.edu>",
		"Tony Asher <asher@udel.edu>",
		NULL
	};

	if (about_dialog == NULL) {
		about_dialog = gnome_about_new
			(_("GNOME Formula Editor"),
			 VERSION,
			 _("Copyright (C) 2000 Helix Code, Inc., "
			   "Rob Wehde, Matt Spilich, Tony Asher"),
			 authors,
			 _("GNOME mathematical formula editor"),
			 NULL);
	}

	gtk_widget_show_all (about_dialog);
}
