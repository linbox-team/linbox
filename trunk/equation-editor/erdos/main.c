/* main.c - 
 *
 * Copyright (C) 1999 Chris Lahey.
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
 */



#include "config.h"

#include "main.h"
#include "color.h"
#include "erdos.h"
#include "erdos-expression.h"
#include "erdos-factory.h"
#include "item-erdos.h"
#include <gnome.h>
#include <gnome-xml/parser.h>

static void destroy_callback(GtkWidget *app, gpointer data)
{
  exit(0);
}

int main( int argc, char *argv[] )
{
  xmlDoc *doc;
  ErdosExpression *exp;
  GtkWidget *app;
  GtkWidget *canvas;
#if 0
  char **args;
  poptContext ctx;
#endif

  bindtextdomain (PACKAGE, GNOMELOCALEDIR);
  textdomain (PACKAGE);

  gnome_init( "Erdös", VERSION, argc, argv);
#if 0
  args = poptGetArgs(ctx);

  for(i = 0; args && args[i]; i++) {
      if ( args[i][0] == '/' )
	{
	  open_file( args[i] );
	}
      else
	{
	  gchar *current_dir = g_get_current_dir();	
	  gchar *fullname = append2( current_dir, TRUE, "/", FALSE );
	  fullname = append2( fullname, TRUE, args[i], FALSE );
	  open_file( fullname );
	  g_free( fullname );
	}
  }

  poptFreeContext(ctx);
#endif

  color_init();

  doc = xmlParseFile("test.mathml");
  exp = process_xml(doc->root);

  app = gnome_app_new("Erdös", NULL);

  canvas = gnome_canvas_new();
  gnome_canvas_item_new( gnome_canvas_root( GNOME_CANVAS( canvas) ),
			 item_erdos_get_type(),
			 "expression", exp,
			 NULL);
  gnome_canvas_set_scroll_region ( GNOME_CANVAS( canvas ),
				   0, 0,
				   100, 100 );

  gnome_app_set_contents( GNOME_APP( app ), canvas );


  /* Connect the signals */
  gtk_signal_connect( GTK_OBJECT( app ), "destroy",
		      GTK_SIGNAL_FUNC( destroy_callback ),
		      ( gpointer ) app );

  gtk_signal_connect( GTK_OBJECT( canvas ), "button_press_event",
		      GTK_SIGNAL_FUNC( about_callback ),
		      ( gpointer ) app );

  gtk_widget_show_all( app );

  gtk_main(); 

  /* Not reached. */
  return 0;
}

void about_callback( GtkWidget *widget, gpointer data )
{
  
  const gchar *authors[] =
  {
    "Christopher James Lahey <clahey@umich.edu>",
    NULL
  };

  GtkWidget *about =
    gnome_about_new ( _( "Erdös" ), VERSION,
		      _( "Copyright (C) 1999, Christopher James Lahey" ),
		      authors,
		      _( "This is named after the mathematician Paul Erdös" ),
		      NULL);
  gtk_widget_show (about);                                            
}
