/* -*- mode: c; style: linux -*- */

/* bonobo-embeddable.c
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

#include <gnome.h>
#include <bonobo.h>
#include <liboaf/liboaf.h>

static void
view_activate_cb (BonoboView *view, gboolean activate, void *closure) 
{
	bonobo_view_activate_notify (view, activate);
}

static BonoboView *
bonobo_text_view_factory (BonoboEmbeddable *embeddable, 
			  const Bonobo_ViewFrame view_frame,
			  void *closure) 
{
	GtkWidget *textedit;
	BonoboView *view;

	/* Create text edit */
	textedit = gtk_text_new (NULL, NULL);
	gtk_text_insert (GTK_TEXT (textedit), NULL, NULL, NULL, 
			 "This is some text", strlen ("This is some text"));
	gtk_text_set_editable (GTK_TEXT (textedit), TRUE);
	gtk_widget_show (textedit);

	view = bonobo_view_new (textedit);
	gtk_signal_connect (GTK_OBJECT (view), "activate",
			    GTK_SIGNAL_FUNC (view_activate_cb), NULL);

	return view;
}

static BonoboObject *
bonobo_text_factory (BonoboGenericFactory *factory, gpointer data) 
{
	BonoboEmbeddable  *embeddable;

	embeddable = bonobo_embeddable_new (bonobo_text_view_factory, NULL);

	return BONOBO_OBJECT (embeddable);
}

static void
bonobo_text_factory_init (void) 
{
	static BonoboGenericFactory *bonobo_textedit_factory_obj = NULL;

	if (bonobo_textedit_factory_obj != NULL)
		return;

	bonobo_textedit_factory_obj =
		bonobo_generic_factory_new
		("OAFIID:embeddable-factory:bonobo-textedit:31ab377b-0743-4777-a6ef-d9e58f556772",
		 bonobo_text_factory, NULL);

	if (bonobo_textedit_factory_obj == NULL)
		g_error ("Could not register factory");
}

int 
main (int argc, char **argv) 
{
	CORBA_Environment ev;
	CORBA_ORB orb;

	CORBA_exception_init (&ev);

	gnome_init_with_popt_table ("text-factory", "0.0",
				    argc, argv, oaf_popt_options,
				    0, NULL);

	orb = oaf_init (argc, argv);

	if (bonobo_init (orb, NULL, NULL) == FALSE)
		g_error (_("Could not initialize Bonobo"));

	bonobo_text_factory_init ();
	bonobo_main ();

	return 0;
}
