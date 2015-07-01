#include <gnome.h>
#include <libgnorba/gnorba.h>

#include "mymodule.h"

int
main (int argc, char **argv) 
{
	CORBA_Environment ev;
	CORBA_ORB orb;
	MyModule_MyInterface object;
	char buffer[2048], *ior;

	CORBA_exception_init (&ev);
	gnome_CORBA_init ("mymodule-client", "1.0", &argc, argv, 0, &ev);
	orb = gnome_CORBA_ORB ();

	fgets (buffer, 2048, stdin);
	buffer[strlen (buffer) - 1] = '\0';
	ior = g_strdup (buffer);

	object = CORBA_ORB_string_to_object (orb, ior, &ev);

	if (ev._major != CORBA_NO_EXCEPTION) {
		g_error ("Exception during query: %s", 
			 CORBA_exception_id (&ev));
	}

	g_message ("Result: %s", MyModule_MyInterface_randomize (object, "Hello, world", &ev));
}
