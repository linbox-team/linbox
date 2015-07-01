#include <gnome.h>
#include <libgnorba/gnorba.h>

#include "mymodule.h"

/*** App-specific servant structures ***/
typedef struct {
   POA_MyModule_MyInterface servant;
   PortableServer_POA poa;

} impl_POA_MyModule_MyInterface;

/*** Implementation stub prototypes ***/
static void impl_MyModule_MyInterface__destroy(impl_POA_MyModule_MyInterface * servant,
					       CORBA_Environment * ev);
static CORBA_char *
 impl_MyModule_MyInterface_randomize(impl_POA_MyModule_MyInterface * servant,
				     CORBA_char * str,
				     CORBA_Environment * ev);

/*** epv structures ***/
static PortableServer_ServantBase__epv impl_MyModule_MyInterface_base_epv =
{
   NULL,			/* _private data */
   NULL,			/* finalize routine */
   NULL,			/* default_POA routine */
};
static POA_MyModule_MyInterface__epv impl_MyModule_MyInterface_epv =
{
   NULL,			/* _private */
   (gpointer) & impl_MyModule_MyInterface_randomize,

};

/*** vepv structures ***/
static POA_MyModule_MyInterface__vepv impl_MyModule_MyInterface_vepv =
{
   &impl_MyModule_MyInterface_base_epv,
   &impl_MyModule_MyInterface_epv,
};

/*** Stub implementations ***/
static MyModule_MyInterface 
impl_MyModule_MyInterface__create(PortableServer_POA poa, CORBA_Environment * ev)
{
   MyModule_MyInterface retval;
   impl_POA_MyModule_MyInterface *newservant;
   PortableServer_ObjectId *objid;

   newservant = g_new0(impl_POA_MyModule_MyInterface, 1);
   newservant->servant.vepv = &impl_MyModule_MyInterface_vepv;
   newservant->poa = poa;
   POA_MyModule_MyInterface__init((PortableServer_Servant) newservant, ev);
   objid = PortableServer_POA_activate_object(poa, newservant, ev);
   CORBA_free(objid);
   retval = PortableServer_POA_servant_to_reference(poa, newservant, ev);

   return retval;
}

static void
impl_MyModule_MyInterface__destroy(impl_POA_MyModule_MyInterface * servant, CORBA_Environment * ev)
{
   PortableServer_ObjectId *objid;

   objid = PortableServer_POA_servant_to_id(servant->poa, servant, ev);
   PortableServer_POA_deactivate_object(servant->poa, objid, ev);
   CORBA_free(objid);

   POA_MyModule_MyInterface__fini((PortableServer_Servant) servant, ev);
   g_free(servant);
}

static CORBA_char *
impl_MyModule_MyInterface_randomize(impl_POA_MyModule_MyInterface * servant,
				    CORBA_char * str,
				    CORBA_Environment * ev)
{
   CORBA_char *retval;

   retval = ORBit_alloc (strlen (str) + 1, NULL, NULL);
   strcpy (retval, str);
   strfry (retval);

   g_message ("randomize called");

   return retval;
}

int
main (int argc, char **argv) 
{
	CORBA_Environment ev;
	CORBA_ORB orb;
	PortableServer_POA poa;
	PortableServer_POA root_poa;
	PortableServer_POAManager root_poa_manager;
	CORBA_PolicyList policy_list = { 0, 0, NULL };
	MyModule_MyInterface mymodule;

	CORBA_exception_init (&ev);
	orb = gnome_CORBA_init ("mymodule", "1.0", &argc, argv,
				GNORBA_INIT_SERVER_FUNC, &ev);
	orb = gnome_CORBA_ORB ();

	root_poa = CORBA_ORB_resolve_initial_references (orb, "RootPOA", &ev);
	root_poa_manager = PortableServer_POA__get_the_POAManager (root_poa,
								   &ev);
	PortableServer_POAManager_activate (root_poa_manager, &ev);
	PortableServer_POA_create_POA (root_poa, "my_own_poa",
				       root_poa_manager, &policy_list, &ev);
	poa = PortableServer_POA_find_POA (root_poa, "my_own_poa", 
					   CORBA_TRUE, &ev);

	mymodule = impl_MyModule_MyInterface__create (poa, &ev);

	g_message ("IOR is %s", CORBA_ORB_object_to_string (orb, mymodule, &ev));

	gtk_main ();

	CORBA_free (root_poa);
	CORBA_exception_free (&ev);

	return 0;
}
