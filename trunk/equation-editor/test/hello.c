#include <gtk/gtkmain.h>
#include <gtk/gtkwindow.h>
#include <gtk/gtktext.h>
#include <gtk/gtkbutton.h>
#include <gtk/gtkvbox.h>

void destroy( GtkWidget *widget, gpointer data ) {
   gtk_main_quit();
}

int main( int argc, char **argv ) {
   GtkWindow *window;
   GtkText *entry;
   GtkButton *button;
   GtkVBox *box;
   
   gtk_init( &argc, &argv );
   window = GTK_WINDOW( gtk_window_new( GTK_WINDOW_TOPLEVEL ) );
   gtk_window_set_title( window, "Hello Buttons!" );
   gtk_widget_set_usize( GTK_WIDGET( window ), 500, 500 );
   gtk_signal_connect( GTK_OBJECT( window ), "delete_event",
		       GTK_SIGNAL_FUNC( destroy ), NULL );
   
/*   gtk_container_border_width( GTK_CONTAINER( window ), 10 ); */
   
   box = GTK_VBOX( gtk_vbox_new( FALSE, 10 ) );
   gtk_container_add( GTK_CONTAINER( window ), GTK_WIDGET( box ) );
   
   entry = GTK_TEXT( gtk_text_new( NULL, NULL ) );
   gtk_text_set_editable( entry, TRUE );
   gtk_box_pack_start( GTK_BOX( box ), GTK_WIDGET( entry ), TRUE, TRUE, 0 );
   gtk_widget_set_usize( GTK_WIDGET( entry ), 480, 400 );
   gtk_widget_show( GTK_WIDGET( entry ) );
   
   button = GTK_BUTTON( gtk_button_new_with_label( "Quit" ) );
   gtk_signal_connect( GTK_OBJECT( button ), "clicked",
		       GTK_SIGNAL_FUNC( destroy ), NULL );
   gtk_box_pack_start( GTK_BOX( box ), GTK_WIDGET( button ), FALSE, FALSE, 0 );
   gtk_widget_set_usize( GTK_WIDGET( button ), 100, 20 );
   gtk_widget_show( GTK_WIDGET( button ) );
   
   gtk_widget_show( GTK_WIDGET( box ) );
   gtk_widget_show( GTK_WIDGET( window ) );
   
   gtk_main();
   return 0;
}
