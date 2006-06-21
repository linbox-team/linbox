import java.io.*;

/** A simple class to install one or both halves of the LinBox computation
  * server.
  */
class InstallServers {
	private static BufferedReader in =
		new BufferedReader( new InputStreamReader( System.in ) );
	
	/** This main method prompts the user for which parts he/she wants to
	  * install, and then goes from there.
	  */
	public static void main( String[] args ) {
	    boolean repeat;
	    do {
		System.out.println( "Linbox Server Installation options:" );
		System.out.println( "\ti: Install Computation Server Interface" );
		System.out.println( "\tc: Install Computing Server" );
		System.out.println( "\tb: Install both" );
		System.out.println( "\tq: Quit" );
		System.out.print( "[i/c/b/r/q]: " );
		String temp = null;
		try { temp = in.readLine(); }
		catch( IOException e ) { System.exit(1); }
		if( temp == null ) System.exit(1);
		boolean both = false;
		repeat = false;
		switch(temp.charAt(0)) {
			case 'q': System.exit(0);
			case 'b': both = true;
			case 'i': installInterface();
			          if( !both ) break;
			case 'c': installComputingServer();
			          break;
			default:  repeat = true;
		}
	    } while( repeat );
	}

	private static void installInterface() {
		System.out.println( "Don't know how to do this yet!" );
	}

	private static void installComputingServer() {
		System.out.println( "Don't know how to do this yet!" );
	}
}
