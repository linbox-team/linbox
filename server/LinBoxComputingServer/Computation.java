package LinBoxComputingServer;

import java.io.File;
import java.net.URL;
import java.io.InputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.IOException;

/** The ComputingServer's view of a computation.  Contains all information about
  * the computation, and has the ability to run the Computation.
  */
public class Computation extends Thread implements ComputationCodes {
	/** The size of the buffer to use when communicating with Computation
	  * processes.
	  */
	public static final int BUFSIZE = 1024;

	private ComputingServer cs;
	private String compName;
	private String matURL;
	private File matFile = null;
	private long prime;
	private int id;
	private Process p = null;

	/** Make (don't run) a new computation
	  * @param w An instance of a running ComputingServer.
	  * @param id The id of this computation
	  * @param code The code indicating what sort of computation this is
	  * @param url The url of the matrix data file
	  * @param p The prime number to compute mod, or -1 if computation is to
	  *          be performed over the integers.
	  */
	public Computation
	    (ComputingServer w, int id, int code, String url, long p ) {
		this.id = id;
		cs = w;
		setCompName( code );
		prime = p;
		matURL = url;
	}

	private void setCompName( int code ) {
		switch( code ) {
			case RANK: compName = "examples/rank"; break;
			case DET: compName = "examples/det"; break;
			default:
			    sendError( "Cannot perform that computation." );
			    break;
		}
	}

	/** Load the data from url(s) and start the computation in a separate
	  * process.
	  */
	public void run() {
		URL loc;
		try { loc = new URL( matURL ); }
		catch( IOException e ) {
			sendError( "Malformed URL.\n" );
			return;
		}
		try {
			InputStream in = loc.openStream();
			matFile = new File( cs.TEMP_DIR,
			                    Long.toString(id) + ".mat" );
			FileOutputStream out = new FileOutputStream(matFile);
			int l;
			byte[] buff = new byte[BUFSIZE];
			while( (l = in.read(buff)) > 0 )
				out.write( buff, 0, l );
			out.close();
			in.close();
		}
		catch( IOException e ) {
			sendError( "Error retrieving data from URL.\n" );
			return;
		}

		File progFile = new File( cs.LINBOX_FOLDER, compName );
		String[] cmd;
		if( prime == -1 ) 
			cmd = new String[] { progFile.getAbsolutePath(),
			                     matFile.getAbsolutePath() };
		else
			cmd = new String[] { progFile.getAbsolutePath(),
			                     matFile.getAbsolutePath(),
					     Long.toString( prime ) };
		try { p = (Runtime.getRuntime()).exec( cmd ); }
		catch( IOException e ) {
			sendError( "Could not run computation.\n" );
			return;
		}
		try { p.waitFor(); }
		catch( InterruptedException e ) {
			p.destroy();
			cleanUp();
			return;
		}
		StringBuffer sbuf = new StringBuffer();
		try {
			InputStreamReader in =
				new InputStreamReader( p.getErrorStream() );
			int l;
			char[] buff = new char[BUFSIZE];
			while( (l = in.read(buff,0,BUFSIZE)) > 0 )
				sbuf.append( buff, 0, l );
			in.close();
			in = new InputStreamReader( p.getInputStream() );
			while( (l = in.read(buff,0,BUFSIZE)) > 0 )
				sbuf.append( buff, 0, l );
			in.close();
		} catch( IOException e ) {
			sendError( "Error retrieving results.\n" );
			return;
		}
		cs.sendResults( id, sbuf.toString() );
		cleanUp();
	}

	private void sendError( String e ) {
		cs.sendResults( id, e );
		cleanUp();
	}

	private void cleanUp() {
		if( matFile != null ) matFile.delete();
	}
}
