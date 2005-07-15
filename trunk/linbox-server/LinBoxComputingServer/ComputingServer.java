package LinBoxComputingServer;

import java.net.Socket;
import java.net.ServerSocket;
import java.net.InetAddress;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.util.Hashtable;

/** An object of this class listens for connections from a (possibly remote)
  * ComputingServerInterface and, once connected, receives and performs
  * computations as requested.
  */
public class ComputingServer implements ComputationCodes, Locations {
	private static ComputingServer instance = null;
	private static boolean keepGoing = true;
	/** Tells the Computing server to bow out relatively gracefully and quit
	  */
	public static void die() {
		keepGoing = false;
		if( instance != null ) instance.killAllAndQuit();
	}
	
	/** This main method starts listening for connections, and once one is
	  * made, performs Computations as requested.
	  * @param args Any command-line arguments will be ignored.
	  */
	public static void main( String[] args ) {
		Socket s = null;
		ServerSocket ss = null;
		while( keepGoing ) {
		    try {
			if( ss == null ) ss = new ServerSocket(port);
		    	s = ss.accept();
			System.err.println( "Connected!" );
		    }
		    catch( IOException e ) {
		    	ss = null;
			s = null;
		    }
		    catch( Exception e ) {
			e.printStackTrace( System.err );
			System.exit(1);
		    }
		    if( s != null ) {
			instance = new ComputingServer(s);
			instance.go();
			instance.killAll();
			instance = null;
		    	s = null;
		    }
		    if( keepGoing ) {
		    	try { Thread.currentThread().sleep(TIMEOUT * 1000); }
		    	catch( InterruptedException e2 ) {}
		    }
		}
		System.exit(0);
	}

	private Socket sock;
	private DataInputStream in;
	private DataOutputStream out;
	private Hashtable computations = new Hashtable();

	/** Kill all running computations and exit the program. */
	public void killAllAndQuit() {
		//TODO: make sure this doesn't happen in the middle of reading
		//      a computation
		killAll();
		System.exit(0);
	}

	private ComputingServer( Socket s ) {
		sock = s;
		try {
			in = new DataInputStream( sock.getInputStream() );
			out = new DataOutputStream( sock.getOutputStream() );
		}
		catch( Exception e ) {
			e.printStackTrace( System.err );
			System.exit(2);
		}
	}

	private void go() {
	    boolean cont = true;
	    while( cont && sock.isConnected() ) {
		try {
		    int code = in.readInt();
		    if( code == NEW ) {
			System.err.println("New computation request received");
		    	int id = in.readInt();
			int type = in.readInt();
			String matURL = in.readUTF();
			long prime = in.readLong();
			Computation comp =
				new Computation(this, id, type, matURL, prime );
			computations.put( new Integer(id), comp );
			comp.start();
		    }
		    else if( code == KILL ) {
		    	int id = in.readInt();
			Object o = computations.get(new Integer(id));
			if( o != null ) ((Thread)o).interrupt();
		    }
		    else System.err.println( "Invalid code: " + code );
		}
		catch( Exception e ) {
			System.err.println( "Disconnected" );
			killAll();
			cont = false;
		}
	    }
	}

	private void killAll() {
		java.util.Enumeration comps = computations.elements();
		while( comps.hasMoreElements() )
			((Computation)comps.nextElement()).interrupt();
		if( sock != null ) {
			try { sock.close(); }
			catch( Exception e ) {}
		}
	}

	/** Send the results of a computation back to its requester.
	  * @param id The id of the completed computation
	  * @param res The result of the computation
	  */
	public synchronized void sendResults( int id, String res ) {
		computations.remove( new Integer(id) );
		if( !sock.isConnected() ) return;
		try {
			out.writeInt(DONE);
			out.writeInt(id);
			out.writeUTF(res);
			System.err.println( "Sent results for computation " + id );
		} catch( IOException e ) {
			System.err.println( "Warning: problems sending results for computation " + id );
		}
	}
}
