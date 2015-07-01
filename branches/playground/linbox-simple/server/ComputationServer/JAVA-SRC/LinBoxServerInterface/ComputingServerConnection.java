package LinBoxServerInterface;

import java.net.Socket;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.io.IOException;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.Enumeration;

/** This is the LinBoxServerInterface's connection to a Computing Server.
  * There should be at most one instance of this class at any given time.  To
  * get an instance that is up and running, use getInstance().  If the
  * connection is dropped or reset, it will keep trying to connect automatically
  */
public class ComputingServerConnection extends Thread
                             implements ComputationCodes, Locations {
	static private ComputingServerConnection instance = null;
	/** Get an instance of this class.  If none has been started yet, a new
	  * connection will be made.  Otherwise, the existing instance will be
	  * returned.
	  * @return An running ComputingServerConnection.
	  */
	public static ComputingServerConnection getInstance() {
		if( instance == null ) {
			instance = new ComputingServerConnection();
			instance.start();
		}
		return instance;
	}
	
	/** The amount of time, in seconds, to wait before attempting a
	  * reconnect when the connection fails.
	  */
	public static final int TIMEOUT = 1;
	private static InetAddress hostAddr;
	private Socket sock = null;

	protected ComputingServerConnection() 
	{
		try { hostAddr = InetAddress.getByName( hostName ); }
		catch( UnknownHostException e )
		{ throw new RuntimeException( e ); }
	}

	/** Attempt to make and maintain a connection to the remote
	  * (or local) ComputingServer.
	  */
	public void run() {		
	    boolean fail = true;
	    while( fail ) {
	    	try { sock = new Socket( hostAddr, hostPort ); }
		catch( IOException e ) {}
		if( sock != null ) {
			resendAll();
			listenForResponses();
			try { sock.close(); }
			catch( Exception e ) {}
			sock = null;
		}
		try { Thread.currentThread().sleep( 1000 * TIMEOUT ); }
		catch( InterruptedException e ) {}
	    }
	}

	private void listenForResponses() {
	    DataInputStream din = null;
	    try { din = new DataInputStream( sock.getInputStream() ); }
	    catch( IOException e ) { return; }
	    boolean keepGoing = true;
	    while( keepGoing && sock.isConnected() ) {
	      try {
	    	int code = din.readInt();
		if( code == DONE ) {
			int id = din.readInt();
			String resp = din.readUTF();
			Computation.saveResponse( id, resp );
		}
		else keepGoing = false;
	      } catch( IOException e ) { keepGoing = false; }
	    }
	}

	private synchronized void resendAll() {
		Enumeration comps = Computation.getAll();
		Computation thisOne;
		while( comps.hasMoreElements() ) {
			thisOne = (Computation)comps.nextElement();
			if( thisOne.getResult() == null ) {
				try { sendNew( thisOne ); }
				catch( IOException e ) {}
			}
		}
	}

	/** Send a new computation to the ComputingServer.
	  * @param comp A new Computation to be performed.
	  */
	public synchronized void sendNew( Computation comp ) 
		throws IOException
	{
	    try {
		while( sock == null ) {
		    try { Thread.currentThread().sleep( 1000 * TIMEOUT ); }
		    catch( InterruptedException e ) {}
		}
		DataOutputStream dout =
			new DataOutputStream( sock.getOutputStream() );
		dout.writeInt( NEW );
		comp.write( dout );
		comp.setSubmittedOn( new java.util.Date() );
		dout.flush();
	    } catch( Exception e ) { throw new IOException(e.getMessage()); }
	}

	/** Send a message to the ComputingServer to kill the specified
	  * Computation, if it is still running.
	  * @param id The id of the Computation to kill.
	  */
	public synchronized void kill( int id ) {
		while( sock == null ) {
		    try { Thread.currentThread().sleep( 1000 * TIMEOUT ); }
		    catch( InterruptedException e ) {}
		}
		try {
			DataOutputStream dout =
				new DataOutputStream( sock.getOutputStream() );
			dout.writeInt( KILL );
			dout.writeInt( id );
			dout.flush();
		} catch( IOException e ) {}
	}
}
