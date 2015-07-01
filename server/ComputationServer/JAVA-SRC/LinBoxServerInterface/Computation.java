package LinBoxServerInterface;

import java.util.*;
import java.io.*;
import java.beans.XMLEncoder;
import java.beans.XMLDecoder;

/** The main server-side representation of a computation.
 *  Holds all pertinent information about computations, as well as a number of
 *  useful static methods pertaining to all computations in memory.
 */
public class Computation implements ComputationCodes, Locations {
	
	static private Hashtable all = new Hashtable();
	
	/** Get a Computation by its id.
	  * @param id The id number of the computation to find and return.
	  * @return The computation with the given id, or null if none is found.
	  */
	static public Computation getByID( int id ) {
		return (Computation)all.get(new Integer(id));
	}

	/** Get an enumeration of all Computations in memory.
	  * @return An enumeration of all Computations in memory.
	  */
	static public Enumeration getAll() {
		return all.elements();
	}

	static private Random generator = new Random();

	/** Destroy this Computation and all record of it.
	  * The Computation will not be destroyed if there is any thread waiting
	  * for it, i.e. a webpage waiting for the computation to finish.
	  */
	public synchronized void tryDestroy() {
		if( waiting.size() == 0 ) {
			if( result == null )
				ComputingServerConnection.getInstance()
					.kill( id );
			all.remove(new Integer(id));
			if( matrixFile != null ) {
			    try {
			    	matrixFile.delete();
			    } catch( Exception e ) {}
			}
			if( archive != null ) {
			    try {
			    	archive.delete();
			    } catch( Exception e ) {}
			}
		}
	}

	/** Store the server's response for a given computation.
	  * Does nothing if no computation with the given id can be found.
	  * @param id The id of the Computation that got a response.
	  * @param resp The server's response
	  */
	public static synchronized void saveResponse( int id, String resp ) {
		Computation comp = getByID( id );
		if( comp != null ) {
			comp.setResult( resp );
			comp.saveArchive();
			comp.interruptWaiting();
		}
	}

	/** Loads any computations that have been saved to the disk.
	  * Useful e.g. when restarting from an unexpected reboot.
	  */
	public static synchronized void loadArchives() {
		File[] list = computationDir.listFiles
			( new ExtensionFilter(".xml") );
		Computation temp;
		XMLDecoder xin;
		for( int i = 0; i < list.length; ++i ) {
			xin = null;
			try {
				xin = new XMLDecoder(
					new FileInputStream( list[i] )
					            );
				temp = (Computation) xin.readObject();
				temp.archive = list[i];
			}
			catch( Exception e ) { }
			finally { if( xin != null ) xin.close(); }
		}
	}

	static { 
		loadArchives();
		ComputationSweeper.scheduleSweeps();
	}
	
	private int type = 0; // Should be one of ComputationCodes
	/** Get the type code of this computation.
	  * @return The type code of this computation.
	  */
	public int getType() { return type; }
	/** Set the type code of this computation.
	  * @param t The type code of this computation.
	  */
	public void setType( int t ) { type = t; }

	private Date createdOn = new Date();
	/** Get the date when this Computation was created (by the user).
	  * @return The creation date
	  */
	public Date getCreatedOn() { return createdOn; }
	/** Set the date when this Computation was first created.
	  * Should only be used by XMLDecoder.
	  * @param c The new created date.
	  */
	public void setCreatedOn( Date c ) { createdOn = c; }

	private Date replyReceivedOn = null;
	/** Get the date when a reply was received from the server.
	  * @return The date a reply was received.
	  */
	public Date getReplyReceivedOn() { return replyReceivedOn; }
	/** Set the date when a reply was received from the server.
	  * @param r The date a reply was received.
	  */
	public void setReplyReceivedOn( Date r ) { replyReceivedOn = r; }

	private Date submittedOn = null;
	/** Get the date when this computation was submitted to the server.
	  * @return The submitted date
	  */
	public Date getSubmittedOn() { return submittedOn; }
	/** Set the date this computation was submitted to the server.
	  * @param s The submitted date
	  */
	public void setSubmittedOn( Date s ) { submittedOn = s; }

	private Date lastAccessed = new Date();
	/** Get the date when this computation was last accessed by the user.
	  * @return The date of last user access
	  */
	public Date getLastAccessed() { return lastAccessed; }
	/** Set the date when this computation was last accessed by the user.
	  * @param d The date of most recent user access.
	  */
	public void setLastAccessed( Date d ) { lastAccessed = d; }
	/** Set the last user access date to the current date. */
	public void updateAccessed() { lastAccessed = new Date(); }
	
	private Vector waiting = new Vector();
	/** Get an Enumeration of the threads waiting for this Computation to
	  * be complete.
	  * @return The threads waiting on this Computation.
	  */
	public Enumeration getWaiting() { return waiting.elements(); }
	/** Add a thread to the list of those waiting on this Computation.
	  * @param t A thread waiting for this Computation to finish.
	  */
	public void addWaiting( Thread t ) { waiting.add(t); }
	/** Indicate a thread is no longer waiting for this Computation.
	  * @param t A thread that is currently waiting for the Computation.
	  */
	public void removeWaiting( Thread t ) { waiting.remove(t); }
	/** Clear the list of threads waiting for this Computation. */
	public void clearWaiting() { waiting.clear(); }
	/** Interrupt the threads waiting for this Computation.
	  * Probably will indicate that the Computation is complete and a result
	  * has been received from the server.
	  */
	public void interruptWaiting() {
		Iterator iter = waiting.iterator();
		while( iter.hasNext() ) {
			((Thread)iter.next()).interrupt();
		}
	}

	private String matrixURL = null;
	/** Get the URL of the matrix data for this Computation.
	  * @return The URL of the matrix data for this Computation.
	  */
	public String getMatrixURL() { return matrixURL; }
	/** Set the URL of the matrix data for this Computation.
	  * @param url The matrix URL.
	  */
	public void setMatrixURL( String url ) { matrixURL = url; }

	private File matrixFile = null;
	/** Get the local file where the matrix is stored (if any)
	  * @return The local file where the matrix is stored, or null if the
	  *         matrix is not stored locally.
	  */
	public File getMatrixFile() { return matrixFile; }
	/** Set the file location of the matrix data, stored locally.
	  * @param f The location of the matrix on this machine.
	  */
	public void setMatrixFile( File f ) {
		matrixFile = f;
		setMatrixURL( baseAddress + "/" + getId() + ".sms" );
	}

	private File archive = null;
	/** Save this Computation to an XML archive in the appropriate folder.
	  */
	public void saveArchive() {
		if( archive == null )
			archive = new File
			    ( computationDir, String.valueOf( id ) + ".xml" );
		XMLEncoder xout = null;
		try {
			xout = new XMLEncoder( new FileOutputStream(archive) );
			xout.writeObject( this );
		} catch( IOException e ) { archive = null; }
		finally { if( xout != null ) xout.close(); }
	}

	private String vectorURL = null;
	/** Get the URL of the vector data for this Computation.
	  * @return The URL of the vector data for this Computation, or null if
	  *         this computation does not use vector data.
	  */
	public String getVectorURL() { return vectorURL; }
	/** Set the URL of the vector data for this Computation.
	  * @param url The vector URL.
	  */
	public void setVectorURL( String url ) { vectorURL = url; }
	
	private long prime = -1;
	/** Get the prime modulus to use for this computation.
	  * -1 Indicates computation over the integers.
	  * @return The prime modulus for this computation.
	  */
	public long getPrime() { return prime; }
	/** Set the prime modulus to use for this computation.
	  * @param p The prime to use, or -1 to compute over the integers.
	  */
	public void setPrime( long p ) { prime = p; }

	private String result = null;
	/** Get the server's response for this Computation.
	  * @return The result of the Computation.
	  */
	public String getResult() { return result; }
	/** Set the server's response for this Computation.
	  * @param r The result of the computation.
	  */
	public void setResult( String r ) {
		if( replyReceivedOn == null ) replyReceivedOn = new Date();
		result = r;
	}
	
	private int id = -1;
	/** Get the id of this Computation
	  * @return The unique nonnegative integer id for this computation.
	  */
	public int getId() {
		if( id == -1 ) setId();
		return id;
	}
	/** Set the id of this Computation to a specific value.
	  * Only actually changes the id value if it has not yet been set.
	  * @param i The id for this Computation.
	  */
	public void setId( int i ) {
		if( id == -1 ) {
		// only set the ID if it has not yet been set.
			id = i;
			Integer idObj = new Integer(id);
			if( !all.containsKey(idObj) ) 
				all.put( idObj, this );
		}
	}
	/** Automatically set the id of this Computation.
	  * Does nothing if the ID has already been set.
	  */
	public void setId() {
	    if( id != -1 ) return;
	    Integer temp = null;
	    boolean fail = true;
	    while(fail) {
	    	temp = new Integer(generator.nextInt(Integer.MAX_VALUE));
		synchronized( all ) {
		    if( !all.containsKey( temp ) ) {
		    	all.put( temp, this );
			fail = false;
		    }
		}
	    }
	    id = temp.intValue();
	}

	/** Write the pertinent information about this Computation.
	  * Most likely this will be connected to the server performing the
	  * Computation.
	  * @param out The stream to which to write.
	  */
	public void write( java.io.DataOutputStream out ) 
		throws java.io.IOException
	{
		out.writeInt( getId() );
		out.writeInt( getType() );
		out.writeUTF( getMatrixURL() );
		if( getType() == SOLVE ) out.writeUTF( getVectorURL() );
		out.writeLong( getPrime() );
	}

	/** Default constructor. */
	public Computation() {}
	
}
