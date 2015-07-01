package LinBoxServerInterface;

import java.io.File;
import java.util.Enumeration;
import java.util.Timer;
import java.util.TimerTask;

/** Objects of this class traverse both main memory and specific temporary
  * folders on the disk in order to find and destroy old Computations and
  * matrices, so that the server can keep running forever (theoretically).
  */
class ComputationSweeper extends TimerTask
			 implements Locations, SweeperConstants
{
	private static float hoursToDelete = defaultHoursToDelete;
	/** Get the number of hours a Computation is idle before it is deleted.
	  * @return The max idle time in hours.
	  */
	public static float getHoursToDelete() { return hoursToDelete; }
	/** Set the number of hours a Computation may be idle before it will be
	  * deleted.
	  * @param h The new max idle time, in hours.
	  */
	public static void setHoursToDelete( float h ) {
		hoursToDelete = h;
		updateMillisToDelete();
	}

	private static long millisToDelete;
	private static void updateMillisToDelete() {
		millisToDelete = (long)(hoursToDelete * (60*60*1000));
	}
	static { updateMillisToDelete(); }

	private static int sweepPeriodMins = defaultSweepPeriodMins;
	/** Get the number of minutes between automatic sweepings.
	  * Note that this does not mean that any sweeps are scheduled.
	  * @return The time between sweeps, in minutes.
	  */
	public static int getSweepPeriodMins() { return sweepPeriodMins; }
	/** Set the number of minutes between automatic sweepings.
	  * Only future schedulings of sweeps will be affected.
	  * @param m The time between sweeps, in minutes.
	  */
	public static void setSweepPeriodMins( int m ) {
		sweepPeriodMins = m;
		updateSweepPeriod();
	}

	private static long sweepPeriod;
	private static void updateSweepPeriod() {
		sweepPeriod = ((long)sweepPeriodMins) * 60 * 1000;
	}
	static { updateSweepPeriod(); }

	private static ComputationSweeper instance = null;
	private static Timer t = new Timer();

	/** Start a new schedule of automatic sweepings.
	  * If a schedule is already started, it will be stopped.
	  */
	public static synchronized void scheduleSweeps() {
		if( instance != null ) instance.cancel();
		instance = new ComputationSweeper();
		t.scheduleAtFixedRate( instance, 0, sweepPeriod );
	}

	/** Cancel all scheduled sweeps, if any. */
	public static synchronized void cancelSweep() {
		if( instance != null ) instance.cancel();
		instance = null;
	}

	/** Perform a sweep operation immediately, and in the current thread. */
	public static void sweepNow() {
		(new ComputationSweeper()).run();
	}

	/** Perform a sweep operation. */
	public void run() {
		// Sweep computations in memory
		Enumeration comps = Computation.getAll();
		Computation thisOne;
		while( comps.hasMoreElements() ) {
			thisOne = (Computation)comps.nextElement();
			if( ( System.currentTimeMillis() -
			      thisOne.getLastAccessed().getTime() ) >
			    millisToDelete )
				thisOne.tryDestroy();
		}

		// Sweep computation (.xml) files
		File[] list = computationDir.listFiles
			( new ExtensionFilter(".xml") );
		String temp;
		int id = -1;
		for( int i = 0; i < list.length; ++i ) {
			temp = list[i].getName();
			try {
			    id = Integer.parseInt
				(temp.substring( 0, temp.lastIndexOf(".xml") ));
			} catch( NumberFormatException e ) { continue; }
			if( Computation.getByID(id) == null )
				list[i].delete();
		}

		// Sweep matrix (.sms) files
		list = matrixDir.listFiles( new ExtensionFilter(".sms") );
		for( int i = 0; i < list.length; ++i ) {
			temp = list[i].getName();
			try {
			    id = Integer.parseInt
				(temp.substring( 0, temp.lastIndexOf(".sms") ));
			} catch( NumberFormatException e ) { continue; }
			if( Computation.getByID(id) == null )
				list[i].delete();
		}
	}
}
