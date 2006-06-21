package LinBoxComputingServer;

/** Use this class to shut down a running ComputingServer */
class QuitComputing {
	/** This main method just finds a running ComputingServer and tells it
	  * to bow out gracefully.
	  * @param args Any command-line arguments will be ignored.
	  */
	public static void main( String[] args ) {
		ComputingServer.die();
	}
}
