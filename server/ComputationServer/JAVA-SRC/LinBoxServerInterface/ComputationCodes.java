package LinBoxServerInterface;

/** An interface to store codes for the various computation operations, as well
  * as codes to perform actions on the server (i.e. NEW, DONE, KILL).
  */
public interface ComputationCodes {
	int NEW = 1;
	int DONE = 2;
		// followed by id, result
	int KILL = 3;
		// followed by id
	int RANK = 4;
	int DET = 5;
	int SOLVE = 6;	// includes vector specification
	int SOLVE0 = 7;
	int CHARPOLY = 8;
	int MINPOLY = 9;
	int SMITH = 10;
}
