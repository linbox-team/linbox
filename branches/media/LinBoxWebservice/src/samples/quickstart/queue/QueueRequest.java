// QueueRequest.java
// Matthew Fendt Summer '07-'08
// Contains the information to hold a web service request

package samples.quickstart.middleman;
 
import java.rmi.Remote;
import java.rmi.RemoteException;

public interface QueueRequest extends Remote
{
    int enqueue(String matrixOp, String matrix, String IPAddress)
	throws RemoteException;

    String pollAnswer(int ID) throws RemoteException;

    double[] getTime(int ID) throws RemoteException;

    int countOps() throws RemoteException;

    int checkID(int ID) throws RemoteException;
}
