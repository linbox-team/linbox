// LinBoxQueue.java
// Matthew Fendt Summer '07-'08 Research
// A PERSISTANT queue that will hold web service requests.

import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;

public class LinBoxQueue
{
    public static QueueEntry head;
    public QueueEntry tail;
    public int operationNumber;

    public LinBoxQueue()
    {
	head = null;
	tail = null;
	operationNumber = 0;
    }

    public int insertEntry(String matrixOperation, String matrix, 
			    String submitterIPAddress)
    {
	// The ID of the request
	operationNumber++;

	// For recording the time when the request is received
	DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
	Date d = new Date();

	// Enqueue the request
	QueueEntry e = new QueueEntry(matrixOperation, matrix, 
				      dateFormat.format(d), submitterIPAddress,
				      operationNumber);
	
	// Set the entries' previous and next pointers
	e.setPrevious(null);


	if (head == null)
	    {
		head = e;
		tail = e;
	    }

	else
	    {
		e.setNext(head);
		head.setPrevious(e);
		head = e;
	    }

	return operationNumber;
    }

    // Finds a given entry, given the ID
    public static QueueEntry findEntry(int ID)
    {
	QueueEntry temp;

	for (temp = head; temp != null; temp = temp.getNext())
	    {
		if (temp.getOperationID() == ID)
		    {
			return temp;
		    }
	    }

	return null;
    }

    // Removes the last entry in the queue
    String[] removeEntry()
    {
	String[] s = new String[2];
	s[0] = tail.getMatrixOperation();
	s[1] = tail.getMatrix();

	if (head != tail)
	    {
		tail = tail.getPrevious();
		tail.setNext(null);
	    }

	// There is only one node left

	else if (head == tail)
	    {
		head = null;
		tail = null;
	    }

	return s;
    }

    // Prints the contents of the queue
    void printQueue()
    {
	QueueEntry temp;

	for (temp = head; temp != null; temp = temp.getNext())
	    {
		System.out.println("Operation " + temp.getOperationID() + 
				   "\n\tSubmitted " + temp.getSubmittedTime() +
				   "\n\tEstimated computation time: " + 
				   temp.getEstimatedComputationTime() + 
				   "\n\tEstimated return time: " + 
				   temp.getEstimatedReturnTime() + 
				   "\n\tOperation: "+temp.getMatrixOperation()+
				   "\n\tSubmitter IP Address: " + 
				   temp.getSubmitterIPAddress() + 
				   "\n\n\n");
	    }
    }

    // Prints the latest addition to the queue
    void printLast()
    {
	QueueEntry temp = head;

	System.out.println("Operation " + temp.getOperationID() + 
			   "\n\tSubmitted " + temp.getSubmittedTime() +
			   "\n\tEstimated computation time: " + 
			   temp.getEstimatedComputationTime() + 
			   "\n\tEstimated return time: " + 
			   temp.getEstimatedReturnTime() + 
			   "\n\tOperation: "+temp.getMatrixOperation()+
			   "\n\tSubmitter IP Address: " + 
			   temp.getSubmitterIPAddress() + 
			   "\n");
    }

    // Checks how many operations are in line in front of a given operation
    int checkID(int ID)
    {
	QueueEntry temp;
	int count = -1;

	for (temp = head; temp != null; temp = temp.getNext())
	    {
		// Once we find the ID in the list, start keeping track of
		// subsequent operations
		if (temp.getOperationID() == ID)
		    count = 0;

		// If we started keeping track of the operations, increment
		else if (count >= 0)
		    count++;

	    }

	return count;
    }

    // Find the total estimated time of all operations in front of a given
    // operation
    double totalComputationTime(int ID)
    {
	QueueEntry temp;
	double total = 0;

	for (temp = findEntry(ID); temp != null; temp = temp.getNext())
	    {
		String length = temp.getEstimatedComputationTime();
		total += Double.parseDouble(length);
	    }

	return total;
    }
}
