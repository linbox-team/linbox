// QueueEntry.java
// Matthew Fendt Summer '07-'08 Research
// A data entry in the LinBoxQueue that contains information about the web
// service request

//package samples.quickstart.queue;

public class QueueEntry
{
    String matrixOperation;
    String matrix;
    String submittedTime;
    String submitterIPAddress;
    String estimatedComputationTime;
    String estimatedReturnTime;
    QueueEntry previous;
    QueueEntry next;
    int operationID;

    public QueueEntry(String mO, String m, String sT,String IP, int ID) 
    {
	matrixOperation = mO;
	matrix = m;
	submittedTime = sT;
	submitterIPAddress = IP;
	estimatedComputationTime = "Unknown";
	estimatedReturnTime = "Unknown";
	previous = null;
	next = null;
	operationID = ID;
    }

    // ---------------------------------------------------

    public void setEstimatedComputationTime(String time) 
    {estimatedComputationTime = time;}

    public void setEstimatedReturnTime(String time)
    {estimatedReturnTime = time;}

    public void setNext(QueueEntry n)  {next = n;}

    public void setPrevious(QueueEntry p)  {previous = p;}

    // ---------------------------------------------------

    public String getMatrixOperation()  {return matrixOperation;}

    public String getMatrix()  {return matrix;}

    public String getSubmittedTime()  {return submittedTime;}

    public String getSubmitterIPAddress()  {return submitterIPAddress;}

    public String getEstimatedComputationTime()  
    {return estimatedComputationTime;}

    public String getEstimatedReturnTime()  {return estimatedReturnTime;}

    public QueueEntry getNext()  {return next;}

    public QueueEntry getPrevious()  {return previous;}

    public int getOperationID()  {return operationID;}
}
