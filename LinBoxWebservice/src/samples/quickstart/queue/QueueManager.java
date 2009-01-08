// QueueManager.java
// Matthew Fendt Summer '07-'08 Research
// Receives calls from the 'Middleman' web service.  Enqueues these requests,
// calls the appropriate 'Transfer Agent' web service, and gives the answers
// back to the 'Middleman' web service

import java.rmi.registry.Registry;
import java.rmi.registry.LocateRegistry;
import java.rmi.RemoteException;
import java.rmi.server.UnicastRemoteObject;
import java.lang.Thread;
import java.io.*;
import java.net.*;

import samples.quickstart.middleman.QueueRequest;
import samples.quickstart.clients.OperationTimeEstimator;

import samples.quickstart.service.adb.TransferAgentStub;
import samples.quickstart.service.adb.TransferAgentCallbackHandler;

public class QueueManager implements QueueRequest
{
    // ***********************************************************************
    // The queue for holding requests
    static LinBoxQueue q;

    // Seconds since an enqueue request
    static int secsSinceCall = 0;

    // The maximum number of answers that can be stored
    static final int SIZE = 100;

    // Stores the answers until the 'Middleman' web service picks them up
    static String[] answers;

    // Stores the estimated time to complete the operation
    static double[] timeEstimate;

    // Keeps track of if we finished the computation
    static int[] finished;
    
    // A connection to the 'Transfer Agent' web service
    static TransferAgentStub stub;

    // The number of operations in the queue;
    static int numOps;

    // The machine specs
    //    static int[] specs;

    // For quickly estimating computation time
    static OperationTimeEstimator ote;

    // Are we currently processing a request
    static boolean processing = false;

    // ***********************************************************************

    // Constructor class
    public QueueManager() 
    {
	// Creates the queue where the request will be stored until they are
	// ready to be processed
	q = new LinBoxQueue();	
	
	numOps = 0;


	try {
	    // Create the stub that will connect to the LinBox web service
	    stub = new TransferAgentStub("http://hmrg.pc.cis.udel.edu:2000/axis2/services/TransferAgent");

	    // Create an operation time estimator object
	    ote = new OperationTimeEstimator();

	/*
	    // For getting the machine specs from the computing machine
	    TransferAgentTransferAgentSOAP11PortStub.GetMachineSpecsResponse 
		res = stub.getMachineSpecs();

	    specs = res.get_return();


	*/
	}

	catch (Exception e) { }

    }
    
    // Returns the number of operations in the queue
    public int countOps()
    {
	return numOps;
    }

    // Returns the number of operations in the queue
    public int checkID(int ID)
    {
	return q.checkID(ID);
    }

    // ***********************************************************************
    public int enqueue(String matrixOp, String matrix, String IPAddress)
    {
	// Insert the request into the queue and get back the ID of the request
	int ID = q.insertEntry(matrixOp, matrix, IPAddress);

	// Reset the time since an enqueue request
	secsSinceCall = 0;

	numOps++;

	finished[ID] = 0;

	// Do a quick time estimate

	String estimate = ote.estimateTime(matrixOp, matrix, false);

	timeEstimate[ID] = 
	    Double.parseDouble(estimate);


	QueueEntry t = q.findEntry(ID);

	t.setEstimatedComputationTime(estimate);

	System.out.println("Estimate for " + ID + ":" + estimate);

	System.out.println("**************************************");
	System.out.println("New entry into queue:\n");
	q.printLast();
	System.out.println("**************************************");


	// Gives back the request ID to the 'Middleman' web service
	return (ID); 
    }

    // ***********************************************************************
    private static void dequeue()
    {
	// Increase the timeout
	long soTime = 60 * 60 * 1000; // 1 hour
	stub._getServiceClient().getOptions().setTimeOutInMilliSeconds(soTime);

	// When there is something in the queue, wait for a little 
	// while and then remove it from the queue

	try {
	    Thread.currentThread().sleep(1000 * 5);

	    System.out.println("**************************************");
	    System.out.println("About to dequeue.  Printing contents of queue before:");
	    q.printQueue();		

	    
	    System.out.println("**************************************");
	    System.out.println("Removing entry " + 
			       q.tail.getOperationID() + " from queue.");
	    System.out.println("**************************************");

	    // Where in the answers array the result will be stored
	    final int answerLocation = q.tail.getOperationID();

	    // Remove the request from the queue and store the operation and
	    // matrix
	    String[] s = q.removeEntry();	
	    numOps--;

	    System.out.println("**************************************");
	    System.out.println("Printing contents of queue after:");
	    q.printQueue();		


	    // If the operation was rank, call the rank method from the
	    // 'Transfer Agent' web service
	    if (s[0].equalsIgnoreCase("rank"))
		{
		    // We are processing a request
		    processing = true;

		    // First, get a time estimate
		    getTimeEstimate(s[1], "rank", answerLocation);

		    // Creates a 'Rank' object
		    TransferAgentStub.Rank req = 
			new TransferAgentStub.Rank();

		    // Set the parameter- the user's matrix
		    req.setMatrix(s[1]);

		    // The callback handler that will listen for the web
		    // service to give back an answer
		    TransferAgentCallbackHandler 
			callBackHandler=
			new TransferAgentCallbackHandler()
			{
			    public void receiveResultrank(TransferAgentStub.RankResponse res)
			    {
				// Store the user's answer
				answers[answerLocation] = res.get_return();

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

				System.out.println("Answer to operation " + 
						   answerLocation + " was " +
						   res.get_return());

				synchronized (this) {
				    this.notifyAll();
				}
			    }
			    
			    public void receiveErrorrank(Exception e)
			    {
				e.printStackTrace();

				answers[answerLocation] = e.getMessage();
				    //"Error in handler";

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

			    }
			    
			};
		    // Call the web service
		    stub.startrank(req, callBackHandler);	
		}

	    // If the operation was determinant, call the determinant method 
	    // from the 'Transfer Agent' web service
	    if (s[0].equalsIgnoreCase("determinant"))
		{
		    // We are processing a request
		    processing = true;

		    // Creates a 'Determinant' object
		    TransferAgentStub.Determinant req = 
		      new TransferAgentStub.Determinant();

		    // Set the parameter- the SOAP11Ps matrix
		    req.setMatrix(s[1]);

		    // The callback handler that will listen for the web
		    // service to give back an answer
		    TransferAgentCallbackHandler 
			callBackHandler=
			new TransferAgentCallbackHandler()
			{
			    public void receiveResultdeterminant(TransferAgentStub.DeterminantResponse res)
			    {
				// Store the user's answer
				answers[answerLocation] = res.get_return();

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

				// We are no longer processing a reuqest
				processing = false;

				System.out.println("Answer to operation " + 
						   answerLocation + " was " +
						   res.get_return());

				synchronized (this) {
				    this.notifyAll();

				}
			    }
			    
			    public void receiveErrordeterminant(Exception e)
			    {
				System.err.println("Error in QueueManager:dequeue:determinant: " + e.toString());
				e.printStackTrace();
				answers[answerLocation] = "Error in handler";

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;
			    }
			    
			};

		    // Call the web service
		    stub.startdeterminant(req, callBackHandler);	
		}

	    // If the operation was trace, call the trace method from the
	    // 'Transfer Agent' web service
	    if (s[0].equalsIgnoreCase("trace"))
		{
		    // We are processing a request
		    processing = true;

		    // Creates a 'Trace' object
		    TransferAgentStub.Trace req = 
			new TransferAgentStub.Trace();

		    // Set the parameter- the user's matrix
		    req.setMatrix(s[1]);

		    // The callback handler that will listen for the web
		    // service to give back an answer
		    TransferAgentCallbackHandler 
			callBackHandler=
			new TransferAgentCallbackHandler()
			{
			    public void receiveResulttrace(TransferAgentStub.TraceResponse res)
			    {
				// Store the user's answer
				answers[answerLocation] = res.get_return();

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

				System.out.println("Answer to operation " + 
						   answerLocation + " was " +
						   res.get_return());

				synchronized (this) {
				    this.notifyAll();


				}
			    }
			    
			    public void receiveErrortrace(Exception e)
			    {
				answers[answerLocation] = "Error in handler";
			    
				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;
			    }
			    
			};

		    // Call the web service
		    stub.starttrace(req, callBackHandler);	
		}

	    // If the operation was valence, call the valence method from the
	    // 'Transfer Agent' web service
	    if (s[0].equalsIgnoreCase("valence"))
		{
		    // We are processing a request
		    processing = true;

		    // Creates a 'Valence' object
		    TransferAgentStub.Valence req = 
			new TransferAgentStub.Valence();

		    // Set the parameter- the user's matrix
		    req.setMatrix(s[1]);

		    // The callback handler that will listen for the web
		    // service to give back an answer
		    TransferAgentCallbackHandler 
			callBackHandler=
			new TransferAgentCallbackHandler()
			{
			    public void receiveResultvalence(TransferAgentStub.ValenceResponse res)
			    {
				// Store the user's answer
				answers[answerLocation] = res.get_return();

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

				System.out.println("Answer to operation " + 
						   answerLocation + " was " +
						   res.get_return());

				synchronized (this) {
				    this.notifyAll();


				}
			    }
			    
			    public void receiveErrorvalence(Exception e)
			    {
				answers[answerLocation] = "Error in handler";

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;
			    }
			    
			};

		    // Call the web service
		    stub.startvalence(req, callBackHandler);	
		}

	    // If the operation was SNF, call the SNF method from the
	    // 'Transfer Agent' web service
	    if (s[0].equalsIgnoreCase("smithNormalForm"))
		{
		    // We are processing a request
		    processing = true;

		    // Creates a 'SNF' object
		    TransferAgentStub.SmithNormalForm 
			req = new 
		      TransferAgentStub.SmithNormalForm();

		    // Set the parameter- the user's matrix
		    req.setMatrix(s[1]);

		    // The callback handler that will listen for the web
		    // service to give back an answer
		    TransferAgentCallbackHandler 
			callBackHandler=
			new TransferAgentCallbackHandler()
			{
			    public void receiveResultsmithNormalForm(TransferAgentStub.SmithNormalFormResponse res)
			    {
				// Store the user's answer
				answers[answerLocation] = res.get_return();

				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;

				System.out.println("Answer to operation " + 
						   answerLocation + " was " +
						   res.get_return());

				synchronized (this) {
				    this.notifyAll();
				}
			    }
			    
			   public void receiveErrorsmithNormalForm(Exception e)
			    {
				answers[answerLocation] = "Error in handler";
				System.out.println("Error in handler:SNF");
			    
				// We finished the computation
				finished[answerLocation] = 1;

				// We are no longer processing a reuqest
				processing = false;
			    }
			    
			};

		    // Call the web service
		    stub.startsmithNormalForm(req, callBackHandler);	
		}	    
	}
	
	catch (Exception e) 
	    {
		System.out.println("Error in dequeue");
		e.printStackTrace();
	    }
    }

    // ***********************************************************************
    // Returns the answer to the user
    public String pollAnswer(int ID)
    {
	return answers[ID];
    }


    // ***********************************************************************
    // Returns the time estimate to the user
    public double[] getTime(int ID)
    {
	// The first number will be the estimate for that one particular 
	// operation.  The second number will be an estimate of how long until
	// that operation will be done.
	double response[] = new double[2];

	// Return the latest time estimate
	response[0] = timeEstimate[ID];

	// Now find the total length of all of the operations in front of this
	// one.
	double totalTime = q.totalComputationTime(ID);

	// If the request is still in the queue, then set the remaining time
	if (totalTime != 0)
	    response[1] = totalTime;

	// If the response is not in the queue and still being worked on, the
	// time left is the time to do the operation
	else if (totalTime == 0 && finished[ID] == 0)
	    response[1] = response[0];

	// If we have already finished the computation, there is no time
	// remaining
	else if (totalTime == 0 && finished[ID] == 1)
	    response[1] = 0;

	return response;
    }

    // ***********************************************************************

    public static void getTimeEstimate(String matrix, String operation, int ID)
    {
	try 
	    {

		if (operation.equalsIgnoreCase("rank"))
		    {
			// Create an EstimateRankTime object
			TransferAgentStub.EstimateRankTime
			    req = new 
			    TransferAgentStub.
			    EstimateRankTime();
			
			// Set the parameters
			req.setMatrix(matrix);
			
			// Create the response and get the estimate
			TransferAgentStub.
			    EstimateRankTimeResponse
			    res = stub.estimateRankTime(req);
			
			double time = res.get_return();
			
			// Store the estimate
			timeEstimate[ID] = time;
		    }
	    }

	catch (Exception e)
	    {
		
	    }
    }


    public static void main(String[] args)
    {
	try{
	    QueueManager m = new QueueManager();

	    // Makes the object available through RMI
	    QueueRequest stub = 
		(QueueRequest) UnicastRemoteObject.exportObject(m, 0);

	    // Creates the answer array
	    answers = new String[SIZE];
	    for (int i = 0; i < SIZE; i++)
		{
		    answers[i] = null;
		}

	    // Creates the time estimate array
	    timeEstimate = new double[SIZE];
	    for (int i = 0; i < SIZE; i++)
		{
		    timeEstimate[i] = -2;
		}


	    // Creates the finished array
	    finished = new int[SIZE];
	    for (int i = 0; i < SIZE; i++)
		{
		    finished[i] = -1;
		}

	    // Bind the remote object's stub in the registry
	    Registry registry = LocateRegistry.getRegistry();
	    registry.bind("QueueRequest", stub);


	    System.err.println("Server ready");

	    // Starts listening for enqueues	    
	    while (true)
	    {				
		// If there is nothing in the queue or we are already 
		// processing,, sleep for a while and then check again
		while (q.tail == null || processing == true)
		    {
			Thread.currentThread().sleep(1000);
			secsSinceCall++;

			if (secsSinceCall % 5 == 0)
			{
			    System.out.println("It has been " + 
				       secsSinceCall + 
				       " seconds since a call to the queue.");
			    continue;
			}
		    }

		// When there are entries in the queue that are ready to be
		// processed, dequeue them
		dequeue();		
	    }
	}
	    
	catch(Exception e) 
	    {
		System.err.println("Server excpetion: " + e.toString());
		e.printStackTrace();
	    }
    }
}
