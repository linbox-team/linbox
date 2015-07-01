// Matthew Fendt, Summer '07-'08 Research
// Middleman for the LinBox web service.  This will take requests from the
// clients, queue them, eventually send the computations to LinBox, and then
// send the answers back to the client

package samples.quickstart.middleman;

import java.io.*;
import java.net.*;
import java.lang.Thread;

import java.rmi.registry.LocateRegistry;
import java.rmi.registry.Registry;
import java.rmi.*;

import org.apache.axis2.AxisFault;
import org.apache.axis2.Constants;
import org.apache.axis2.addressing.EndpointReference;
import org.apache.axis2.client.Options;
import org.apache.axis2.client.ServiceClient;
import org.apache.axis2.client.async.AsyncResult;
import org.apache.axis2.client.async.Callback;
import org.apache.axis2.description.AxisService;

public class TransferAgentMiddlemanSkeleton implements TransferAgentMiddlemanSkeletonInterface{

    public GetTimeEstimateResponse getTimeEstimate(GetTimeEstimate param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Get the answer back
	    double[] answer = r.getTime(requestID);

	    // Create a new response that can be called by the client
	    GetTimeEstimateResponse res = new GetTimeEstimateResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    GetTimeEstimateResponse res = new GetTimeEstimateResponse();

	    // Set the return answer value to be the answer that was computed
	    double[] answer = {-2, -2};
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    GetTimeEstimateResponse res = new GetTimeEstimateResponse();

	    // Set the return answer value to be the answer that was computed
	    double[] answer = {-2, -2};
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    GetTimeEstimateResponse res = new GetTimeEstimateResponse();

	    // Set the return answer value to be the answer that was computed
	    double[] answer = {-2, -2};
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}
    }


    public EnqueueAndGetIDResponse enqueueAndGetID(EnqueueAndGetID param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");


	    String IP = param0.getRequestIP();
	    String op = param0.getOperation();
	    String matrix = param0.getMatrix();
	    int requestID;

	    if (op.equalsIgnoreCase("rank"))
		requestID = r.enqueue("Rank", matrix, IP);

	    else if (op.equalsIgnoreCase("determinant"))
		requestID = r.enqueue("Determinant", matrix, IP);

	    else if (op.equalsIgnoreCase("valence"))
		requestID = r.enqueue("Valence", matrix, IP);

	    else if (op.equalsIgnoreCase("trace"))
		requestID = r.enqueue("Trace", matrix, IP);

	    else if (op.equalsIgnoreCase("smithNormalForm"))
		requestID = r.enqueue("SmithNormalForm", matrix, IP);

	    else 
		requestID = -11;
		
	    // Create a new response that can be called by the client
	    EnqueueAndGetIDResponse res = 
		new EnqueueAndGetIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(requestID);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    EnqueueAndGetIDResponse res = 
		new EnqueueAndGetIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-2);

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	   EnqueueAndGetIDResponse res = 
		new EnqueueAndGetIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-3);

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    EnqueueAndGetIDResponse res = 
		new EnqueueAndGetIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-4);

	    // Return the response object
	    return res;
	}
    }





    public RequestrankResponse requestrank(Requestrank param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Wait for the answer back
	    String answer = r.pollAnswer(requestID);
	    while (answer == null)
		{
		    Thread.currentThread().sleep(1000);
		    answer = r.pollAnswer(requestID);
		}

	    // Create a new response that can be called by the client
	    RequestrankResponse res = new RequestrankResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    RequestrankResponse res = new RequestrankResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Registry reference could not be created. Possible problem: 'rmiregistry' command not executed on middleman side.");

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    RequestrankResponse res = new RequestrankResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Error in registry lookup or unbind. Possible problem: 'QueueManager' class not executed on middleman side");

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    RequestrankResponse res = new RequestrankResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Unknown error in web service");

	    // Return the response object
	    return res;
	}
    }


    public RequestdeterminantResponse requestdeterminant(Requestdeterminant param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Wait for the answer back
	    String answer = r.pollAnswer(requestID);
	    while (answer == null)
		{
		    Thread.currentThread().sleep(1000);
		    answer = r.pollAnswer(requestID);
		}

	    // Create a new response that can be called by the client
	    RequestdeterminantResponse res = new RequestdeterminantResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    RequestdeterminantResponse res = new RequestdeterminantResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Registry reference could not be created. Possible problem: 'rmiregistry' command not executed on middleman side.");

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    RequestdeterminantResponse res = new RequestdeterminantResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Error in registry lookup or unbind. Possible problem: 'QueueManager' class not executed on middleman side");

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    RequestdeterminantResponse res = new RequestdeterminantResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Unknown error in web service");

	    // Return the response object
	    return res;
	}
    }

    public RequesttraceResponse requesttrace(Requesttrace param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Wait for the answer back
	    String answer = r.pollAnswer(requestID);
	    while (answer == null)
		{
		    Thread.currentThread().sleep(1000);
		    answer = r.pollAnswer(requestID);
		}

	    // Create a new response that can be called by the client
	    RequesttraceResponse res = new RequesttraceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    RequesttraceResponse res = new RequesttraceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Registry reference could not be created. Possible problem: 'rmiregistry' command not executed on middleman side.");

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    RequesttraceResponse res = new RequesttraceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Error in registry lookup or unbind. Possible problem: 'QueueManager' class not executed on middleman side");

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    RequesttraceResponse res = new RequesttraceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Unknown error in web service");

	    // Return the response object
	    return res;
	}
    }

    public RequestvalenceResponse requestvalence(Requestvalence param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Wait for the answer back
	    String answer = r.pollAnswer(requestID);
	    while (answer == null)
		{
		    Thread.currentThread().sleep(1000);
		    answer = r.pollAnswer(requestID);
		}

	    // Create a new response that can be called by the client
	    RequestvalenceResponse res = new RequestvalenceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    RequestvalenceResponse res = new RequestvalenceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Registry reference could not be created. Possible problem: 'rmiregistry' command not executed on middleman side.");

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    RequestvalenceResponse res = new RequestvalenceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Error in registry lookup or unbind. Possible problem: 'QueueManager' class not executed on middleman side");

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    RequestvalenceResponse res = new RequestvalenceResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Unknown error in web service");

	    // Return the response object
	    return res;
	}
    }



    public RequestsmithNormalFormResponse requestsmithNormalForm
	(RequestsmithNormalForm param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    int requestID = param0.getRequestID();

	    // Wait for the answer back
	    String answer = r.pollAnswer(requestID);
	    while (answer == null)
		{
		    Thread.currentThread().sleep(1000);
		    answer = r.pollAnswer(requestID);
		}

	    // Create a new response that can be called by the client
	    RequestsmithNormalFormResponse res = 
		new RequestsmithNormalFormResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    RequestsmithNormalFormResponse res = 
		new RequestsmithNormalFormResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Registry reference could not be created. Possible problem: 'rmiregistry' command not executed on middleman side.");

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    RequestsmithNormalFormResponse res = 
		new RequestsmithNormalFormResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Error in registry lookup or unbind. Possible problem: 'QueueManager' class not executed on middleman side");

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    RequestsmithNormalFormResponse res = 
		new RequestsmithNormalFormResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return("Unknown error in web service");

	    // Return the response object
	    return res;
	}
    }

    public CountOperationsInQueueResponse countOperationsInQueue()
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    // Get the number of operations in the queue
	    int count = r.countOps();

	    // Create a new response that can be called by the client
	    CountOperationsInQueueResponse res = 
		new CountOperationsInQueueResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(count);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    CountOperationsInQueueResponse res = 
		new CountOperationsInQueueResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    CountOperationsInQueueResponse res = 
		new CountOperationsInQueueResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    CountOperationsInQueueResponse res = 
		new CountOperationsInQueueResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}
    }


    public CheckIDResponse checkID(CheckID param0)
    {
	try {
	    // Look up the queue in the registry
	    Registry registry = LocateRegistry.getRegistry();

	    // Formulate a request
	    QueueRequest r = (QueueRequest)registry.lookup("QueueRequest");

	    // Get the number of operations in front of this request
	    int count = r.checkID(param0.getID());

	    // Create a new response that can be called by the client
	    CheckIDResponse res = 
		new CheckIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(count);

	    // Return the response object
	    return res;
	}

	catch(RemoteException e) {
	    // Create a new response that can be called by the client
	    CheckIDResponse res = 
		new CheckIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}

	catch(NotBoundException e) {
	    // Create a new response that can be called by the client
	    CheckIDResponse res = 
		new CheckIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}

	catch(Exception e) {
	    // Create a new response that can be called by the client
	    CheckIDResponse res = 
		new CheckIDResponse();

	    // Set the return answer value to be the answer that was computed
	    res.set_return(-1);

	    // Return the response object
	    return res;
	}
    }
}
