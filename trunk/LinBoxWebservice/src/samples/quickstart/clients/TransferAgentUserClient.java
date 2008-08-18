// Matthew Fendt, Summer '07 - '08 Research
// Calls the intermidate client that the user will utilize to make LinBox web 
// service calls

package samples.quickstart.clients;


import samples.quickstart.middleman.TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub;
import samples.quickstart.service.adb.TransferAgentTransferAgentHttpportStub;
import samples.quickstart.middleman.TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler;
import java.io.*;
import java.net.*;

import org.apache.axis2.AxisFault;
import org.apache.axis2.Constants;
import org.apache.axis2.addressing.EndpointReference;
import org.apache.axis2.client.Options;
import org.apache.axis2.client.ServiceClient;
import org.apache.axis2.client.async.AsyncResult;
import org.apache.axis2.client.async.Callback;
import org.apache.axis2.description.AxisService;

public class TransferAgentUserClient {

    public static void main(String[] args)
    {
	try
	    {
		// Connects to the web service
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub stub = new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub("http://hmrg.pc.cis.udel.edu:2000/axis2/services/TransferAgentMiddleman");
		
		long soTime = 60 * 60 * 1000; // 1 hour
		stub._getServiceClient().getOptions().
		    setTimeOutInMilliSeconds(soTime);

		// Call the appropriate web service
		request(stub, args);
		
	    }
	
	
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
	
    }

    public static void request(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub tastub, String[] args)
    {
	try{
	    String str = null;

	    try{
		// The matrix has to be sent to LinBox not as a file but as a 
		// stream. This code opens the matrix file that the user 
		// specifies and converts it to a UTF8 stream
		StringBuffer buffer = new StringBuffer();
		FileInputStream fis = new FileInputStream(args[1]);

		InputStreamReader isr = new InputStreamReader(fis, "UTF8");
		Reader in = new BufferedReader(isr);
		int ch;
		
		while ((ch = in.read()) > -1)
		    {
			buffer.append((char)ch);
		    }
		in.close();
		
		// Optional, displays the matrix to the screen
		str = buffer.toString();
		// System.out.println(str);
		
		// Converts the matrix to UTF8 encoding
		str = URLEncoder.encode(str, "UTF-8");
		
		// Optional, displays what the matrix file looks like in UTF8
		//System.out.println(str);
	    }

	    catch (Exception e) 
		{ 
		    if (!(args[0].equalsIgnoreCase("queuestatus") ||
			  (args[0].equalsIgnoreCase("checkID")) ||
			  (args[0].equalsIgnoreCase("getTimeEstimate")))) 
			{
			    System.out.println("Couldn't open file");
			    System.exit(0);
			}
		}

	    // ***************************************************************
	    // ************* RANK ********************************************
	    // ***************************************************************
	    if (args[0].equalsIgnoreCase("rank")){

	       // Creates an ID_Request object
	       TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID ereq = new 
		   TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID();

	      // Set the ID_Request parameters
	      ereq.setMatrix(str);
	      ereq.setOperation("rank");

	      // Create the response and get the ID
	      TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		  EnqueueAndGetIDResponse eres = tastub.enqueueAndGetID(ereq);
	      int ID = eres.get_return();

	      System.out.println("Your ID is: " + ID + "\n");


	      if (ID > 0) {
		// Creates a 'Rank' object
	      final TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		    Requestrank req = new 
		    TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		    Requestrank();

		// Set the parameter of the 'Rank' object- the ID
		req.setRequestID(ID);
		
		// Extends the default 'Callback Handler'.  This will be
		// listening for the result to be asynchronously sent back
		// from the web service and will print the result to the
		// screen
		final TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler 
		    callBackHandler = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler() {
			public void receiveResultrequestrank(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestrankResponse res)
			{
			    System.out.println("Rank is:    " + 
					       res.get_return());

			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrequestrank(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		// Submit the asynchronous request; the answer will be sent
		// back some time later
		tastub.startrequestrank(req, callBackHandler);

		// We will just wait for the request to get sent back to us
		// because this is a simple program, but we could continue
		// execution and eventually the Callback Handler will get the
		// callback from the web service
		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) 
			{System.out.println("Error in callback");}
		}
	      }
	    }

	    // ***************************************************************
	    // If the desired operation is 'determinant,' calls the Determinant
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("determinant")){

	       // Creates an ID_Request object
	       TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID ereq = new 
		   TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID();

	      // Set the ID_Request parameters
	      ereq.setMatrix(str);
	      ereq.setOperation("determinant");

	      // Create the response and get the ID
	      TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		  EnqueueAndGetIDResponse eres = tastub.enqueueAndGetID(ereq);
	      int ID = eres.get_return();

	      System.out.println("Your ID is: " + ID + "\n");


	      if (ID != -1) {

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requestdeterminant req = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requestdeterminant();
		req.setRequestID(ID);

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler 
		    callBackHandler = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler() {
			public void receiveResultrequestdeterminant(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestdeterminantResponse res)
			{
			    System.out.println("Determinant is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrequestdeterminant(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startrequestdeterminant(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (Exception e) 
			{	
			    e.printStackTrace();
			}
		}
	      }
	    }
	    
	    // ***************************************************************
	    // If the desired operation is 'valence,' calls the Valence
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("valence")){

	       // Creates an ID_Request object
	       TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID ereq = new 
		   TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID();

	      // Set the ID_Request parameters
	      ereq.setMatrix(str);
	      ereq.setOperation("valence");

	      // Create the response and get the ID
	      TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		  EnqueueAndGetIDResponse eres = tastub.enqueueAndGetID(ereq);
	      int ID = eres.get_return();

	      System.out.println("Your ID is: " + ID + "\n");


	      if (ID != -1) {
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requestvalence req = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requestvalence();
		req.setRequestID(ID);

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler 
		    callBackHandler = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler() {
			public void receiveResultrequestvalence(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestvalenceResponse res)
			{
			    System.out.println("Valence is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrequestvalence(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startrequestvalence(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) 
			{
			    e.printStackTrace();
			}
		}
	      }
	    }

	    // **************************************************************  
	    // If the desired operation is 'trace,' calls the Trace
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("trace")){

	       // Creates an ID_Request object
	       TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID ereq = new 
		   TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID();

	      // Set the ID_Request parameters
	      ereq.setMatrix(str);
	      ereq.setOperation("determinant");

	      // Create the response and get the ID
	      TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		  EnqueueAndGetIDResponse eres = tastub.enqueueAndGetID(ereq);
	      int ID = eres.get_return();

	      System.out.println("Your ID is: " + ID + "\n");


	      if (ID != -1) {
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requesttrace req = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.Requesttrace();

		req.setRequestID(ID);

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler 
		    callBackHandler = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler() {
			public void receiveResultrequesttrace(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequesttraceResponse res)
			{
			    System.out.println("Trace is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrequesttrace(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startrequesttrace(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) 
			{
			    e.printStackTrace();
			}
		}
	      }
		       
	    }


	    // ***************************************************************
	    // If the desired operation is 'SNF,' calls the SNF
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("smithNormalForm")){

	       // Creates an ID_Request object
	       TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID ereq = new 
		   TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		   EnqueueAndGetID();

	      // Set the ID_Request parameters
	      ereq.setMatrix(str);
	      ereq.setOperation("smithNormalForm");

	      // Create the response and get the ID
	      TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.
		  EnqueueAndGetIDResponse eres = tastub.enqueueAndGetID(ereq);
	      int ID = eres.get_return();

	      System.out.println("Your ID is: " + ID + "\n");


	      if (ID != -1) {
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestsmithNormalForm req = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestsmithNormalForm();

		req.setRequestID(ID);

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler 
		    callBackHandler = 
		    new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortCallbackHandler() {
			public void receiveResultrequestsmithNormalForm(TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.RequestsmithNormalFormResponse res)
			{
			    System.out.println("SNF is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrequesttrace(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startrequestsmithNormalForm(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) 
			{
			    e.printStackTrace();
			}
		}

	      }
	    }
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("queueStatus")){

		// Get the answer back
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.CountOperationsInQueueResponse res = tastub.countOperationsInQueue();

		System.out.println("There are " + res.get_return() +
				   " operations currently in the queue");


	    }

	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("checkID")){
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.CheckID req = new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.CheckID();
		req.setID(Integer.parseInt(args[1]));

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.CheckIDResponse res = tastub.checkID(req);

		System.out.println("There are " + res.get_return() +
				   " operations in line in front of yours.");

	    }


	    
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("getTimeEstimate")){
		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.GetTimeEstimate req = new TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.GetTimeEstimate();
		req.setRequestID(Integer.parseInt(args[1]));

		TransferAgentMiddlemanTransferAgentMiddlemanSOAP11PortStub.GetTimeEstimateResponse res = tastub.getTimeEstimate(req);

		if (res.get_return()[0] < 0)
		    {
			System.out.println("There was a problem getting the estimate");
		    }

		else
		    {
			System.out.println("Estimated time to do computation: "
					   + res.get_return()[0] + "\n" + 
					   "Estimated total wait time: " + 
					   res.get_return()[1]);
		    }
	    }

	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("estimateLength")){
		OperationTimeEstimator ote = new OperationTimeEstimator();

		String answer = ote.estimateTime(args[2], str, true);

		System.out.println("Estimated time to compute " + args[2] + 
				   ": " + answer);
	    }
	    // ***************************************************************


	    // If the user asks for anything else, reply that that operation
	    // is not supported
	    else
		{
		    System.out.println("That operation is not currently supported");
		}
	}
	
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
    }
}

