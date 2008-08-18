// Matthew Fendt, Summer '07 Research
// Client that accesses the LinBox Webservice.  A copy of this class will need 
// to be on the user end of the service.

package samples.quickstart.clients;
import samples.quickstart.service.adb.TransferAgentTransferAgentHttpportStub;
import samples.quickstart.service.adb.TransferAgentTransferAgentHttpportCallbackHandler;
import java.io.*;
import java.net.*;

import javax.naming.Context;
import org.apache.axis2.context.ConfigurationContext;
import org.apache.axis2.context.ConfigurationContextFactory;
import org.apache.axis2.AxisFault;
import org.apache.axis2.Constants;
import org.apache.axis2.addressing.EndpointReference;
import org.apache.axis2.client.Options;
import org.apache.axis2.client.ServiceClient;
import org.apache.axis2.client.async.AsyncResult;
import org.apache.axis2.client.async.Callback;
import org.apache.axis2.description.AxisService;
import org.apache.axis2.transport.jms.JMSConstants;

// ONLY FOR A TEST
import java.util.Date;
import java.text.DateFormat;
import java.text.SimpleDateFormat;


public class TransferAgentStandAloneClient {

    public static void main(String[] args)
    {
	try
	    {
		// Creates the stub that will connect to the Linbox web service
		// running on port 2000 of hmrg
		TransferAgentTransferAgentHttpportStub stub = new TransferAgentTransferAgentHttpportStub("http://hmrg.pc.cis.udel.edu:2000/axis2/services/TransferAgent");

		// Calls the web service with the stub and the command line 
		// arguments, namely, the operation and the matrix file
		request(stub, args);
	    }

	catch(Exception e)
	    {
		e.printStackTrace();
	    }
    }

    // Performs the correct opertion of the web service and displays the answer
    // to the user
    public static void request(TransferAgentTransferAgentHttpportStub stub, String[] args)
    {
	try{
	    // The matrix has to be sent to hmrg not as a file, but as a stream
	    // This code opens the matrix file that the user specifies and 
	    // converts it to a UTF8 stream
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
	    String str = buffer.toString();
	    System.out.println(str);

	    // Converts the matrix to UTF8 encoding
	    str = URLEncoder.encode(str, "UTF-8");

	    // Optional, displays what the matrix file looks like in UTF8
	    System.out.println(str);


	    // ***************************************************************
	    // If the desired operation is 'rank,' calls the Rank function
	    // of the stub class asynchronously
	    // ***************************************************************
	    if (args[0].equalsIgnoreCase("rank")){

		// ONLY FOR A TEST
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date d1 = new Date();
		System.out.println("Time before call: " + df.format(d1));
		// END ONLY FOR A TEST


		// Creates a 'stub' object
		TransferAgentTransferAgentHttpportStub tastub= new
		    TransferAgentTransferAgentHttpportStub();

		// Creates a 'Rank' object
		TransferAgentTransferAgentHttpportStub.Rank req = 
		    new TransferAgentTransferAgentHttpportStub.Rank();

		// Set the parameter of the 'Rank' object- the matrix
		req.setMatrix(str);

		// Extends the default 'Callback Handler'.  This will be
		// listening for the result to be asynchronously sent back
		// from the web service and will print the result to the
		// screen
		TransferAgentTransferAgentHttpportCallbackHandler 
		    callBackHandler = 
		    new TransferAgentTransferAgentHttpportCallbackHandler() {
			public void receiveResultrank(TransferAgentTransferAgentHttpportStub.RankResponse res)
			{
			    System.out.println("Rank is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorrank(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		// Submit the asynchronous request; the answer will be sent
		// back some time later
		tastub.startrank(req, callBackHandler);

		// ONLY FOR A TEST
		Date d2 = new Date();
		System.out.println("Time after call: " + df.format(d2));
		// END ONLY FOR A TEST

		// We will just wait for the request to get sent back to us
		// because this is a simple program, but we could continue
		// execution and eventually the Callback Handler will get the
		// callback from the web service
		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) {}
		}
	    }

	    // ***************************************************************
	    // If the desired operation is 'determinant,' calls the Determinant
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("determinant")){

		// ONLY FOR A TEST
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date d1 = new Date();
		System.out.println("Time before call: " + df.format(d1));

		TransferAgentTransferAgentHttpportStub tastub= new
		    TransferAgentTransferAgentHttpportStub();

		TransferAgentTransferAgentHttpportStub.Determinant req = 
		    new TransferAgentTransferAgentHttpportStub.Determinant();
		req.setMatrix(str);

		TransferAgentTransferAgentHttpportCallbackHandler 
		    callBackHandler = 
		    new TransferAgentTransferAgentHttpportCallbackHandler() {
			public void receiveResultdeterminant(TransferAgentTransferAgentHttpportStub.DeterminantResponse res)
			{
			    System.out.println("Determinant is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrordeterminant(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startdeterminant(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) {}
		}
		
		// ONLY FOR A TEST
		Date d2 = new Date();
		System.out.println("Time after call: " + df.format(d2));

	    }
	    
	    // ***************************************************************
	    // If the desired operation is 'valence,' calls the Valence
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("valence")){

		// ONLY FOR A TEST
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date d1 = new Date();
		System.out.println("Time before call: " + df.format(d1));

		TransferAgentTransferAgentHttpportStub tastub= new
		    TransferAgentTransferAgentHttpportStub();

		TransferAgentTransferAgentHttpportStub.Valence req = 
		    new TransferAgentTransferAgentHttpportStub.Valence();
		req.setMatrix(str);

		TransferAgentTransferAgentHttpportCallbackHandler 
		    callBackHandler = 
		    new TransferAgentTransferAgentHttpportCallbackHandler() {
			public void receiveResultvalence(TransferAgentTransferAgentHttpportStub.ValenceResponse res)
			{
			    System.out.println("Valence is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorvalence(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startvalence(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) {}
		}
		       
		// ONLY FOR A TEST
		Date d2 = new Date();
		System.out.println("Time after call: " + df.format(d2));

	    }

	    // **************************************************************  
	    // If the desired operation is 'trace,' calls the Trace
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("trace")){

		// ONLY FOR A TEST
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date d1 = new Date();
		System.out.println("Time before call: " + df.format(d1));

		TransferAgentTransferAgentHttpportStub tastub= new
		    TransferAgentTransferAgentHttpportStub();

		TransferAgentTransferAgentHttpportStub.Trace req = 
		    new TransferAgentTransferAgentHttpportStub.Trace();
		req.setMatrix(str);

		TransferAgentTransferAgentHttpportCallbackHandler 
		    callBackHandler = 
		    new TransferAgentTransferAgentHttpportCallbackHandler() {
			public void receiveResulttrace(TransferAgentTransferAgentHttpportStub.TraceResponse res)
			{
			    System.out.println("Trace is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrortrace(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.starttrace(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) {}
		}
		       
		// ONLY FOR A TEST
		Date d2 = new Date();
		System.out.println("Time after call: " + df.format(d2));

	    }


	    // ***************************************************************
	    // If the desired operation is 'SNF,' calls the SNF
	    // function of the stub class asynchronously
	    // ***************************************************************
	    else if (args[0].equalsIgnoreCase("smithNormalForm")){

		// ONLY FOR A TEST
		DateFormat df = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date d1 = new Date();
		System.out.println("Time before call: " + df.format(d1));

		TransferAgentTransferAgentHttpportStub tastub= new
		    TransferAgentTransferAgentHttpportStub();

		TransferAgentTransferAgentHttpportStub.SmithNormalForm req = 
		    new TransferAgentTransferAgentHttpportStub.SmithNormalForm();
		req.setMatrix(str);


		TransferAgentTransferAgentHttpportCallbackHandler 
		    callBackHandler = 
		    new TransferAgentTransferAgentHttpportCallbackHandler() {
			public void receiveResultsmithNormalForm(TransferAgentTransferAgentHttpportStub.SmithNormalFormResponse res)
			{
			    System.out.println("SNF is:    " + 
					       res.get_return());
			    synchronized (this) {
				this.notifyAll();
			    }
			}

			public void receieveErrorsmithNormalForm(Exception e) {
			    System.out.println("Error has occured");
			}
		    };

		tastub.startsmithNormalForm(req, callBackHandler);

		synchronized (callBackHandler) {
		    try{
			// Wait for response
			callBackHandler.wait();
		    }  catch (InterruptedException e) {}
		}
		       
		// ONLY FOR A TEST
		Date d2 = new Date();
		System.out.println("Time after call: " + df.format(d2));

	    }

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
