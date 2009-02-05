// Matthew Fendt, Summer '07 Research
// Client that accesses the LinBox Webservice.  A copy of this class will need 
// to be on the user end of the service.

package samples.quickstart.clients;

import samples.quickstart.service.adb.TransferAgentStub;
import java.io.*;
import java.net.*;

public class TransferAgentDependentClient {

    public static String request(String op, String matrix, File f, 
				 int sentAsFile)
    {
	try
	    {
		// Creates the stub that will connect to the Linbox web service
		// running on port 2000 of hmrg
		TransferAgentStub stub = new TransferAgentStub("http://hmrg.pc.cis.udel.edu:2000/axis2/services/TransferAgent");

		long soTime = 60 * 60 * 1000; // 1 hour
		stub._getServiceClient().
		    getOptions().setTimeOutInMilliSeconds(soTime);

	    String str;


	    if (sentAsFile == 1) {

	    // The matrix has to be sent to hmrg not as a file, but as a stream
	    // This code opens the matrix file that the user specifies and 
	    // converts it to a UTF8 stream
	    StringBuffer buffer = new StringBuffer();
	    FileInputStream fis = new FileInputStream(f);
	    InputStreamReader isr = new InputStreamReader(fis, "UTF8");
	    Reader in = new BufferedReader(isr);
	    int ch;

	    while ((ch = in.read()) > -1)
		{
		    buffer.append((char)ch);
		}
	    in.close();

	    str = buffer.toString();
	    }

	    // If the matrix was inputted as a string
	    else 
		{
		    str = matrix;
		}


	    // Converts the matrix to UTF8 encoding
	    str = URLEncoder.encode(str, "UTF-8");


	    // If the desired operation is 'rank,' calls the Rank function
	    // of the stub class
	    if (op.equalsIgnoreCase("rank")){
		TransferAgentStub.Rank req = 
		    new TransferAgentStub.Rank();
		req.setMatrix(str);
		
		// Asks for the answer back. Note that this call is not blocked
		// so the program could run for a while before it asks for the 
		// answer.
		TransferAgentStub.RankResponse res = 
		    stub.rank(req);

		return res.get_return();
	    }
	    
	    // If the desired operation is 'determinant,' calls the Determinant
	    // function of the stub class
	    else if (op.equalsIgnoreCase("determinant")){
		TransferAgentStub.Determinant req = 
		    new TransferAgentStub.Determinant();
		req.setMatrix(str);
		
		// Asks for the anwser back. Note that this call is not blocked
		// so the program could run for a while before it asks for the 
		// answer.
		TransferAgentStub.DeterminantResponse res=
		    stub.determinant(req);

		return res.get_return();
	    }
	    
	    // If the desired operation is 'valence,' calls the Valence
	    // function of the stub class
	    else if (op.equalsIgnoreCase("valence")){
		TransferAgentStub.Valence req = 
		    new TransferAgentStub.Valence();
		req.setMatrix(str);
		
		// Asks for the answer back. Note that this call is not blocked
		// so the program could run for a while before it asks for the 
		// answer.
		TransferAgentStub.ValenceResponse res = 
		    stub.valence(req);

		return res.get_return();
	    }
	    
	    // If the desired operation is 'trace,' calls the Trace
	    // function of the stub class
	    else if (op.equalsIgnoreCase("trace")){
		TransferAgentStub.Trace req = 
		    new TransferAgentStub.Trace();
		req.setMatrix(str);
		
		// Asks for the answer back. Note that this call is not blocked
		// so the program could run for a while before it asks for the 
		// answer.
		TransferAgentStub.TraceResponse res = 
		    stub.trace(req);

		return res.get_return();
	    }

	    // If the desired operation is 'smithNormalForm,' calls the
	    // SmithNormalForm function of the stub class
	    else if (op.equalsIgnoreCase("smithNormalForm")){
		TransferAgentStub.SmithNormalForm req = 
		    new TransferAgentStub.SmithNormalForm();
		req.setMatrix(str);
		
		// Asks for the answer back. Note that this call is not blocked
		// so the program could run for a while before it asks for the 
		// answer.
		TransferAgentStub.SmithNormalFormResponse res = 
		  stub.smithNormalForm(req);

		return res.get_return();
	    }


	    // If the user asks for anything else, reply that that operation
	    // is not supported
	    else
		{
		    return "That operation is not currently supported";
		}
	    }
	
	catch(Exception e)
	    {
		e.printStackTrace();
		return "There was an error in your computation.  Some common problems may include: Linbox library not loaded correctly, client could not load file, matrix format is incorrect (eg no 0 0 0... line at the end), hmrg web service is down, etc";
	    }
    }
}


