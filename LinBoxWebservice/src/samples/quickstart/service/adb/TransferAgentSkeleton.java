// Matthew Fendt, Summer '07 Research
// This is the version of TransferAgent that works with the service.  It is 
// based on an auto generated skeleton that was built from TransferAgent.java
// This class is called by the TransferAgentClient.java that is on the 
// client's end


package samples.quickstart.service.adb;

import samples.quickstart.service.adb.xsd.*;
import java.io.*;
import java.net.*;

import java.lang.Thread;

public class TransferAgentSkeleton implements TransferAgentSkeletonInterface
{
    static 
    {
	// Loads the linbox functions that will be used
	java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");
    }

    public RankResponse rank(Rank param0)
    {
	try{
	    // ONLY FOR A TEST
	    //	    Thread.currentThread().sleep(1000 * 15);


	    // Gets the matrix from the 'Rank' object and converts it from
	    // UTF8 back to ASCII
	    String str = param0.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    String r;
	    
	    // The answer determined by the linbox function library
	    r = linboxfunctions.rankFiles(str);

	    // Create a new response that can be called by the client
	    RankResponse res = new RankResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(r);

	    // Return the response object
	    return res;
	    }

	
	catch(Exception e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		String r = "Unsupported";
		RankResponse res = new RankResponse();
		res.set_return(r);
		return res;
	    }
    }






    public EstimateRankTimeResponse estimateRankTime(EstimateRankTime param0)
    {
	try{
	    String str = param0.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    double answer;
	    
	    // The answer determined by the linbox function library
	    answer = linboxfunctions.estimateRankTime(str);

	    // Create a new response that can be called by the client
	    EstimateRankTimeResponse res = new EstimateRankTimeResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(answer);

	    // Return the response object
	    return res;
	    }

	
	catch(UnsupportedEncodingException e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		int answer = -3;
		EstimateRankTimeResponse res = new EstimateRankTimeResponse();
		res.set_return(answer);
		return res;
	    }
    }















    public DeterminantResponse determinant(Determinant param1)
    {
	try{
	    // Gets the matrix from the 'Determinant' object and converts it 
	    // from UTF8 back to ASCII
	    String str = param1.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    String d;
	    
	    // The answer determined by the linbox function library
	    d = linboxfunctions.detFiles(str);

	    // Create a new response that can be called by the client
	    DeterminantResponse res = new DeterminantResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(d);

	    // Return the response object
	    return res;
	    }

	
	catch(UnsupportedEncodingException e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		String d = "Unsupported";
		DeterminantResponse res = new DeterminantResponse();
		res.set_return(d);
		return res;
	    }
    }

    public ValenceResponse valence(Valence param2)
    {
	try{
	    // Gets the matrix from the 'Valence' object and converts it from
	    // UTF8 back to ASCII
	    String str = param2.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    String v;
	    
	    // The answer determined by the linbox function library
	    v = linboxfunctions.valFiles(str);

	    // Create a new response that can be called by the client
	    ValenceResponse res = new ValenceResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(v);

	    // Return the response object
	    return res;
	    }

	
	catch(UnsupportedEncodingException e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		String v = "Unsupported";
		ValenceResponse res = new ValenceResponse();
		res.set_return(v);
		return res;
	    }
    }

    public TraceResponse trace(Trace param3)
    {
	try{
	    // Gets the matrix from the 'Trace' object and converts it from
	    // UTF8 back to ASCII
	    String str = param3.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    String tr;
	    
	    // The answer determined by the linbox function library
	    tr = linboxfunctions.traceFiles(str);

	    // Create a new response that can be called by the client
	    TraceResponse res = new TraceResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(tr);

	    // Return the response object
	    return res;
	    }

	
	catch(UnsupportedEncodingException e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		String tr = "Unsupported";
		TraceResponse res = new TraceResponse();
		res.set_return(tr);
		return res;
	    }
    }

    public SmithNormalFormResponse smithNormalForm(SmithNormalForm param3)
    {
	try{
	    // Gets the matrix from the 'SmithNormalForm' object and converts it from
	    // UTF8 back to ASCII
	    String str = param3.getMatrix();
	    str = URLDecoder.decode(str, "US-ASCII");

	    String snf;
	    
	    // The answer determined by the linbox function library
	    snf = linboxfunctions.smithNormalFormFiles(str);

	    // Create a new response that can be called by the client
	    SmithNormalFormResponse res = new SmithNormalFormResponse();

	    // Set the return value to be the answer that was comoputed
	    res.set_return(snf);

	    // Return the response object
	    return res;
	    }

	
	catch(UnsupportedEncodingException e)
	    {
		// If the parameter is in an unsupported form, cannot perform
		// the computation
		String snf = "Unsupported";
		SmithNormalFormResponse res = new SmithNormalFormResponse();
		res.set_return(snf);
		return res;
	    }
    }










    public GetMachineSpecsResponse getMachineSpecs()
    {
	int[] specs = {2, 445, 3};


	// Create a new response that can be called by the client
	GetMachineSpecsResponse res = new GetMachineSpecsResponse();
	
	// Set the return value to be the answer that was comoputed
	res.set_return(specs);
	
	// Return the response object
	return res;
    }
}

