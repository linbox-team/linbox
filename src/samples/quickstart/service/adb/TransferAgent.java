// TransferAgent.java
// Matthew Fendt, Summer '07 Research
// A POJO that will take a matrix method and a matrix file, compute the answer,
// and return the data when requested.  Used with the LinBox web service.

// @@@@@@@@@@@@ Currently supported methods: 
// 1) determinant
// 2) rank
// 3) trace

package samples.quickstart.service.adb;

// The LinBox shared library
import samples.quickstart.service.adb.xsd.*;

public class TransferAgent
{
    //-------------------------------------------------------------------------
    
    public String rank(String matrix)
    {
	// Load the LinBox shared library
	try {java.lang.System.loadLibrary("linboxfunctions");}
	
	// If the library cannot be loaded, then that will cause an error
	catch (UnsatisfiedLinkError e)
	    {
		return "Error in loading linbox library";
	    }

	//	return linboxfunctions.rankFiles("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/" + matrix);
	return linboxfunctions.rankFiles(matrix);

    }

    //-------------------------------------------------------------------------
    public String determinant(String matrix)
    {
	// Load the LinBox shared library
	try {java.lang.System.loadLibrary("linboxfunctions");}
	
	// If the library cannot be loaded, then that will cause an error
	catch (UnsatisfiedLinkError e)
	    {
		return "Error in loading linbox library";
	    }

	String answer = linboxfunctions.detFiles(matrix);

	// If the number is too big to store
	// @@ This should not be an issue
	/*	if (answer > 2147483647)
	    {
		System.out.println("Answer is too large");
		
	    }
	*/
	return answer;

    }
    //-------------------------------------------------------------------------

    public String valence(String matrix)
    {
	// Load the LinBox shared library
	try {java.lang.System.loadLibrary("linboxfunctions");}
	
	// If the library cannot be loaded, then that will cause an error
	catch (UnsatisfiedLinkError e)
	    {
		return "Error in loading linbox library";
	    }

	return linboxfunctions.valFiles(matrix);

    }

    //-------------------------------------------------------------------------

    public String trace(String matrix)
    {
	// Load the LinBox shared library
	try {java.lang.System.loadLibrary("linboxfunctions");}
	
	// If the library cannot be loaded, then that will cause an error
	catch (UnsatisfiedLinkError e)
	    {
		return "Error in loading linbox library";
	    }

	return linboxfunctions.traceFiles(matrix);

    }

    //-------------------------------------------------------------------------

 } // End TransferAgent class 
