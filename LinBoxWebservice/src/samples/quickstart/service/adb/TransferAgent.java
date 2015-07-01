// TransferAgent.java
// Matthew Fendt, Summer '07 Research
// A POJO that interacts with the actual linbox library.  Takes a matrix in the
// form of a string and returns with answer in the form of a string

// @@@@@@@@@@@@ Currently supported methods: 
// 1) determinant
// 2) rank
// 3) trace
// 4) valence
// 5) Smith normal form

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

	return linboxfunctions.rankFiles(matrix);
    }
    //-------------------------------------------------------------------------

    public double estimateRankTime(String matrix)
    {
	return -1;
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

	return linboxfunctions.detFiles(matrix);
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

    public String smithNormalForm(String matrix)
    {
	// Load the LinBox shared library
	try {java.lang.System.loadLibrary("linboxfunctions");}
	
	// If the library cannot be loaded, then that will cause an error
	catch (UnsatisfiedLinkError e)
	    {
		return "Error in loading linbox library";
	    }

	return linboxfunctions.smithNormalFormFiles(matrix);
    }

    //-------------------------------------------------------------------------

    public int[] getMachineSpecs()
    {
	int[] answer = {-1, -1};
	return answer;
    }







    //-------------------------------------------------------------------------

 } // End TransferAgent class 
