import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;
import java.net.*;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/25/01
 * Last Updated: 11/14/01
 * File: ResultServlet.java
 * 
 * This file contains a servlet which is called when the user types a
 * number into Result.html and hits "submit" to retrieve the answer(s) to
 * the GAP/Homology commands.  This servlet opens the correct output file,
 * reads in the data, and calls Message.DisplayGapResult to display the 
 * results to the user in HTML.
 */
public class ResultServlet extends HttpServlet{

    // Data member to hold the random number in a string from Result.html
    // See JavaDOC below for description of String contents
    String stringRandomInt;
    

    /**
     * @param HttpServletRequest request
     * @param HttpServletResponse response
     * @return none
     *
     * This method processes the Result.html file and displays the results
     * to the user, as mentioned above.
     * If the key is only a number, it is from the Homology server.
     * If, however, the key begins with "SM" and then has a number,
     * it is from the Smith Form Server.
     */
    public void doPost(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{
	
	// the User inputted number for the filename
	String stringRandomInt = request.getParameter("randomint");
	
	//System.out.println("Trying to input key: " + stringRandomInt);

	// The conversion of stringRandomInt to an actual Integer type
	int randomInt;
	
	// In this try block, only the Homology key will work.  This
	// is because the Smith Form Server key is not a straight 
	// integer (with "SM") so this try block contains the code
	// to process and display the results of the Homology server
	try{
	    // Convert the string to an int
	    randomInt = Integer.parseInt(stringRandomInt);
	    
	    // Create a variable with the filename and open the file to read
	    String temp = "Out" + randomInt + ".rslt";
	    File fileResults = new File("Gap/Output", temp);

	    // Call method to get the results stored in to a string
	    String results = getResultsFromFile(fileResults);

	    // Display the Homology Results
	    Message.DisplayGapResult(response, results, "");
	    
	} // end of try block
	
	// This catch block is where the Smith Form results will fall,
	// since during the attempt to convert the string to an Integer
	// it will fail and end up here.  It differentiates between
	// a correct Smith Form key and invalid keys.
	//	catch (NumberFormatException nfe){
	catch (FileNotFoundException fnfe){
	    try{
		// Check the string to make sure it starts with SM
		// and find the correct output file
		StringBuffer sbTemp = new StringBuffer("");
		sbTemp.append("SMOut").append(stringRandomInt.substring(0) );
		sbTemp.append(".rslt");
		//System.out.println("KEY: " + sbTemp.toString() );

		// Create a variable with the filename and open the file
		// to read
		File fileResults = new File("/home/schrag/SmithForm/Output",
					    sbTemp.toString() );
		
		// Call method to get the results stored in a string
		String results = getResultsFromFile(fileResults);

		// Display the Smith Form results to the user
		Message.DisplayGapResult(response, results, "");
	    }

	    // If it is an invalid key, a FileNotFoundException will
	    // occur, and an error message will be displayed to the user
	    catch(FileNotFoundException fnf){
		Message.DisplayResultError(response);
	    }
	}

	// If the number is invalid (not found), display an error message
	//	catch (FileNotFoundException fnfe){
	    // Message.DisplayResultError(response);
	    //	}
	
    } // end of method doPost

    /**
     * @param File fileResults
     * @return String
     *
     * This method takes a File which contains the GAP output and 
     * opens it and stores all the contents (the results) into a String
     * which is returned to the doPost method for later display
     */
    private String getResultsFromFile(File fileResults) throws IOException{
	FileReader outResults = new FileReader(fileResults);
	
	// Create a buffer to read in the file
	BufferedReader br = new BufferedReader(outResults);
	String s;  // a temp string to read in each line
	
	// A buffer to hold all the data read in
	StringBuffer sbData = new StringBuffer("");
	
	// Read in the whole file (all the results)
	while ((s = br.readLine()) != null){
	    sbData.append(s).append("<br>");
	}

	// return the string containing the results
	return sbData.toString();

    } // end of method getResultsFromFile
    
} // end of class ResultServlet
