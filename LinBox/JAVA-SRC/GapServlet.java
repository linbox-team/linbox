import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;
import java.util.Random;
import java.util.Date;
import java.net.*;


/**
 * @author Brian Gold
 * @modified by Eric Schrag
 * @version 1.05
 * 
 * Date: 7/25/01
 * Last Updated: 7/3/02
 * File: GapServlet.java
 * 
 * This file contains the method doPost which reads the client submission
 * of the form Gap.html, and the method computeResult which processes
 * the commands found through the Gap.html form and displays the result to the
 * user.  computeResult is also responsible for the logging of the Gap/
 * Homology Server.
 */
public class GapServlet extends HttpServlet{
    
    // Data member functions, explained below in doPost
    String IPAddress;
    String EmailAddress;
    String email_bool;
    String Results_option;
    String stringMatrixURL;
    String comments;
    String timing;
    boolean boolMatrixOK;  // for the optional URL of a matrix on Gap.html
    
    static Process proc;

    /**
     * @param HttpServletRequest request
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * contains the results of the Gap computations at the current time.
     * See above comments (preceding the class definition) for more
     * information
     */
    public void doGet(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{

	// Get the string containing the number corresponding to 
	// the current file from Message.DisplayRollingPage
	String stringRandom = request.getParameter("randomint");

	try{
	    // Parse the string to obtain the number
	    // If the string does not contain a valid number,
	    // the catch block is entered, and an error message is shown
	    int randomInt = Integer.parseInt(stringRandom);

	    // Open the correct output file for reading
	    String fname = "Out" + randomInt + ".rslt";
	    File file = new File("Gap/Output", fname);
	    FileReader frOutput = new FileReader(file);

	    // Create a buffer to read in data from the output file
	    BufferedReader br = new BufferedReader(frOutput);

	    // Create a temporary string to hold each line
	    String s;

	    // Set up the printing mechanism to display Dynamic HTML
	    PrintWriter pwOutput;
	    response.setContentType("text/html");
	    pwOutput = response.getWriter();

	    // Create three StringBuffers, each to hold a different
	    // chunk of data: sbBeg holds the HTML before the current
	    // results; sbData holds the current data; and sbEnd holds
	    // the ending data of the HTML
	    StringBuffer sbBeg = new StringBuffer();
	    StringBuffer sbData = new StringBuffer();
	    StringBuffer sbEnd = new StringBuffer();
	    
	    // Set up the dynamic HTML
	    sbBeg.append("<HTML><HEAD><TITLE>\n");
	    sbBeg.append("Homology Server - Calculations Currently Running\n");
	    sbBeg.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	    sbBeg.append("<a name=\top\">");
	    sbBeg.append("<font size=\"4\"><strong>");
	    sbBeg.append("Your GAP/Homology commands are currently being ");
	    sbBeg.append("calculated.</strong></size></font><br><br>");
	    sbBeg.append("<meta http-equiv=\"refresh\" content=\"30;>");
	    sbBeg.append("<meta name=\"keywords\" content=\"automatic redirection\">");
	    sbBeg.append("<TABLE CELLSPACING=\"5\"><TR><TH></strong>Update Partial Results*</TH>");
	    sbBeg.append("<TH></strong>Stop Current Calculations</TH><TH></TH><TH>");
	    sbBeg.append("</strong>Current memory load of Server**</TH>");
	    sbBeg.append("</TR><TR><TH>");
	    sbBeg.append("<form action=GapServlet method=get>");
	    sbBeg.append("<INPUT TYPE=hidden name=\"randomint\" value=");
	    sbBeg.append(randomInt);
	    sbBeg.append(">");
	    sbBeg.append("<INPUT TYPE=submit value=\"Update Partial Results\">");
	    sbBeg.append("</form></TH><TH>");
	    sbBeg.append("<form action=KillServlet method=post>");
	    sbBeg.append("<INPUT TYPE=hidden name=\"serverID\" value=Simplicial>");
	    sbBeg.append("<INPUT TYPE=submit value=\"Kill Calculations\">");
	    sbBeg.append("</form></TH><TH></TH>");
	    sbBeg.append("<TH><form action=TopServlet method=get>");
	    sbBeg.append("<A HREF=\"TopServlet\" TARGET=\"main\">Memory Load</A></form>");
	    sbBeg.append("</TH></TR></strong><TR><font size=\"2\">");
	    sbBeg.append("<TH></strong>*Automatically updates every");
	    sbBeg.append(" 30 seconds</TH><TH></TH><TH></TH><TH></strong>");
	    sbBeg.append("**Opens new browser");
	    sbBeg.append(" window</TH></size></font></TR></TABLE>");
	    sbBeg.append("<hr><font size=\"4\"><strong>Here is the output");
	    sbBeg.append(" so far:</strong></size></font><br><br>");
	    
	    // Read in the file until #donegap is read, or the file ends
	    // (either way it reads the whole file, it just needs to
	    // differentiate between the completed results and partial results)
	    while ( (s = br.readLine()) != null){
		if (s.equals("#donegap"))
		    break;
		sbData.append(s).append("<br>");
	    }
	    sbData.append("<br><br>");

	    // If the last line is #donegap, the calculations have completed
	    // Therefore, display the final answers by calling the appropriate
	    // method, Message.DisplayGapResult
	    if (s!=null && s.equals("#donegap") ){
		Message.DisplayGapResult(response, sbData.toString(), "");
	    }

	    // Else, the results are only partial, and so therefore, 
	    // setup the rest of the page, including the automatic refresh
	    // and the submit button
	    else{

		sbEnd.append("<br><br><br>");
		sbEnd.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
		sbEnd.append("<br><a href=SmithForm.html>To the Smith Form Server</a>");
		sbEnd.append("<br><br><br>Thank you for using the Homology Server.");
		sbEnd.append("<br><br>top of <a href=\"#top\">page</a>");
		sbEnd.append("</BODY></HTML>");
		
		sbBeg.append(sbData).append(sbEnd);
		
		// Display the partial results to the user
		pwOutput.println(sbBeg.toString() );
		pwOutput.flush();
		pwOutput.close();
	    }
	    
	    
	    
	}  // end of try block
	
	// If for some reason the string does not contain a valid number
	// which should NOT occur at all, this will be displayed on the server
	// (not to the user)  But again, I note, that if this occurs, something
	// is wrong elsewhere in the code, very wrong and should be fixed.
	catch(NumberFormatException nfe){
	    System.out.println("Number Format Exception in GeneratedServlet");
	}
    } // end of method doGet

    
    /**
     * @param HttpServletRequest request
     * @param HttpServletResponse response
     * @return none
     *
     * This method processes the form Gap.html submitted from the user
     * and updates and makes the necessary log files.  It also runs
     * GAP and feeds it the commands and takes the results from GAP
     */
    public void doPost(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{
	
	// the User inputted url of the GAP commands
	String stringURL = request.getParameter("URLcommands");

	// user input of actions to perform
	String stringCommands = request.getParameter("gapcommands");

	// User selection of how to handle the comments
	comments = request.getParameter("comments");

        // User selection of the timing
        timing = request.getParameter("timing");

	// get the IP Address of the user
	IPAddress = request.getRemoteAddr();

	// get the E-mail address from the user to send a notification to
	EmailAddress = request.getParameter("email");

	// Does the user want the results included in the email? Yes or no
	email_bool = request.getParameter("emailbool");

	// does the user want the results, or wait for email?
	Results_option = request.getParameter("results_option");

	// the User inputted url of the matrix, optional
	stringMatrixURL = request.getParameter("URLmatrix");

	// default setting is false, the optional matrix textfield will
	// be checked shortly
	boolMatrixOK = false;
	

	URLHandler urlhURL = new URLHandler(stringURL);
	URLHandler urlhMatrix = new URLHandler(stringMatrixURL);
	

	// BEGIN THE ERROR CHECKING OF THE FORM

	// Only one of both the URL text field and the text box
	// may be filled out.  Otherwise it is an error.

	// If neither are filled out, display an error message
	if (ErrorCheck.CheckCommands(stringCommands) == false &&
	    ErrorCheck.CheckURLError(urlhURL) == false){  
	    Message.DisplayBothError(response);
	}
	
	// If both are filled out, don't know which to process
	// so display an error message
	else if ( (ErrorCheck.CheckURLError(urlhURL) == true) &&
		  (ErrorCheck.CheckCommands(stringCommands) == true) ){
	    Message.DisplayNoneError(response);
	}
	
	// If the email notification button has been checked and 
	// the email address has an error, display an error message
	else if (Results_option.equals("email") &&
		 ErrorCheck.CheckEmailError(EmailAddress) == true){
	    Message.DisplayEmailError(response);
	}
	
	// The user has submitted a URL of GAP commands
	else if (ErrorCheck.CheckURLError(urlhURL) == false){
	    // if the user has also submitted a matrix URL and that
	    // has errors, display an error message
	    if (ErrorCheck.CheckUrlDataError(urlhURL) == true){
		Message.DisplayURLErrorMessage(response);
	    }
	    
	    // check the optional matrix URL for errors
	    else if (ErrorCheck.CheckURLError(urlhMatrix) == false){
		// If the matrix URL has errors, display an error message
		if (ErrorCheck.CheckMatrixDataError(urlhMatrix) == true)
		    Message.DisplayURLErrorMessage(response);
		else
		    // set the flag to true
		    boolMatrixOK = true;
	    }
	    
	    // If the matrix URL is invalid, display an error message
	    else if ( !(stringMatrixURL.equals("")==false ||
			stringMatrixURL != null ||
			(stringMatrixURL.equals("http://my-site/~me-public/matrix-file") == false) ) &&
		      (ErrorCheck.CheckURLError(urlhMatrix) == true) ){

		Message.DisplayURLErrorMessage(response);
	    }
	    
	    // else, compute the results of the commands in the URL
	    computeResult(response, urlhURL.data, urlhMatrix.data);
	}
	

        // Check user commands for improper use of commands
        else if (ErrorCheck.commandCheck(stringCommands) == true){
	    Message.DisplayInvalidCommand(response);
	}


	// If the user has entered commands into the text box, process them
	else if (ErrorCheck.CheckCommands(stringCommands) == false){
	    computeResult(response, stringCommands, urlhMatrix.data);
	}
	
	// Should never get to this point, something is missing if it does
	else{
	    System.out.println("Something is wrong here.");
	}
	
    } // end of method doPost
    

    /**
     * @param HttpServletResponse response
     * @param String data
     * @param String matrixdata
     * @return none
     *
     * This method is called by doPost() above and actually creates and writes
     * the log files and runs GAP and processes the results.  Basically,
     * this is the main event, doing all the processing, handing off
     * to utility functions when necessary (like Message class, etc)
     */
    public void computeResult(HttpServletResponse response, String data,
			      String matrixfiledata)
	throws IOException{
	
        StringBuffer tempData = new StringBuffer("RequirePackage(\"homology\"); ");
        tempData.append("\n");        
        tempData.append(GapCalculations.getCommentaryCommand(comments) );
        tempData.append("\n");
        tempData.append(GapCalculations.getTimingCommand(timing) );
        tempData.append("\n");
        tempData.append(data);
        tempData.append("\n");
        tempData.append("quit;");
        tempData.append("\n");

	// append "quit;" to end of data so that Gap exits when done
	// processing/calculating all the commands
	//data += "\nquit;\n";
	
        data = tempData.toString();

	if (ErrorCheck.commandCheck(data) == true)
	    Message.DisplayInvalidCommand(response);


	try{
	    
	    // If there is a URL with a matrix entered into the optional
	    // text field, save it as a local file so that it can be read
	    // correctly by GAP.  It is saved merely as the last part of the
	    // URL, everything after the last '/'
// 	    if (boolMatrixOK){
// 		FileWriter outMatrix = new FileWriter(
// 		      stringMatrixURL.substring( 
// 				stringMatrixURL.lastIndexOf('/')+1) );
// 		outMatrix.write(matrixdata);
// 		outMatrix.flush();
// 		outMatrix.close();
// 	    }
	    
	    // Get the current date and time
	    Date date = new Date();
	    Random random = new Random();
	    int randomInt;
	    String fname;
	    File fileCommands;

	    // Get a random Integer and see if a log file has already
	    // been created based off of that integer.  If it has been
	    // used before, get a different random Integer;  When
	    // one has not been used before, make a new file based on that
	    // random Integer, as each GAP submission receives its own
	    // entry into the master log, along with its very own Input
	    // and Output files, for the commands and results, respectively.
	    // Input files are of the format "In" + randomInt + ".gap"
	    // Output files are of the format "Out" + randomInt + ".rslt"
	    do{
		randomInt = random.nextInt(999999);
		fname = "In" + randomInt + ".gap";
		fileCommands = new File("Gap/Input", fname);
	    }while(fileCommands.exists() == true);
	    
	    



	    // Create and write to the input log file
	    FileWriter outCommands = new FileWriter(fileCommands);
	    // write timestamp, ipaddress, email address, and commands
	    // to the input log file
	    outCommands.write(date.toString() );
	    outCommands.write("\n");
	    outCommands.write(IPAddress);
	    outCommands.write("\n");
	    outCommands.write(EmailAddress);
	    outCommands.write("\n");
	    outCommands.write(data);
	    outCommands.close();

	    String fresults = "Out" + randomInt + ".rslt";

	    // Append to master log file, gap.log
	    // Write the following: date/time, IPAddress, Email address,
	    // filename for commands (input file), and filename for results
	    // (output file)
	    FileWriter outSystemLog = new FileWriter("Gap/gap.log", true);
	    outSystemLog.write(date.toString() );
	    outSystemLog.write("\n");
	    outSystemLog.write(IPAddress);
	    outSystemLog.write("\n");
	    outSystemLog.write(EmailAddress);
	    outSystemLog.write("\n\t");
	    outSystemLog.write(fname);
	    outSystemLog.write("\n\t");
	    outSystemLog.write(fresults);
	    outSystemLog.write("\n\n");
	    outSystemLog.close();

	    // create output log file
	    FileWriter outResult = new FileWriter("Gap/Output/Out" +
						  randomInt + ".rslt", true);
	    // write timestamp and email address to output logfile
	    outResult.write(date.toString() );
	    outResult.write("\n");
	    outResult.write(EmailAddress);
	    outResult.write("\n");

	    // Run the shell script "gap.sh" which runs GAP
	    // the "-b" option suppresses the GAP banner
	    proc = Runtime.getRuntime().exec("gap.sh -b");
            
            // Display the appropriate HTML page to the user depending
	    // on if the user wants to wait for the results or 
	    // receive an email notification (only difference is one
	    // states that an email will be sent, and the other,
	    // obviously (since there is no email) does not)
	    if (Results_option.equals("email") ){
		Message.DisplayRollingPage(response, data, EmailAddress,
					   randomInt);
	    }
	    else{
		Message.DisplayRollingPage(response, data, "", randomInt);
	    }


	    // Feed the commands into GAP
	    BufferedWriter outWriter = new BufferedWriter(new OutputStreamWriter(proc.getOutputStream()));
	    outWriter.write(data.toString() );
	    outWriter.write("\n");
	    outWriter.write("quit;");
	    outWriter.flush();


	    // temporary string for reading line by line
	    String s;
	    // output buffer for GAP results
	    StringBuffer sb = new StringBuffer("");
	    // output buffer for GAP comments
	    StringBuffer sbC = new StringBuffer();

	    StringBuffer sbFileComments = new StringBuffer("");


	    // Make a buffer to read in the comments for if the user wants
	    // the comments separate from the results
	    BufferedReader br;

	    // Get any error messages and comments and write them to
	    // the correct places in the correct order
	    if (comments.equals("none") == false &&
		timing.equals("none") == false){
		br = new BufferedReader (new InputStreamReader
		    (proc.getErrorStream() ) );
		while ((s = br.readLine()) != null){
		    
		    // write line to the output file
		    sb.append(s).append("<br>");
		    outResult.write(s);
		    outResult.write("\n");
		    outResult.flush();
		    
		    // Else, user does not want comments, so ignore
		}  // end of while loop

	    }
		outResult.flush();

	    // Set up a buffer to read the output returned by GAP
	    br = new BufferedReader (new InputStreamReader
		(proc.getInputStream() ) );

	    // More buffers to get the results of the GAP computations
	    char buf[] = new char[data.toString().length()];
	    data.toString().getChars(0, data.toString().length(), buf, 0);
	    CharArrayReader f = new CharArrayReader(buf);
	    BufferedReader brData = new BufferedReader (f);
	    String sInput = "";
	    
	    // Read all of the results outputted by GAP
	    while ((s = br.readLine()) != null){

		// Fix the format of the results if needed and output to file
	        while (s.startsWith("gap> gap>") ){
		    sInput = brData.readLine();
		    sb.append("<br>COMMAND> ").append(sInput);
		    sb.append("<br>");
		    outResult.write("<br>COMMAND> ");
		    outResult.write(sInput);
		    outResult.write("\n");
		    outResult.flush();

		    s = s.substring(5);
		}

		// Fix the format of the results if needed
		if (s.startsWith("gap> >") ){
		    String temp = s.substring(6);
		    while(temp.startsWith(" >") ){
			temp = temp.substring(2);
		    }
		    s = s.substring(0,4) + temp;
		}

		// Process the line, to a file and the StringBuffer
		if ( (s.startsWith("gap>") ) &&
		     (sInput != null) ){

		    sInput = brData.readLine();
		    sb.append("<br>COMMAND> ").append(sInput);
		    sb.append("<br>");
		    outResult.write("<br>COMMAND> ");
		    outResult.write(sInput);
		    outResult.write("\n");
		    outResult.flush();
		}
		
		// Process the line by outputting to the output file
		outResult.write(s);
		outResult.write("\n");
		outResult.flush();
		sb.append(s).append("<br>");
	    }

	    // Flush the output buffer
	    outResult.flush();


	    // If there are comments, write them to the output file
	    //if (sbFileComments.equals("")==false){
		//outResult.write("<br><br><font size=\"24\"> ");
		//outResult.write("Comments: </font>");
		//outResult.write(sbFileComments.toString() );
		//outResult.flush();
	    //}

	    // Write a string to the output file so GeneratedServlet
	    // knows that the results are complete
	    outResult.write("\n#donegap");
	    outResult.flush();
	    outResult.close();
	    
	    // If the user has chosen to be notified by email
	    // create the message and send it
	    if (Results_option.equals("email") ){
		
		StringBuffer sbEmailNotify = new StringBuffer("Your Gap commands have completed.  Please go to the following website, http://linbox.pc.cis.udel.edu:8080/gap/Result.html, and enter the following number into the text box found on that page to retrieve your results.\n\nNumber: ");
		sbEmailNotify.append(randomInt).append("\n\nThank you. --The Linbox/Gap Interface team.");
		
		if (email_bool.equals("yes") ){
		    String sTempA = sb.toString();
		    String sTempB = sbC.toString();
		    int tempindex = 0;
		    while( (tempindex=sTempA.indexOf("<br>"))!=(-1) ){
			sTempA = sTempA.substring(0, tempindex) +
			    "\n" + sTempA.substring(tempindex+4);
		    }
		    while( (tempindex=sTempB.indexOf("<br>"))!=(-1) ){
			sTempB = sTempB.substring(0, tempindex-1) +
			    "\n" + sTempA.substring(tempindex+4);
		    }
			    
		    
		    sbEmailNotify.append("\n\nThe results are listed below:\n\n").append(sTempA).append(sTempB);
		}
		
		try{
		    // Send the message
		    //System.out.println("About to send mail.");
		    SendMail.sendMail(EmailAddress, "gold@cis.udel.edu",
				      "AUTO-NOTIFY - Calculations Completed",
				      sbEmailNotify.toString() );
		    //System.out.println("mail sent???");
		}

		// Various errors in failure to send email
		catch(UnknownHostException uhe){
		    System.out.println("uhe");
		    uhe.printStackTrace();
		}
		catch(ProtocolException pe){
		    System.out.println("pe");
		}
		catch(IOException ioe){
		    System.out.println("ioe");
		}
		
	    }

	    // If the user has opted to wait for the results, 
	    // display them, now that they are ready
	    else{	
		String comments = sbC.toString();
		String result = sb.toString();
		
		Message.DisplayGapResult(response, result, comments);
	    }
	    
	} // end of try block
	
	// Catch any general errors that may occur
	catch(IOException e){
	    System.out.println("IOException");
	    e.printStackTrace();
	}

    }  // end of method computeResult()



    public static synchronized void killProcess(){
        proc.destroy();
        //System.out.println("Simplicial Process has been killed");
    }


} // end of class GapServlet


