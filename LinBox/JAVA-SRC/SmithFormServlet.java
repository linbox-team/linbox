import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;
import java.util.Random;
import java.util.Date;
import java.net.*;
import java.sql.*;




/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/31/01
 * Last Updated: 2/1/02
 * File: SmithFormServlet.java
 * 
 * This file contains the method doPost which reads the client submission
 * of the form SmithForm.html, and the method computeResult which processes
 * the commands found through the SmithForm.html form and displays the result
 * to the user.  computeResult is also responsible for the logging of the 
 * Smith Form server.
 */
public class SmithFormServlet extends HttpServlet{
    
    // Data member functions, explained below in doPost
    String IPAddress;
    
    String EmailAddress;
    //String email_bool;
    String Results_option;
    //String show_comments = "Output";
    boolean boolMatrixOK;  // for the optional URL of a matrix on Gap.html
    
    String stringMatrixURL;
    String smithChoice;
    String optionalModulusM;
    boolean boolOptionalMod;
    
    String comments;
    String timing;
    String algorithm;
    
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
	    String fname = "SMOut" + randomInt + ".rslt";
	    File file = new File("SmithForm/Output", fname);
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
	    sbBeg.append("SmithForm Server - Calculations Currently Running\n");
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
	    sbBeg.append("<form action=SmithFormServlet method=get>");
	    sbBeg.append("<INPUT TYPE=hidden name=\"randomint\" value=");
	    sbBeg.append(randomInt);
	    sbBeg.append(">");
	    sbBeg.append("<INPUT TYPE=submit value=\"Update Partial Results\">");
	    sbBeg.append("</form></TH><TH>");
	    sbBeg.append("<form action=KillServlet method=post>");
	    sbBeg.append("<INPUT TYPE=hidden name=\"serverID\" value=Smith>");
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

	    sbEnd.append("<br><br><br>");
	    sbEnd.append("<a href=\"http://linalg.org/servers.html\">Back to the LinBox Server</a>");
	    sbEnd.append("<br><br><br>Thank you for using the LinBox Team Servers.");
	    sbEnd.append("<br><br>top of <a href=\"#top\">page</a>");
	    sbEnd.append("</BODY></HTML>");

	    sbBeg.append(sbData).append(sbEnd);

	    // Display the partial results to the user
	    pwOutput.println(sbBeg.toString() );
	    pwOutput.flush();
	    pwOutput.close();
	    
	    
	    
	    
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
     * This method processes the form SmithForm.html submitted from the user
     * and updates and makes the necessary log files.  It also runs
     * GAP and feeds it the commands and takes the results from GAP
     */
    public void doPost(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{
	
	// the User inputted url of the matrix
	stringMatrixURL = request.getParameter("URLmatrix");
	
        // the user's choice of calculations to perform
        smithChoice = request.getParameter("Smithchoice");
	
        // the optional modulus (will be null or "" if not entered)
        optionalModulusM = request.getParameter("modulusm");
	
	// get the IP Address of the user
	IPAddress = request.getRemoteAddr();
	
	// the user's choice of commenting
	comments = request.getParameter("comments");
	
	// the user's choice of timing
	timing = request.getParameter("timing");
	
        // the user's choice of algorithm
        algorithm = request.getParameter("algorithm");

	// the user's email address
	EmailAddress = request.getParameter("email");
	
	// selection to wait or to receive email notification
	Results_option = request.getParameter("results_option");
	
	// does the user want the results included in the email?
	//email_bool = request.getParameter("emailbool");
	
	// Check to see if a mod has been entered
        if (optionalModulusM.equals("") || optionalModulusM==null)
	    boolOptionalMod = false;
        else
	    boolOptionalMod = true;
	
        //Error check of the Matrix URL
        URLHandler urlhMatrix = new URLHandler(stringMatrixURL);
	
        if (ErrorCheck.CheckURLError(urlhMatrix) == true){
	    System.out.println("Checking URL");
	    Message.DisplayURLErrorMessage(response);
	}
	
        else if (ErrorCheck.CheckMatrixDataError(urlhMatrix) == true){
	    System.out.println("Checking data in URL");
	    Message.DisplayURLErrorMessage(response);
	}

	else if (Results_option.equals("wait")==false &&
		  ErrorCheck.CheckEmailError(EmailAddress) == true) {
	    Message.DisplayEmailError(response);
	}
	
	// If the matrix is ok (if the URL exists, basically)
        else{

	    System.out.println("No Errors in submitted form");
	    // If the user has NOT entered a modulus
	    if (boolOptionalMod==false){
		System.out.println("Computing Smith, with no modulus");
		ComputeSmithResult(response, urlhMatrix.data, false, 0);
	    }
	  
	    // If the user has entered a modulus, get the Integer and 
	    // see if it is valid, if not, catch the exception and display
	    // an error message to the user
	    else{
		
		// attempt to get the modulus
		int mod;
		
		try{
		    mod = Integer.parseInt(optionalModulusM);
		    System.out.println("Computing Smith, with mod");
		    ComputeSmithResult(response, urlhMatrix.data, true, mod);
		}
		catch(NumberFormatException nfe){
		    Message.DisplayInvalidNumberMessage(response);
		}

	    } // end of else

        } // end of else

    } // end of method doPost
    

    /**
     * @param HttpServletResponse response
     * @param String stringMatrixData
     * @param boolean boolMod
     * @int indMod
     * @return none
     * 
     * This method creates the log files through the generation of the
     * random numbers, and calls the utility functions to generate
     * the GAP commands to calculate.  This method then calls
     * the GAP shell, gap.sh to run the commands, and calls the appropriate
     * function to display the output
     */
    private void ComputeSmithResult(HttpServletResponse response,
				    String stringMatrixData, boolean boolMod,
				    int intMod) throws IOException{
	
	// Get the current date and time
	Date date = new Date();
	Random random = new Random();
	int randomInt;
	String fname;
	File fileCommands;
	Connection conn;
	Statement stmt;


	
	try{

	    Class.forName("org.gjt.mm.mysql.Driver").newInstance();

	    conn = DriverManager.getConnection("jdbc:mysql://linbox.pc.cis.udel.edu/linbox?user=root");


	    stmt = conn.createStatement();

	    stmt.executeQuery("UPDATE servers SET in_use=1 WHERE ip = '" + InetAddress.getLocalHost().getHostName() + "'");
	}

	catch(ClassNotFoundException E) {
	        System.out.println("Error in ClassCreation of driver in Servlet");
	}

	catch(SQLException E) {
	      System.out.println(E.getMessage());
	}

	catch(Exception e) {
	                System.out.println("General Error in Servlet");
	    }






	// Generate the random number
	System.out.println("Generating filename for result");
	do{
	    randomInt = random.nextInt(999999);
	    fname = "SMIn" + randomInt + ".gap";
	    fileCommands = new File("/home/schrag/SmithForm/Input", fname);
	}while(fileCommands.exists() == true);
	System.out.println("Filename Generated");
	fileCommands.createNewFile();
	System.out.println("File Created");

	
	try{
	// Create and write to the input log file
	    System.out.println("attempting to generate file");
	FileWriter outToto = new FileWriter(fileCommands);
	System.out.println("make FileWriter");
	


	// write timestamp, ipaddress, email address, and commands
	// to the input log file
	outToto.write(date.toString() );
	System.out.println("put string in buffer");
	outToto.write("\n");
	System.out.println("add newline");
	outToto.flush();
	System.out.println("flush buffer");
	outToto.close();
	System.out.println("Input log created and written");
        }

	catch(Exception e){

	    System.out.println("Something is wrong with file");
	    System.out.println(e.getMessage());

	}


	// If the user wants the results by email, display a page to the user
	//if (Results_option.equals("wait")!=false )
	//    Message.DisplayEmailPage(response, "", randomInt);
        //else
	System.out.println("Reaches strange email thing");
            Message.DisplayEmailPage(response, EmailAddress, randomInt);
	    System.out.println("Leaves strange email thing");


	// Create the file to store the actual matrix into
        String matrixFilename = stringMatrixData;              
	
	// Create and write to the input log file
	FileWriter outCommands = new FileWriter(fileCommands);
	// write timestamp, ipaddress, email address, and commands
	// to the input log file
	outCommands.write(date.toString() );
	outCommands.write("\n");
	outCommands.write(IPAddress);
	outCommands.write("\n");
	outCommands.write(stringMatrixURL);
	outCommands.write("\n");
	outCommands.write(EmailAddress);
	outCommands.write("\n");
	outCommands.flush();
	
	String fresults = "SMOut" + randomInt + ".rslt";
	
	// Append to master log file, smithform.log
	// Write the following: date/time, IPAddress, Matrix URL, Email
	// address,
	// filename for commands (input file), and filename for results
	// (output file)
	FileWriter outSystemLog = new FileWriter("/home/schrag/SmithForm/smithform.log", true);
	outSystemLog.write(date.toString() );
	outSystemLog.write("\n");
	outSystemLog.write(IPAddress);
	outSystemLog.write("\n");
	outSystemLog.write(stringMatrixURL);
	outSystemLog.write("\n");
	outSystemLog.write(EmailAddress);
	outSystemLog.write("\n\t");
	outSystemLog.write(fname);
	outSystemLog.write("\n\t");
	outSystemLog.write(fresults);
	outSystemLog.write("\n\n");
	outSystemLog.close();
	
	// create output log file
	FileWriter outResult = new FileWriter("/home/schrag/SmithForm/Output/SMOut" +
					      randomInt + ".rslt", true);
	// write timestamp and email address to output logfile
	outResult.write(date.toString() );
	outResult.write("\n");
	outResult.write(stringMatrixURL);
	outResult.write("\n");
	outResult.write(EmailAddress);
	outResult.write("\n");
	
	// Call the utility methods to create each individual
	// GAP command from the user's form submission
	StringBuffer inputCommands = GapCalculations.getCommentaryCommand(comments);
	inputCommands.append("\n");
	inputCommands.append(GapCalculations.getTimingCommand(timing) );
	inputCommands.append("\n");
	inputCommands.append(GapCalculations.getSmithCommand(matrixFilename,
					     boolMod, intMod, smithChoice, algorithm) );
	inputCommands.append("\n");
	inputCommands.append("quit;");
	inputCommands.append("\n");
	
	// Copy the commands to the input log file
	outCommands.write(inputCommands.toString() );
	outCommands.flush();
	outCommands.close();
	
	// Run the shell script "gap.sh" which runs GAP
	// the "-b" option suppresses the GAP banner
	System.out.println("Starting GAP for Smith result");
	proc = Runtime.getRuntime().exec("gap.sh -b");
	
	// Feed the commands into GAP from StringBuffer inputCommands
	BufferedWriter outWriter = new BufferedWriter(new                              OutputStreamWriter(proc.getOutputStream()));
	outWriter.write("RequirePackage(\"homology\");");
	outWriter.write("\n");
	outWriter.write(inputCommands.toString() );
	outWriter.flush();
	

	// Make a buffer to read in the comments for if the user wants
	// the comments separate from the results
	StringBuffer sbFileComments = new StringBuffer("");
	BufferedReader br;
	    
	// temporary string for reading line by line
	String s;
	// output buffer for GAP results
	StringBuffer sb = new StringBuffer("");
	// output buffer for GAP comments
	StringBuffer sbC = new StringBuffer();	
	
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
	
	// Read all of the results outputted by GAP
	while ((s = br.readLine()) != null){
	    
	    // Fix the format of the results if needed and output to file
	    while (s.startsWith("gap> gap>") ){
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
	    
	    // Process the line by outputting to the output file
	    outResult.write(s);
	    outResult.write("\n");
	    outResult.flush();
	    sb.append(s).append("<br>");
	}
	
	// Flush the output buffer
	outResult.flush();
	

	// Write a string to the output file so GeneratedServlet
	// knows that the results are complete
	outResult.write("\n#donegap");
	outResult.flush();
	outResult.close();
	
	// If the user had chosen to wait for the results, display them now
	if (Results_option.equals("wait") )
	    Message.DisplayGapResult(response, sb.toString(),
				     sbC.toString() );

	// else the user wanted email notification, so set that up
	else{
	    	StringBuffer sbEmailNotify = new StringBuffer("Your GAP commands have completed.  Please go to the following website, http://linbox.pc.cis.udel.edu:8080/gap/Result.html, and enter the following key into the text box found on that page to retrieve your results.\n\nKey: ");
		sbEmailNotify.append("SM").append(randomInt);
		sbEmailNotify.append("\n\nThank you. --The Linbox/GAP Interface team.");
		
		// If the user wants the results included in the email, do so
		if (Results_option.equals("longemail") ){
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
		
		// Either way (if the user wants the results in the email
		// or not) attempt to send the message, by calling sendMail
		// in the file SendMail.java
		try{
		    // Send the message
		    //System.out.println("About to send mail.");
		    SendMail.sendMail(EmailAddress, "schrag@cis.udel.edu",
				      "AUTO-NOTIFY - Calculations Completed",
				      sbEmailNotify.toString() );
		    //System.out.println("Mail sent");
		}

		// Various errors in failure to send email, so catch them all
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
		
	} // end of else
	

	try{

	    Class.forName("org.gjt.mm.mysql.Driver").newInstance();

	    conn = DriverManager.getConnection("jdbc:mysql://linbox.pc.cis.udel.edu/linbox?user=root");


	    stmt = conn.createStatement();

	    stmt.executeQuery("UPDATE servers SET in_use=0 WHERE ip = '" + InetAddress.getLocalHost().getHostName() + "'");
	}

	catch(ClassNotFoundException E) {
	        System.out.println("Error in ClassCreation of driver in Servlet");
	}

	catch(SQLException E) {
	      System.out.println(E.getMessage());
	}

	catch(Exception e) {
	                System.out.println("General Error in Servlet");
	    }





    }  // end of method computeResult
    
    public static synchronized void killProcess(){
        proc.destroy();
        //System.out.println("Smith Process has been killed");
    }

} // end of class SmithFormServlet
