import java.net.*;
import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;

/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 7/25/01
 * Last Updated: 2/1/02
 * File: Message.java
 * 
 * This file contains only static methods which display dynamic HTML
 * pages to the user
 */
class Message{

    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page which notifies the
     * user that his GAP calculations have been successfully halted.
     * This is called by KillServlet
     */
    public static void DisplayKilledProcess(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Homology Server - Calculations Killed\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("You have successfully halted the calculations of your ");
        sbBuf.append("GAP submission.<br><br>");
        sbBuf.append("<form action=TopServlet method=get>");
        sbBuf.append("To view the current memory load (the 'top' command) ");
        sbBuf.append("of the Homology Server, click&nbsp;");
        sbBuf.append("<A HREF=\"TopServlet\" TARGET=\"main\">here</A>.</form>");
        sbBuf.append("<br><br>");
	sbBuf.append("<a href=\"http://linbox.pc.cis.udel.edu/LinBox.jsp\">Back to the LinBox Server</a>");
	sbBuf.append("<br><br><br>Thank you for using the LinBox Team Servers.");
	sbBuf.append("<br><br>top of <a href=\"#top\">page</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }	// end of method DisplayKilledProcess



    /**
     * @param HttpServletResponse response
     * @param String EmailAddress
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * tells the user that an email will be sent when the calculations
     * are completed.  This is used by the Smith Form Server
     */
    public static void DisplayEmailPage(HttpServletResponse response,
					  String EmailAddress, int randomInt)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Smith Form Server - Calculations Currently Running\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=\"4\"><strong>");
	sbBuf.append("Your request is currently being calculated.  ");
	sbBuf.append("</strong></size></strong><br><br>");
        if (EmailAddress.equals("") == false){
   	    sbBuf.append("When it is finished, an e-mail will be sent to ");
	    sbBuf.append(EmailAddress);
	    sbBuf.append(" informing you how to retrieve the results.  ");
        }


	sbBuf.append("<TABLE CELLSPACING=\"5\"><TR><TH></strong></TH>");
	sbBuf.append("<TH></strong>Stop Current Calculations</TH><TH>");
	sbBuf.append("</TH><TH></strong>Current memory load of Server**</TH>");
	sbBuf.append("</TR><TR><TH>");
	sbBuf.append("<form action=SFLBServlet method=get>");
	sbBuf.append("<INPUT TYPE=hidden name=\"randomint\" value=");
	sbBuf.append(randomInt);
	sbBuf.append(">");
	sbBuf.append("<INPUT TYPE=submit value=\"View Results\">");
	sbBuf.append("</form></TH><TH>");
	sbBuf.append("<form action=KillServlet method=post>");
	sbBuf.append("<INPUT TYPE=hidden name=\"serverID\" value=Smith>");
	sbBuf.append("<INPUT TYPE=submit value=\"Kill Calculations\">");
	sbBuf.append("</form></TH><TH></TH>");
	sbBuf.append("<TH><form action=TopServlet method=get>");
	sbBuf.append("<A HREF=\"TopServlet\" TARGET=\"main\">Memory Load</A></form>");
	sbBuf.append("</TH></TR></strong><TR><font size=\"2\">");
	sbBuf.append("<TH></TH><TH></TH><TH></TH><TH></strong>");
	sbBuf.append("**Opens new browser");
	sbBuf.append(" window</TH></size></font></TR></TABLE>");
	sbBuf.append("<br> If you want to leave and come back for your results later,<br>");
	sbBuf.append("bookmark <A HREF=\"http://" + InetAddress.getLocalHost().getHostName() + ":8080/gap/Result.html\">this link</A> and enter the key number below at that page.<br>");
	sbBuf.append("Your key is " + randomInt);
	sbBuf.append("<br><br><br>");
	sbBuf.append("<a href=\"http://linbox.pc.cis.udel.edu:8080/gap/LinBox.jsp\">Back to the LinBox");
	sbBuf.append(" Server</a>");
	sbBuf.append("<br><br><br>");
	sbBuf.append("<br><br><br>Thank you for using the LinBox Team Servers.");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }	// end of method DisplayEmailPage


    /**
     * @param HttpServletResponse response
     * @param String result
     * @param String comments
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * contains the results of the Gap computations, found in the result 
     * string.  If the user has selected to keep the comments and results
     * together, the comments will also be found in this string.  If the
     * user has selected for the comments to be separate, the comments
     * will be in the comments string.  If the user does not want the
     * comments, the comments string will be "" (NOTE: Not the same as NULL)
     * If the results are from the Homology server, boolHomology will be
     * true, otherwise, if the results are from the SmithForm server,
     * boolHomology will be false
     * This method is called by ResultServlet, SmithFormServlet, GapServlet.
     */
    public static void DisplayGapResult(HttpServletResponse response,
					String result, String comments)
	throws IOException{

	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Server Result\n");

	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<a name=\"top\">");
	sbBuf.append("<font size=36> Results: </font><br>");
	sbBuf.append(result);
	sbBuf.append("<br><br><br>");
	sbBuf.append("<a href=\"http://linalg.org/servers.html\">Back to the Linalg Computational Page</a>");
	sbBuf.append("<br><br><br>");
	sbBuf.append("<br><br><br>Thank you for using the LinBox Team Servers.");
	sbBuf.append("<br><br> <a href=\"#top\">to top</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }	// end of method DisplayGapResult


    /**
     * @param HttpServletResponse response
     * @param String commands
     * @param String EmailAddress
     * @param int randomInt
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * tells the user that the Gap calculations are running, and, if
     * "EmailAddress" is not "", tells the user that notification will be
     * emailed when the calculations are complete.  The string "commands"
     * contains all the commands submitted by the user, which are displayed
     * on the generated page.  It also contains a
     * "Submit" button which, when pressed, loads the doGet method of
     * GapServlet, in order to display the current results of the calculations.
     * calculations are still running.  The integer randomInt is sent through
     * a hidden HTML parameter to the doGet method of GapServlet.
     * This method is called from GapServlet
     */
    public static void DisplayRollingPage(HttpServletResponse response,
					String commands, String EmailAddress, 
					int randomInt)
	throws IOException{

	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Homology Server - Calculations Currently Running\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	
	sbBuf.append("<TABLE CELLSPACING=\"5\"><TR><TH></strong></TH>");
	sbBuf.append("<TH></strong>Stop Current Calculations</TH><TH></TH>");
	sbBuf.append("<TH></strong>Current memory load of Server**</TH>");
	sbBuf.append("</TR><TR><TH>");
	sbBuf.append("<form action=GapServlet method=get>");
	sbBuf.append("<INPUT TYPE=hidden name=\"randomint\" value=");
	sbBuf.append(randomInt);
	sbBuf.append(">");
	sbBuf.append("<INPUT TYPE=submit value=\"View Results\">");
	sbBuf.append("</form></TH><TH>");
	sbBuf.append("<form action=KillServlet method=post>");
	sbBuf.append("<INPUT TYPE=hidden name=\"serverID\" value=Simplicial>");
	sbBuf.append("<INPUT TYPE=submit value=\"Kill Calculations\">");
	sbBuf.append("</form></TH><TH></TH>");
	sbBuf.append("<TH><form action=TopServlet method=get>");
	sbBuf.append("<A HREF=\"TopServlet\" TARGET=\"main\">Memory Load</A></form>");
	sbBuf.append("</TH></TR></strong><TR><font size=\"2\">");
	sbBuf.append("<TH></TH><TH></TH><TH></TH><TH></strong>");
	sbBuf.append("**Opens new browser");
	sbBuf.append(" window</TH></size></font></TR></TABLE>");

	sbBuf.append("<hr>The following commands are currently being ");
	sbBuf.append("calculated.");

	if (EmailAddress.equals("")==false){
	    sbBuf.append("  When they are finished, an ");
	    sbBuf.append("e-mail will be sent to ");
	    sbBuf.append(EmailAddress);
	    sbBuf.append(" informing you how to retrieve the results.");
	}

	sbBuf.append("<br><br>");
	sbBuf.append("<font size=\"4\"> Commands: </font><br>");
	int index;
	while ( (index = commands.indexOf(';') ) != (-1) ){
	    sbBuf.append(commands.substring(0, index+1) );
	    sbBuf.append("<br>");
	    commands = commands.substring(index+1);
	}
	sbBuf.append(commands);
	sbBuf.append("<br><hr>");
	sbBuf.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
	sbBuf.append("<br><br><br>");
	sbBuf.append("<a href=Result.html>Retrieve the Results</a>");
	sbBuf.append("<br><br><br>Thank you for using the Homology Server.");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }	// end of method DisplayRollingPage


    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message page to the user, telling the user that
     * the key entered into Result.html's text field is invalid.
     * This method is called by ResultServlet
     */
    public static void DisplayResultError(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Invalid Key\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Please enter a valid key");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("<a href=Result.html>Enter Your Key Again</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of message DisplayResultError


    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message page to the user, telling the user that
     * the e-mail address that was entered is incorrectly formatted or
     * that it is missing when it should have been entered.
     * This method is called by GapServlet.
     */
    public static void DisplayEmailError(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Missing Field(s)\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Please enter a valid e-mail address");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of method DispalyEmailError


    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message to the user telling the user that 
     * neither a URL of gap commands was entered, nor was a list of gap
     * commands in the supplied text box entered into the form Gap.html
     * This method is called by GapServlet.
     */
    public static void DisplayNoneError(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Missing Field(s)\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Please fill in one of the two ");
	sbBuf.append("fields, either a URL or a string of commands");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of method DispalyNoneError


    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message to the user telling the user that 
     * both a URL of gap commands was entered, AND a list of gap
     * commands was supplied in the text box in the form Gap.html.  As
     * a result, the Server does not know which calculations to perform,
     * so it performs neither, and generates this error message.
     * This method is called by GapServlet.
     */
    public static void DisplayBothError(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Extra Field(s)\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Please fill in <i>ONLY ONE</i> ");
	sbBuf.append("of the two ");
	sbBuf.append("fields, either a URL or a string of commands, ");
	sbBuf.append("but not both");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of method DisplayBothError



    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message to the user telling the user that 
     * the URL entered for the gap commands was invalid.  The same applies
     * for the URL entered in the optional matrix-textfield in Gap.html
     * This method is called by SmithFormServlet and GapServlet.
     */
    public static void DisplayURLErrorMessage(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Page Not Found\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<STRONG><font size=36>ERROR: Invalid URL\n");
	sbBuf.append("</font></STRONG><br><br>");
	sbBuf.append("<a href=SimplicialHomologyForm.html>Back to the Homology Server</a>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of method DisplayURLErrorMessage


    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message page to the user, telling the user that
     * the number entered in the optional modulus box is invaid.
     * This method is called by SmithFormServlet.
     */
    public static void DisplayInvalidNumberMessage(HttpServletResponse response)
	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Invalid Number\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Please enter a valid modulus integer.");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of message DisplayInvalidNumberMessage



    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message page to the user, telling the user that
     * the number entered in the optional modulus box is invaid.
     * This method is called by SmithFormServlet.
     */
    public static void DisplayInvalidCommand(HttpServletResponse response)

	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Invalid Command\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: Your code contained an invalid Exec(*) command. Please remove it.\n");
	sbBuf.append("<font size=36> If it is truly necessary, please email the administrators of the site.\n");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of message DisplayInvalidCommand






    /**
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * displays an error-message page to the user, telling the user that
     * the number entered in the optional modulus box is invaid.
     * This method is called by SmithFormServlet.
     */
    public static void DisplayModulusMissing(HttpServletResponse response)

	throws IOException{
	
	PrintWriter pwOutput;
	response.setContentType("text/html");
	pwOutput = response.getWriter();
	StringBuffer sbBuf = new StringBuffer();
	  
	sbBuf.append("<HTML><HEAD><TITLE>\n");
	sbBuf.append("Error - Missing Data\n");
	sbBuf.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBuf.append("<font size=36> ERROR: You have not entered a prime modulus, which is necessary for LinBox.\n");
	sbBuf.append("<font size=18> If you want integer calculations, try the GAP SmithForm page.\n");
	sbBuf.append("</size> </font> <br><br><br>");
	sbBuf.append("</BODY></HTML>");
	pwOutput.println(sbBuf.toString() );
	pwOutput.flush();
	pwOutput.close();
    }  // end of message DisplayModulusMissing




}  // end of Message class

