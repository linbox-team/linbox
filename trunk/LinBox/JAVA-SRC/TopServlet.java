import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 11/15/01
 * Last Updated: 1/4/02
 * File: TopServlet.java
 * 
 * This file contains a doGet method which is called when the user clicks on
 * a link to see the memory load of the server (the unix "top" command).  
 * The results of this top command are displayed in a new browser window
 * for the user and are automatically updated every 60 seconds, or whenever the
 * user hits the "submit" button
 
 */
public class TopServlet extends HttpServlet{


    /**
     * @param HttpServletRequest request
     * @param HttpServletResponse response
     * @return none
     *
     * This method generates and displays an HTML page to the user which
     * contains the memory load of the server (the unix "top" command)
     * at the current time.
     */
    public void doGet(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{

	// Run the "top" command in unix once
        Process proc = Runtime.getRuntime().exec("top -n");

        String s;
	StringBuffer sb = new StringBuffer();

	// Set up a buffer to read the output returned by server
	BufferedReader br = new BufferedReader (new InputStreamReader
	    (proc.getInputStream() ) );
	
        int tempindex = 0;

	// Read all of the results outputted by the top command
	while ((s = br.readLine()) != null){
	    
	    tempindex = 0;
	    while( (tempindex=s.indexOf(" "))!=(-1) ){
		s = s.substring(0, tempindex) +
		    "&nbsp;" + s.substring(tempindex+1);
	    }

	    sb.append(s).append("<br>");
	}

        // Set up the printing mechanism to display Dynamic HTML
	PrintWriter pwOutput;
        response.setContentType("text/html");
	pwOutput = response.getWriter();
	    
        StringBuffer sbBeg = new StringBuffer();

	// Set up the dynamic HTML
	sbBeg.append("<HTML><HEAD><TITLE>\n");
	sbBeg.append("Memory Load of the Server\n");
        sbBeg.append("</TITLE></HEAD><BODY bgcolor=white>\n");
	sbBeg.append("<meta http-equiv=\"refresh\" content=\"60;>");
	sbBeg.append("<meta name=\"keywords\" content=\"automatic redirection\">");
        sbBeg.append("<br>");
        sbBeg.append("<font face=\"Courier\"> ");
        sbBeg.append(sb);
        sbBeg.append("<br><br><hr>");
	sbBeg.append("To update the GAP/Homology output manually ");
	sbBeg.append("click the SUBMIT button.  Otherwise, it will ");
	sbBeg.append("be automatically updated every 60 seconds.");
	sbBeg.append("<br><br>");
	sbBeg.append("<form action=TopServlet method=get><br>");
	sbBeg.append("<INPUT TYPE=submit value=\"SUBMIT\">");
        sbBeg.append("</form></font>");
        sbBeg.append("</BODY></HTML>");

        // Display the "top" results to the user
	pwOutput.println(sbBeg.toString() );
	pwOutput.flush();
	pwOutput.close();

    } // end of method doGet

}  // end of class TopServlet
