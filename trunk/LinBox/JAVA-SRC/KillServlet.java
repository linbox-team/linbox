import javax.servlet.*;
import javax.servlet.http.*;
import java.io.*;


/**
 * @author Brian Gold
 * @version 1.0
 * 
 * Date: 10/15/01
 * Last Updated: 1/4/02
 * File: KillServlet.java
 * 
 * This file contains a doPost method which is called when the user hits
 * the button to halt the GAP calculations.  It determines which server was
 * running at the time, and calls the appropriate servlet (GapServlet or
 * SmithFormServlet) to stop the process.  It then notifies the user of
 * its success in stopping the calculations.
 */
public class KillServlet extends HttpServlet{


    /**
     * @param HttpServletRequest request
     * @param HttpServletResponse response
     * @return none
     *
     * This method determines which servlet was running GAP calculations
     * at the time the user decided to halt his calculations, and calls
     * the appropriate servlet to halt the process.
     */
    public void doPost(HttpServletRequest request,
		       HttpServletResponse response)
	throws ServletException, IOException{

        String stringServerToKill = request.getParameter("serverID");
        
	// Determine on which servlet the calculation to be halted was running
        if (stringServerToKill.equals("Simplicial") )
          GapServlet.killProcess();
        else if (stringServerToKill.equals("Smith") )
          SmithFormServlet.killProcess();

	// Display a message to the user; the calculations have been halted
        Message.DisplayKilledProcess(response);

    } // end of method doPost

}  // end of class KillServlet
