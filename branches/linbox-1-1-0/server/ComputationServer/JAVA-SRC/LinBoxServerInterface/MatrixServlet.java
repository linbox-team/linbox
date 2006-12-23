package LinBoxServerInterface;

import javax.servlet.*;
import javax.servlet.http.HttpServletRequest;
import java.io.*;

/** A servlet to handle requests for matrix files stored on the local machine.
  * Requests will only be serviced if the path is of the form
  * {dir}/{id number}.sms.  If using Apache Tomcat, be sure to add the
  * appropriate entries to WEB-INF/web.xml to make this servlet work.
  */
public class MatrixServlet extends GenericServlet {
	/** Service a request.
	  * If the request is valid, and a matrix with the given id is found,
	  * that matrix file will be returned to the requester.  Otherwise, an
	  * empty file is.
	  * @param req The request object; should be of type HttpServletRequest.
	  * @param resp The response object.
	  */
	public void service( ServletRequest req, ServletResponse resp )
		throws IOException
	{
		HttpServletRequest hr = (HttpServletRequest)req;
		String sp = hr.getServletPath();
		int id = -1;
		try {
		    id = Integer.parseInt(
			sp.substring( sp.lastIndexOf( '/' ) + 1,
			              sp.lastIndexOf( ".sms" ) )
				         );
		} catch( IndexOutOfBoundsException e ) { return; }
		catch( NumberFormatException e ) { return; }
		Computation comp = Computation.getByID(id);
		if( comp == null ) return;
		File matFile = comp.getMatrixFile();
		if( matFile == null ) return;
		BufferedReader in = 
			new BufferedReader( new FileReader( matFile ) );
		PrintWriter out = resp.getWriter();
		String line;
		while( (line = in.readLine()) != null )
			out.println( line );
		in.close();
		out.close();
	}
}
