<html>
<head><title>LinBox Computation Results</title>
<%@ page import="LinBoxServerInterface.Computation" %>
<%!
	private static String br = "<br>";
	protected String addHTMLBreaks( String s ) {
		StringBuffer sb = new StringBuffer(s);
		for( int i = 0; i < sb.length(); ++i ) {
			if( sb.charAt(i) == '\n' ) {
				sb.insert(i,br);
				i += 4;
			}
		}
		return sb.toString();
	}
%>
</head>
<body>
<strong>Server responded with:</strong><br>
<font face="courier">
<%
    Computation comp = Computation.getByID(Integer.parseInt(request.getParameter("id")));
    if( comp == null ) out.println( "Sorry, your computation has been lost." );
    else {
	comp.updateAccessed();
	Thread t = Thread.currentThread();
	comp.addWaiting( t );
    	while( comp.getResult() == null ) {
		try{ t.sleep( 10000 ); }
		catch( InterruptedException e ) {}
	}
	comp.removeWaiting(t);
	out.println( addHTMLBreaks( comp.getResult() ) );
	comp.tryDestroy();
    }
%>
</font>
<br>
<br>
<a href="computation.jsp">Back to the computation form</a><br>
<a href="http://www.linalg.org">Back to LinBox home page</a>
</body></html>
