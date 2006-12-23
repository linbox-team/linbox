<html>
<head><title>Submitting computation...</title>
<%@ page import="LinBoxServerInterface.*" %>
<%!
    static private ComputingServerConnection csc =
    	ComputingServerConnection.getInstance();
%>
</head>
<body>
<h3>Submitting your computation...</h3>
<%
    Computation comp = new Computation();
    comp.setId();
    if( request.getParameter("usePrime").equals( "yes" ) )
	comp.setPrime( Long.parseLong( request.getParameter("prime") ) );
    comp.setType( Integer.parseInt( request.getParameter("type") ) );
    if( request.getParameter("generateMatrix").equals( "yes" ) ) {
	if( request.getParameter("matrixType").equals( "random" ) )
	    comp.setMatrixFile( MatrixGenerator.newRandom(
	    	comp.getId(),
		Integer.parseInt( request.getParameter("rows") ),
		Integer.parseInt( request.getParameter("cols") ),
		Float.parseFloat( request.getParameter("sparsity") ) ) );
	else if( request.getParameter("matrixType").equals( "identity" ) )
	    comp.setMatrixFile( MatrixGenerator.newIdentity(
	    	comp.getId(),
		Integer.parseInt( request.getParameter("rows") ) ) );
	else if( request.getParameter("matrixType").equals( "zero" ) )
	    comp.setMatrixFile( MatrixGenerator.newZero(
	    	comp.getId(),
		Integer.parseInt( request.getParameter("rows") ),
		Integer.parseInt( request.getParameter("cols") ) ) );
    }
    else comp.setMatrixURL( request.getParameter( "matrixURL" ) );
    if( comp.getType() == ComputationCodes.SOLVE )
    	comp.setVectorURL( request.getParameter( "vectorURL" ) );
    comp.saveArchive();
    csc.sendNew( comp );
    out.println( "<script type=\"text/javascript\">" );
    out.println( "window.location.replace(\"get-results.jsp?id=" + comp.getId() + "\");" );
    out.println( "</script>" );
%>
</body></html>
