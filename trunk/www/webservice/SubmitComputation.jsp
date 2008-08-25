<!-- Matthew Fendt Summer '07 to Winter'08
     Takes the user information from "Linalg Computation Client," and performs
     a Linbox computation
-->


<html>
<head>
<title> Submit Computation </title>
<link rel="stylesheet" href="linbox.css" title="linbox" type="text/css">
</head>

<h1> Submit Computation </h1>

<%@page import="java.util.*"%>
<%@page import="java.io.*"%>
<%@page import="java.net.*"%>
<%@page import="org.apache.commons.fileupload.servlet.ServletFileUpload"%>
<%@page import="org.apache.commons.fileupload.disk.DiskFileItemFactory"%>
<%@page import="org.apache.commons.fileupload.FileItem"%>
<%@page import="samples.quickstart.clients.TransferAgentDependentClient"%>
<%@page import="samples.quickstart.service.adb.TransferAgentTransferAgentHttpportStub"%>

<% 
	// File Information
	long fileSize = 0;
	String fileName = null;
	String contentType = null;

	// The client that will be called to perform the operaton
	TransferAgentDependentClient tac = new TransferAgentDependentClient();

	// The matrix inputted by the user as a string
	String matByHand = "";

	// The desired operation
	String op = null;

	// If the matrix is in a file, this is where the file will temporarily
	// reside
	File f = new File("/home/www/usa/linbox/apache-tomcat-6.0.13/webapps/ROOT/WEB-INF/classes/tempMat.txt");


	// The requests will always be multi-part because the form had a
	// "file upload" option
	if (ServletFileUpload.isMultipartContent(request))
		{
			// For dealing with files
			ServletFileUpload servletFileUpload = new
				ServletFileUpload(new DiskFileItemFactory());
			List fileItemsList = 
				servletFileUpload.parseRequest(request);

			Iterator it = fileItemsList.iterator();

			// Get the posted data
			while (it.hasNext()) {
			FileItem fileItem = (FileItem)it.next();
			if (fileItem.isFormField())
				{
				// The file item contains a simple name- value
				// pair of a form field

				String parameterName = fileItem.getFieldName();
				String parameterValue = fileItem.getString();

				if (parameterName.equals("matrixByHand"))
					matByHand = parameterValue;

				if (parameterName.equals("operation"))
					op = parameterValue;
				
				}

			else 
				{
				// The file item contains an uploaded file
				fileSize = fileItem.getSize();
				fileName = fileItem.getName();
				contentType = fileItem.getContentType();
				fileItem.write(f);
				}
			}
		}

	// This is the answer to the matrix computation
	String answer = null;


	// Compute the answer from the input
	 if (fileName.equals("") && 
		!(matByHand.equals("") || matByHand.equals(" "))) {
		answer = tac.request(op,matByHand, null, -1); %>
		<br><p>Computing answer from inputted matrix... <br>
	<% } %>


	<!-- Compute the answer from the file -->
	<% if ((matByHand.equals("") || matByHand.equals(" ")) && 
		!fileName.equals("")) {
		answer = tac.request(op, null, f, 1); %>
		<br><p>Computing answer from matrix file... <br>
	<% } %>
	

	<!-- Neither a file nor an input was given -->
	<% if ((matByHand.equals("") || matByHand.equals(" "))
		&& fileName.equals(""))
		answer = "Please input a matrix or upload a matrix file and try again";


	// Both a file and an input was given
	if (!(matByHand.equals("") || matByHand.equals(" ")) && 
		!(fileName.equals("")))
		answer = "Both a matrix input and a matrix file were detected.  Please only enter one or the other and try again.";


	// Delete the temporarily saved matrix file
	f.delete();
		
	%>

	<h2><u>Linbox's answer</u>:</h2>
	<%=op%> is: <%=answer%> <br>


 	<br><br>
	<h2><u>Uploaded File Information:</u>:</h2>

	<b>fileSize: </b><%=fileSize%> <br>
	<b>fileName: </b><%=fileName%> <br>
	<b>contentType: </b><%=contentType%> <br>

</body>
</html>
