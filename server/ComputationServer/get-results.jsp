<html>
<head><title>Retrieving your computation results...</title>
</head>
<body>
<h2>Please wait while your computation is being performed...</h2><br>
<form name="dumbForm" method="GET" action="results.jsp">
<% out.println( "<input type=\"hidden\" name=\"id\" value=\"" + request.getParameter("id") + "\">" ); %>
</form>
<h4>If your computation is taking a while, you can leave this page and return later to check on the results.<br>Just go to <% out.print( "<a href=\"http://linbox.pc.cis.udel.edu:8080/ComputationServer/get-results.jsp?id=" + request.getParameter("id") + "\">" + "http://linbox.pc.cis.udel.edu:8080/ComputationServer/get-results.jsp?id=" + request.getParameter("id") + "</a>"); %>
<script type="text/javascript">
	document.dumbForm.submit();
</script>
</body></html>
