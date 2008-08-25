<!-- Matthew Fendt Summer '07 to Winter'08
     Submits the user information to "SubmitComputation.jsp," which will
     perform a Linbox computation
-->


<html>
<head>
<title> Linalg Computing Client </title>
<link rel="stylesheet" href="linbox.css" title="linbox" type="text/css">
</head>
<body>

<h1> Linalg Computing Client </h1>

This is the ever-growing online computing client situated on linalg.org.
Check back frequently for more developments!


<!-- ******** The service ********* -->

<h3> Select a service </h3>

<form action="SubmitComputation.jsp" method="post" 
	enctype="multipart/form-data">
<p><select name="operation">
<option>Rank</option>
<option>Determinant</option>
<option>Valence</option>
<option>Trace</option>
</select></p>



<!-- ******** Upload the matrix ********* -->

<h3> Upload your matrix file </h3>

<input type="file" name="matrixFile" />



<!-- ******** Input matrix by hand ********* -->

<h3> Input your matrix file by hand </h3><br>

<textarea cols="60" rows="10" name="matrixByHand"> </textarea> <br><br>

<!-- ******** Submit matrix ********* -->

<input type="submit" value="Submit"/>


</form>
</body>
</html>
