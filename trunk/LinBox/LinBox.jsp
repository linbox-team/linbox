<%-- JavaServer Page to dynamically select a server based on servers in use. --%>

<%@ page import="mypackage.*,java.io.*" %>

<%! private serverSearch server = new serverSearch(); %>

<% if( server.checkServer().booleanValue() ){ %>

<html><head><title>LinBox Server</title></head><body bgcolor="white">
<a name="top">
</a><a href="#form"> to form </a>

<center><strong><font size="6">Sparse Integer Matrix Calculation with LinBox
</font></strong></center>
<p>
This page offers experimental access to the LinBox computational software, presently for access to rank, determinant, and minpoly calculations.


</p><p>
Please provide your matrix in a file using the following sparse matrix
format.  
</p><p>
<em>
Each line of the file contains a triple.
The first triple is number of rows, number of columns, and the letter ``M''.
The following triples signify non-zero entries and have the form
row-index, col-index, value.
They must be in increasing lexicographical order.
All values must be machine integers (ints).
The last triple is ``0 0 0'' used as an end marker.
For example:
</em>
</p><pre>32 3200 M 
1 1 345 
1 2 -23 
2 1 77 
2 200 31 
30 3101 11 
0 0 0 
</pre>
<em>
denotes a 32 by 3200 matrix with 5 non-zero entries, three of which
are clustered at the top left.
Optionally, the file may be gzipped.
</em>
<p>
<a name="form">
</a></p>



<hr>
<form action="http://<%= server.getServer() %>:8080/gap/SFLBServlet" method="post">

<a name="form"><strong><font size="4">URL of Matrix File<br>
<input type="text" name="URLmatrix" value="http://www.cis.udel.edu/~linbox/dumbexample.matrix" size="64"> <br> </font></strong>

<strong><font size="4">
Computation to Perform<br></font></strong>
<input type="radio" name="Smithchoice" value="rankLB" checked="checked">
Rank calculation performed on LinBox software
<br>

<input type="radio" name="Smithchoice" value="detLB">
Determinant Calculations performed on LinBox software
<br>

<input type="radio" name="Smithchoice" value="minpolyLB">
MinPoly Calculations performed on LinBox software
<br>

</a><p>

<a name="form"><strong><font size="4"><i>Non-Optional</i> Prime Modulus m.
</font></strong> 
<input type="text" name="modulusm" value="101" size="10">
<br>

<!-- Removed functionality 

</a></p><p>
<a name="form">For details on the following three options see the </a><a href="http://www.cis.udel.edu/%7Edumas/Homology/manual/">manual</a>.
<br>
<strong><font size="4">Algorithm:
</font></strong>
<input type="radio" name="algorithm" value="elim" checked="checked"> 
Via elimination 

<input type="radio" name="algorithm" value="BB"> 
Via blackboxes 

</p><p>

<strong><font size="4">Commentary:
<!--The code is equipped to offer some commentary as the computation proceeds.  Please choose how much to see.

</font></strong>
<input type="radio" name="comments" value="none" checked="checked"> 
No Comments 

<input type="radio" name="comments" value="some"> 
Some Comments 

<input type="radio" name="comments" value="verbose"> 
Verbose 

</p><p>
<strong><font size="4">Timing info: 
</font></strong>
<input type="radio" name="timing" value="none" checked="checked"> 
No timings 

<input type="radio" name="timing" value="high level"> 
Major steps timed

<input type="radio" name="timing" value="detail"> 
Detailed timing

--!>

</p><p>
<strong><font size="4">Results: 
</font></strong>
The computation may take a while.  After you submit this form,
options for monitoring the computation and seeing the results will
be offered.  You may elect to also receive email:
<br>
<input type="radio" name="results_option" value="wait" checked="checked"> 
No email.
<br>
<input type="radio" name="results_option" value="shortemail"> 
Send email containing web address of results.
<br>
<input type="radio" name="results_option" value="longemail"> 
Send the results in email (may be long).
<br>
E-mail address (user@host.ext): <font size="24"> </font>
<input type="text" name="email" size="32">
<br>

<input type="submit" value="Begin Computation">

<input type="reset" value="Reset Form">
</form>

</p><hr>



<!--
<form action="TopServlet" method="get">
To view the current memory load (the 'top' command) of the computation server, 
select 
<a href="http://algol.cis.udel.edu:8080/gap/TopServlet" target="main">load</a>. 
</form>
-->

<p>
For other online computations using the Simplicial Homology or Smith Form GAP Package, go to the <a href="http://linalg.org">Linalg.org</a> Main Page and select the choice you would like.

</p><p>

Please send any questions or comments about the server to <a href="mailto:%20schrag@mail.eecis.udel.edu">Eric Schrag</a> and any questions or comments about the GAP Homology Package (including Smith Normal form) to <a href="mailto:%20saunders@cis.udel.edu">David Saunders</a>.  We would greatly appreciate any feedback you might have.

<BR><BR>


</p><p>
<tiny>Version 1.10 </tiny>
<br><tiny>Revised 27Jan2003</tiny>
<br>

<a href="#top">to top</a>
<a href="#form">to form</a></p></body></html>

<%  }  else { %>

The LinBox team apologizes, but all servers in this LinBox cluster are either currently in use or down.
<BR>
Either someone else is running calculations already on all available servers, or they are down for maintenance.
<BR>
Please come back soon.

<% } %>
