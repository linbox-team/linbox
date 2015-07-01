<html>
<head><title>LinBox Computation Server</title>
<script type="text/javascript">
function disableRHS() { document.compInfo.rhs.disabled="disabled"; }
function enableGenerateMatrix() {
	with( document.compInfo ) {
		matrixType[0].disabled="";
		matrixType[1].disabled="";
		matrixType[2].disabled="";
		rows.disabled="";
		matrixURL.disabled="disabled";
		if( matrixType[1].checked ) identityClicked();
		if( matrixType[2].checked ) zeroClicked();
		if( matrixType[0].checked ) randomClicked();
	}
}
function disableGenerateMatrix() {
	with( document.compInfo ) {
		matrixType[0].disabled="disabled";
		matrixType[1].disabled="disabled";
		matrixType[2].disabled="disabled";
		rows.disabled="disabled";
		cols.disabled="disabled";
		sparsity.disabled="disabled";
		matrixURL.disabled="";
	}
}
function identityClicked() {
	with( document.compInfo ) {
		cols.disabled="disabled";
		sparsity.disabled="disabled";
	}
	rowsChanged();
}
function zeroClicked() {
	with( document.compInfo ) {
		cols.disabled="";
		sparsity.disabled="disabled";
	}
}
function randomClicked() {
	with( document.compInfo ) {
		cols.disabled="";
		sparsity.disabled="";
	}
}

function rowsChanged() {
	with( document.compInfo ) {
		rowParse = parseInt( rows.value );
		if( isNaN( rowParse ) || rowParse <= 0 ) {
			alert( "Please enter a positive integer row value." );
			return false;
		}
		else {
			rows.value = rowParse;
			if( matrixType[1].checked )
				cols.value = rows.value;
			return true;
		}
	}
}

function colsChanged() {
	with( document.compInfo ) {
		colParse = parseInt( cols.value );
		if( isNaN( colParse ) || colParse <= 0 ) {
			alert( "Please enter a positive integer column value.");
			return false;
		}
		else {
			cols.value = colParse;
			return true;
		}
	}
}

function sparsityChanged() {
	with( document.compInfo ) {
		sparParse = parseFloat( sparsity.value );
		if( isNaN( sparParse ) || sparParse < 0 || sparParse > 1 ) {
			alert( "Please enter a sparsity value between 0 and 1, inclusive.\n(0 means zero matrix, 1 means dense matrix.)" );
			return false;
		}
		else {
			sparsity.value = sparParse;
			return true;
		}
	}
}

function validate() {
	return( rowsChanged() && colsChanged() && sparsityChanged() );
}

</script>
</head>
<body>
<form action="submit-computation.jsp" method="post" name="compInfo">
<strong><font size=4>Choose a Matrix to Use<br></font></strong>
<table border="0" cellpadding="5"><tr><td>
<input type="radio" name="generateMatrix" value="yes" checked="checked" onClick="enableGenerateMatrix();">Generate a new matrix.
</td><td>
<input type="radio" name="matrixType" value="random" checked="checked" onClick="randomClicked();">Random matrix<br>
<input type="radio" name="matrixType" value="identity" onClick="identityClicked();">Identity matrix<br>
<input type="radio" name="matrixType" value="zero" onClick="zeroClicked();">Zero matrix
</td>
<td>
<table border="0">
<tr><td align="right">Rows:</td>
<td><input type="text" name="rows" value="100" size="5" onChange="rowsChanged();"></td></tr>
<tr><td align="right">Columns:</td>
<td><input type="text" name="cols" value="100" size="5" onChange="colsChanged();"></td></tr>
<tr><td align="right">Sparsity:</td>
<td><input type="text" name="sparsity" value=".10" size="5" onChange="sparsityChanged();"></td></tr>
</table>
</td>
</tr><tr><td colspan="3">
<input type="radio" name="generateMatrix" value="no" onClick="disableGenerateMatrix();">
Use matrix at the following URL:<input type="text" name="matrixURL" disabled="disabled" value="" size="40">
</td></tr></table>
<br>
<strong><font size=4>Computation to Perform<br></font></strong>
<input type="radio" name="type" value="4" checked="checked"
 onClick="disableRHS();">
Rank calculation
<br>

<input type="radio" name="type" value="5"
 onClick="disableRHS();">
Determinant computation
<br>

<input type="radio" name="type" value="6" disabled="disabled"
 onClick='document.compInfo.rhs.disabled=""'>
Solve Ax=b for x.  URL of vector b:
<input type="text" name="rhs" size="40" disabled="disabled">
<br>

<input type="radio" name="type" value="7" disabled="disabled"
 onClick="disableRHS();">
Solve Ax=<strong>0</strong> for x.
<br>

<input type="radio" name="type" value="8" disabled="disabled"
 onClick="disableRHS();">
Find the characteristic polynomial of a matrix.
<br>

<input type="radio" name="type" value="9" disabled="disabled"
 onClick="disableRHS();">
Find the minimum polynomial of a matrix.
<br>

<input type="radio" name="type" value="10" disabled="disabled"
 onClick="disableRHS();">
Compute the Smith normal form of a matrix.
<br>
<br>
<input type="radio" name="usePrime" value="no" checked="checked"
 onClick='document.compInfo.prime.disabled="disabled"'>
Compute over the integers.
<br>
<input type="radio" name="usePrime" value="yes" 
 onClick='document.compInfo.prime.disabled=""'>
Compute over a finite field with the following prime modulus:
<input type="text" name="prime" value="65521" size="10"i disabled="disabled">
<br>
<br>
<input type="button" value="Begin Computation"
 onClick="if(validate()) document.compInfo.submit();">
<input type="reset" value="Reset Form"
 onClick='disableRHS(); document.compInfo.prime.disabled="disabled"; enableGenerateMatrix(); randomClicked();'>
</form>
<hr>
LinBox has much more to offer than the limited capabilities of this server!  If you can't find what you're looking for here, or for further help, contact the LinBox team at <a href="mailto:linbox@yahoogroups.com">linbox@yahoogroups.com</a>.
<br>Server designed and maintained by <a href="mailto:roche@cis.udel.edu">Dan Roche</a>.
</body></html>
