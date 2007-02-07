<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
<head>
 <title>/trunk/linbox/examples/bigmat.C - linalg - Trac</title><link rel="start" href="/projects/linalg/wiki" /><link rel="search" href="/projects/linalg/search" /><link rel="help" href="/projects/linalg/wiki/TracGuide" /><link rel="stylesheet" href="/projects/linalg/chrome/common/css/trac.css" type="text/css" /><link rel="stylesheet" href="/projects/linalg/chrome/common/css/code.css" type="text/css" /><link rel="stylesheet" href="/projects/linalg/chrome/common/css/browser.css" type="text/css" /><link rel="icon" href="/projects/linalg/chrome/common/trac.ico" type="image/x-icon" /><link rel="shortcut icon" href="/projects/linalg/chrome/common/trac.ico" type="image/x-icon" /><link rel="up" href="/projects/linalg/browser/trunk/linbox/examples?rev=2574" title="Parent directory" /><link rel="alternate" href="/projects/linalg/browser/trunk/linbox/examples/bigmat.C?rev=2574&amp;format=txt" title="Plain Text" type="text/plain" /><link rel="alternate" href="/projects/linalg/browser/trunk/linbox/examples/bigmat.C?rev=2574&amp;format=raw" title="Original Format" type="text/x-c++src" /><style type="text/css">
</style>
 <script type="text/javascript" src="/projects/linalg/chrome/common/js/trac.js"></script>
</head>
<body>


<div id="banner">

<div id="header"><a id="logo" href="http://trac.edgewall.com/"><img src="/projects/linalg/chrome/common/trac_banner.png" width="236" height="73" alt="Trac" /></a><hr /></div>

<form id="search" action="/projects/linalg/search" method="get">
 <div>
  <label for="proj-search">Search:</label>
  <input type="text" id="proj-search" name="q" size="10" accesskey="f" value="" />
  <input type="submit" value="Search" />
  <input type="hidden" name="wiki" value="on" />
  <input type="hidden" name="changeset" value="on" />
  <input type="hidden" name="ticket" value="on" />
 </div>
</form>



<div id="metanav" class="nav"><ul><li class="first"><a href="/projects/linalg/login">Login</a></li><li><a href="/projects/linalg/settings">Settings</a></li><li><a href="/projects/linalg/wiki/TracGuide" accesskey="6">Help/Guide</a></li><li class="last"><a href="/projects/linalg/about">About Trac</a></li></ul></div>
</div>

<div id="mainnav" class="nav"><ul><li class="first"><a href="/projects/linalg/wiki" accesskey="1">Wiki</a></li><li><a href="/projects/linalg/timeline" accesskey="2">Timeline</a></li><li><a href="/projects/linalg/roadmap" accesskey="3">Roadmap</a></li><li class="active"><a href="/projects/linalg/browser">Browse Source</a></li><li><a href="/projects/linalg/report">View Tickets</a></li><li><a href="/projects/linalg/newticket" accesskey="7">New Ticket</a></li><li class="last"><a href="/projects/linalg/search" accesskey="4">Search</a></li></ul></div>
<div id="main">




<div id="ctxtnav" class="nav">
 <ul>
  <li class="last"><a href="/projects/linalg/log/trunk/linbox/examples/bigmat.C?rev=2574">Revision Log</a></li>
 </ul>
</div>

<div id="content" class="browser">
 <h1><a class="first" title="Go to root directory" href="/projects/linalg/browser/?rev=2574">root</a><span class="sep">/</span><a title="View trunk" href="/projects/linalg/browser/trunk?rev=2574">trunk</a><span class="sep">/</span><a title="View linbox" href="/projects/linalg/browser/trunk/linbox?rev=2574">linbox</a><span class="sep">/</span><a title="View examples" href="/projects/linalg/browser/trunk/linbox/examples?rev=2574">examples</a><span class="sep">/</span><a title="View bigmat.C" href="/projects/linalg/browser/trunk/linbox/examples/bigmat.C?rev=2574">bigmat.C</a></h1>

 <div id="jumprev">
  <form action="" method="get"><div>
   <label for="rev">View revision:</label>
   <input type="text" id="rev" name="rev" value="2574" size="4" />
  </div></form>
 </div>

 
  <table id="info" summary="Revision info"><tr>
    <th scope="row">
     Revision <a href="/projects/linalg/changeset/2572">2572</a>
     (checked in by youse, 2 weeks ago)
    </th>
    <td class="message"><p>
Update of parallel CRA code using MPI. <br />
</p>
</td>
   </tr></tr>
  </table>
  <div id="preview"><table class="code"><thead><tr><th class="lineno">Line</th><th class="content">&nbsp;</th></tr></thead><tbody><tr><th id="L1"><a href="#L1">1</a></th>
<td>#include &lt;iostream&gt;</td>
</tr><tr><th id="L2"><a href="#L2">2</a></th>
<td>using namespace std;</td>
</tr><tr><th id="L3"><a href="#L3">3</a></th>
<td></td>
</tr><tr><th id="L4"><a href="#L4">4</a></th>
<td></td>
</tr><tr><th id="L5"><a href="#L5">5</a></th>
<td>// Makes an n by n matrix whose determinant is considerably less than the </td>
</tr><tr><th id="L6"><a href="#L6">6</a></th>
<td>// Hadamard bound.</td>
</tr><tr><th id="L7"><a href="#L7">7</a></th>
<td></td>
</tr><tr><th id="L8"><a href="#L8">8</a></th>
<td>int main(int argc, char* argv[])</td>
</tr><tr><th id="L9"><a href="#L9">9</a></th>
<td>{</td>
</tr><tr><th id="L10"><a href="#L10">10</a></th>
<td>&nbsp; &nbsp; if (argc != 2 ) {</td>
</tr><tr><th id="L11"><a href="#L11">11</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; cerr &lt;&lt; &#34;Usage: bigmat &lt;n&gt;, where &lt;n&gt; is the size you like.&#34; &lt;&lt; endl;</td>
</tr><tr><th id="L12"><a href="#L12">12</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; return -1;</td>
</tr><tr><th id="L13"><a href="#L13">13</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; }</td>
</tr><tr><th id="L14"><a href="#L14">14</a></th>
<td></td>
</tr><tr><th id="L15"><a href="#L15">15</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; int n = atoi(argv[1]);</td>
</tr><tr><th id="L16"><a href="#L16">16</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; n &lt;&lt; &#34; &#34; &lt;&lt; n &lt;&lt; &#34; M&#34; &lt;&lt; endl;</td>
</tr><tr><th id="L17"><a href="#L17">17</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; &#34;1 1 3&#34; &lt;&lt; endl;</td>
</tr><tr><th id="L18"><a href="#L18">18</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; &#34;1 &#34; &lt;&lt; n &lt;&lt; &#34; 2&#34; &lt;&lt; endl;</td>
</tr><tr><th id="L19"><a href="#L19">19</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; for (int i = 2; i &lt;=n; ++i) </td>
</tr><tr><th id="L20"><a href="#L20">20</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; {</td>
</tr><tr><th id="L21"><a href="#L21">21</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; i &lt;&lt; &#34; &#34; &lt;&lt; i-1 &lt;&lt; &#34; &#34; &lt;&lt; 2 &lt;&lt; endl;</td>
</tr><tr><th id="L22"><a href="#L22">22</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; i &lt;&lt; &#34; &#34; &lt;&lt; i &lt;&lt; &#34; &#34; &lt;&lt; 3 &lt;&lt; endl;</td>
</tr><tr><th id="L23"><a href="#L23">23</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; }</td>
</tr><tr><th id="L24"><a href="#L24">24</a></th>
<td>&nbsp; &nbsp; &nbsp; &nbsp; cout &lt;&lt; &#34;0 0 0&#34; &lt;&lt; endl;</td>
</tr><tr><th id="L25"><a href="#L25">25</a></th>
<td>}</td>
</tr></tbody></table>
  </div>

 <div id="help">
  <strong>Note:</strong> See <a href="/projects/linalg/wiki/TracBrowser">TracBrowser</a> for help on using the browser.
 </div>

</div>
<script type="text/javascript">searchHighlight()</script>
<div id="altlinks"><h3>Download in other formats:</h3><ul><li class="first"><a href="/projects/linalg/browser/trunk/linbox/examples/bigmat.C?rev=2574&amp;format=txt">Plain Text</a></li><li class="last"><a href="/projects/linalg/browser/trunk/linbox/examples/bigmat.C?rev=2574&amp;format=raw">Original Format</a></li></ul></div>

</div>

<div id="footer">
 <hr />
 <a id="tracpowered" href="http://trac.edgewall.com/"><img src="/projects/linalg/chrome/common/trac_logo_mini.png" height="30" width="107"
   alt="Trac Powered"/></a>
 <p class="left">
  Powered by <a href="/projects/linalg/about"><strong>Trac 0.9.6</strong></a><br />
  By <a href="http://www.edgewall.com/">Edgewall Software</a>.
 </p>
 <p class="right">
  Visit the Trac open source project at<br /><a href="http://trac.edgewall.com/">http://trac.edgewall.com/</a>
 </p>
</div>



 </body>
</html>

