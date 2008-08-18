Written by Matthew Fendt, Summer '07 to Winter '08

This directory contains the source files for the server side of the Linbox
web service.  It consists of:

- TransferAgent.java: This is the server side of the TransferAgent class that
 interacts directly with the Linbox functions

- TransferAgentSkeleton.java: This performs the necessary web service 
 manipulation that is necessary so that the client can retrieve the results

- xsd: This directory contains the Linbox functions that will be used by the
 server to compute the answer.  Interesting to note that this code is mainly
 in C++, with the necessary wrapper classes to make it able to interact with
 the rest of the server and client, both written in Java
