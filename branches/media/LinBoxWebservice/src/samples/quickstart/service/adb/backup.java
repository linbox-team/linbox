package samples.quickstart.service.adb;

/*
import samples.quickstart.service.adb.xsd.Determinant;
import samples.quickstart.service.adb.xsd.DeterminantResponse;
import samples.quickstart.service.adb.xsd.Rank;
import samples.quickstart.service.adb.xsd.RankResponse;
import samples.quickstart.service.adb.xsd.Trace;
import samples.quickstart.service.adb.xsd.TraceResponse;
import samples.quickstart.service.adb.xsd.Valence;
import samples.quickstart.service.adb.xsd.ValenceResponse;
*/

import samples.quickstart.service.adb.xsd.*;

//import linboxSharedLib.*;

public class TransferAgentSkeleton implements TransferAgentSkeletonInterface
{
    public RankResponse rank(Rank param0)
    {
	int r;

	//	try{java.lang.System.loadLibrary("linboxfunctions");}
	try{java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");}
	catch (UnsatisfiedLinkError e) 
	    {
		r = -13;
		RankResponse res = new RankResponse();
		res.set_return(r);
		return res;
	    }
	
	r = linboxfunctions.rankFiles(param0.getMatrix());
	RankResponse res = new RankResponse();
	res.set_return(r);
	return res;
    }

    public DeterminantResponse determinant(Determinant param1)
    {
	long d;
	try{java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");}
	catch (UnsatisfiedLinkError e) 
	    {
		d = -22;
		DeterminantResponse res = new DeterminantResponse();
		res.set_return(d);
		return res;
	    }
	
	d = linboxfunctions.detFiles(param1.getMatrix());
	DeterminantResponse res = new DeterminantResponse();
	res.set_return(d);
	return res;
    }

    public ValenceResponse valence(Valence param2)
    {
	int v;
	try{java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");}
	catch (UnsatisfiedLinkError e) 
	    {
		v = -33;
		ValenceResponse res = new ValenceResponse();
		res.set_return(v);
		return res;
	    }
	
	v = linboxfunctions.valFiles(param2.getMatrix());
	ValenceResponse res = new ValenceResponse();
	res.set_return(v);
	return res;
    }

    public TraceResponse trace(Trace param3)
    {
	int tr;

	try{java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");}
	catch (UnsatisfiedLinkError e) 
	    {
		tr = -44;
		TraceResponse res = new TraceResponse();
		res.set_return(tr);
		return res;
	    }
	
	tr = linboxfunctions.traceFiles(param3.getMatrix());
	TraceResponse res = new TraceResponse();
	res.set_return(tr);
	return res;
    }
}
