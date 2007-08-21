package samples.quickstart.service.adb;

import samples.quickstart.service.adb.xsd.*;

public class TransferAgentSkeleton implements TransferAgentSkeletonInterface
{
    static 
    {
	java.lang.System.load("/home/fendt/apache-tomcat-6.0.13/webapps/axis2/WEB-INF/services/liblinboxfunctions.so");
    }

    public RankResponse rank(Rank param0)
    {
	String r;
	
	r = linboxfunctions.rankFiles(param0.getMatrix());
	RankResponse res = new RankResponse();
	res.set_return(r);
	return res;
    }

    public DeterminantResponse determinant(Determinant param1)
    {
	String d;
	
	d = linboxfunctions.detFiles(param1.getMatrix());
	DeterminantResponse res = new DeterminantResponse();
	res.set_return(d);
	return res;
    }

    public ValenceResponse valence(Valence param2)
    {
	String v;

	v = linboxfunctions.valFiles(param2.getMatrix());
	ValenceResponse res = new ValenceResponse();
	res.set_return(v);
	return res;
    }

    public TraceResponse trace(Trace param3)
    {
	String tr;

	tr = linboxfunctions.traceFiles(param3.getMatrix());
	TraceResponse res = new TraceResponse();
	res.set_return(tr);
	return res;
    }
}
