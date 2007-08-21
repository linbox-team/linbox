// Matthew Fendt, Summer '07 Research
// Client that accesses the LinBox Webservice.


package samples.quickstart.clients;

import samples.quickstart.service.adb.TransferAgentTransferAgentHttpportStub;

public class TransferAgentClient {

    public static void main(String[] args)
    {
	try
	    {
		TransferAgentTransferAgentHttpportStub stub = new TransferAgentTransferAgentHttpportStub("http://hmrg.pc.cis.udel.edu:2007/axis2/services/TransferAgent");

		request(stub, args);
	    }

	catch(Exception e)
	    {
		e.printStackTrace();
	    }
    }


    public static void request(TransferAgentTransferAgentHttpportStub stub, String[] args)
    {
	try{
	    if (args[0].equalsIgnoreCase("rank")){
		TransferAgentTransferAgentHttpportStub.Rank req = 
		    new TransferAgentTransferAgentHttpportStub.Rank();
		req.setMatrix(args[1]);
		
		TransferAgentTransferAgentHttpportStub.RankResponse res = 
		    stub.rank(req);
		System.out.println("Rank is:    " + res.get_return());
	    }
	    
	    else if (args[0].equalsIgnoreCase("determinant")){
		TransferAgentTransferAgentHttpportStub.Determinant req = 
		    new TransferAgentTransferAgentHttpportStub.Determinant();
		req.setMatrix(args[1]);
		
		TransferAgentTransferAgentHttpportStub.DeterminantResponse res=
                    stub.determinant(req);
		System.out.println("Determinant is:    " + res.get_return());
	    }
	    
	    else if (args[0].equalsIgnoreCase("valence")){
		TransferAgentTransferAgentHttpportStub.Valence req = 
		    new TransferAgentTransferAgentHttpportStub.Valence();
		req.setMatrix(args[1]);
		
		TransferAgentTransferAgentHttpportStub.ValenceResponse res = 
		    stub.valence(req);
		System.out.println("Valence is:    " + res.get_return());
	    }
	    
	    else if (args[0].equalsIgnoreCase("trace")){
		TransferAgentTransferAgentHttpportStub.Trace req = 
		    new TransferAgentTransferAgentHttpportStub.Trace();
		req.setMatrix(args[1]);
		
		TransferAgentTransferAgentHttpportStub.TraceResponse res = 
		    stub.trace(req);
		System.out.println("Trace is:    " + res.get_return());
	    }
	    else
		{
		    System.out.println("Invalid method");
		}
	}
	
	catch(Exception e)
	    {
		e.printStackTrace();
	    }
    }
}
