package samples.quickstart.clients;

import samples.quickstart.service.adb.TransferAgentTransferAgentSOAP12PortStub;

public class TransferAgentClient {

    public static void main(String[] args)
    {
	try
	    {
		TransferAgentTransferAgentSOAP12PortStub stub = new TransferAgentTransferAgentSOAP12PortStub("http://hmrg.pc.cis.udel.edu:2007/axis2/services/TransferAgent");

		request(stub, args);
	    }

	catch(Exception e)
	    {
		e.printStackTrace();
	    }
    }


    public static void request(TransferAgentTransferAgentSOAP12PortStub stub, String[] args)
    {
	try
	    {
		if (args[0].equalsIgnoreCase("rank"))
		    {
			TransferAgentTransferAgentSOAP12PortStub.Rank req = 
			    new TransferAgentTransferAgentSOAP12PortStub.Rank();
			req.setMatrix(args[1]);

			TransferAgentTransferAgentSOAP12PortStub.RankResponse res = 
			    stub.rank(req);
			System.out.println("Rank is:    " + res.get_return());
		    }

		else if (args[0].equalsIgnoreCase("determinant"))
		    {
			TransferAgentTransferAgentSOAP12PortStub.Determinant req = 
			    new TransferAgentTransferAgentSOAP12PortStub.Determinant();
			req.setMatrix(args[1]);

			TransferAgentTransferAgentSOAP12PortStub.DeterminantResponse res = 
			    stub.determinant(req);
			System.out.println("Determinant is:    " + res.get_return());
		    }

		else if (args[0].equalsIgnoreCase("valence"))
		    {
			TransferAgentTransferAgentSOAP12PortStub.Valence req = 
			    new TransferAgentTransferAgentSOAP12PortStub.Valence();
			req.setMatrix(args[1]);

			TransferAgentTransferAgentSOAP12PortStub.ValenceResponse res = 
			    stub.valence(req);
			System.out.println("Valence is:    " + res.get_return());
		    }

		else if (args[0].equalsIgnoreCase("trace"))
		    {
			TransferAgentTransferAgentSOAP12PortStub.Trace req = 
			    new TransferAgentTransferAgentSOAP12PortStub.Trace();
			req.setMatrix(args[1]);

			TransferAgentTransferAgentSOAP12PortStub.TraceResponse res = 
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
