// TransferAgentMiddleman.java
// Matthew Fendt, Summer '07-'08 Research
// The intermediate web service that the users call

package samples.quickstart.middleman;

public class TransferAgentMiddleman
{
    public int enqueueAndGetID(String matrix, String operation, 
			       String requestIP)
    {
	return -1;
    }

    public String requestrank(int requestID, String requestIP)
    {
	return "Rank answer";
    }

    public double[] getTimeEstimate(int requestID)
    {
	double[] answer = {-1.0,-1.0};
	return answer;
    }

    public String requestdeterminant(int requestID, String requestIP)
    {
	return "Det answer";
    }

    public String requestvalence(int requestID, String requestIP)
    {
	return "Valence answer";
    }

    public String requesttrace(int requestID, String requestIP)
    {
	return "Trace answer";
    }

    public String requestsmithNormalForm(int requestID, String requestIP)
    {
	return "SNF answer";
    }

    public int countOperationsInQueue()
    {
	return 0;
    }

    public int checkID(int ID)
    {
	return 0;
    }
}
