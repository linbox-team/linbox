package cvhelper.cv;

public class Name {
    private String _first, _preLast, _last, _jr;

    // every name must have a first & last name
    public Name(String First, String Last);

    // get every field in the name
    // return null if these fields are not set
    public String getFirst();
    public String getPreLast();
    public String getLast();
    public String getJr();

    // set every field in the name
    public void setFirst(String first);
    public void setPreLast(String preLast);
    public void setLast(String last);
    public void setJr(String jr);


}
