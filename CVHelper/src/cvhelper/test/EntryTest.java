package cvhelper.test;

import junit.framework.*;
import cvhelper.cv.*;

public class EntryTest extends TestCase {

    CV cv;
    Entry e1, e2, e3;

    public void setUp() {
        cv = new CV();

        e1 = cv.createEntry("rjs-1","article", "Stuff I wrote in Crayon");

        e2 = cv.createEntry("rjs-2", "inproceedings", "Stuff I wrote with power tools");

        e3 = cv.createEntry("rjs-1", "article", "Stuff I wrote with Crayon");

    }

    public EntryTest() {
        super();
    }

    public EntryTest(String name) {
        super(name);
    }

    public void testEquals() {
        Assert.assertEquals(e1, e3);
        Assert.assertTrue(! e1.equals(e2));
    }

    public void testType() {
        Assert.assertEquals(e1.getType(), "article");

        e2.setType("conference");
        Assert.assertEquals(e2.getType(), "conference");

        Assert.assertTrue(e3.hasType());
        e3.setType(null);
        Assert.assertTrue(!e3.hasType());

    }

    public void testCategory() {
        Assert.assertEquals(e1.getCategory(), "Stuff I wrote in Crayon");

        e2.setCategory("Stuff I wrote in Emacs");
        Assert.assertEquals(e2.getCategory(), "Stuff I wrote in Emacs");

        Assert.assertTrue(e3.hasCategory());
        e3.setCategory(null);
        Assert.assertTrue(!e3.hasCategory());

    }

    public void tesID() {
        Assert.assertEquals(e1.getID(), "rjs-1");

        e2.setCategory("2-sjr");
        Assert.assertEquals(e2.getID(), "2-sjr");

        Assert.assertTrue(e3.hasID());
        e3.setID(null);
        Assert.assertTrue(!e3.hasID());
    }

    public void testValid() {
        Assert.assertTrue(e1.isValid());
        Assert.assertTrue(e2.isValid());
        Assert.assertTrue(!e3.isValid());
    }


    public static Test suite() {
        TestSuite suite = new TestSuite();



        suite.addTest(new EntryTest("testEquals"));
        suite.addTest(new EntryTest("testType"));
        suite.addTest(new EntryTest("testCategory"));
        suite.addTest(new EntryTest("testID"));
        suite.addTest(new EntryTest("testValid"));
    }
}
