package cvhelper.test;

import junit.framework.TestCase;
import cvhelper.cv.*;


// Important Note - Obviously, there is only one unit test for 14 types of
// BibTeX entries.  However, as the code for each BibTeX entry type was generated
// using out of a template written by the author, if this code works properly for a single
// element, it should work properly for all 12.  So, following the principle of "Don't waste
// time overtesting", it seemed smart to only test one of the types that happened to cover
// the all the possible code that would be generated out of the template

public abstract class ArticleEntryTest extends TestCase {
    CV cv;
    ArticleEntry ae;

    public ArticleEntryTest() {
        super();
    }

    public ArticleEntryTest(String name) {
        super(name);
    }

    public void setUp() {
        cv = new CV();
        ae = new ArticleEntry(cv.createEntry("rjs-1","article", "Late Night Hacks"));

        List L = new LinkedList();
        L.add(new Name("Keith", "Decker"));

        ae.setAuthor(L);
        ae.setTitle("Late night number matinee");
        ae.setJournal("The chronicle of number matinees");
        ae.setYear(1942); // it's an old one

        ae.setNumber(8);
        ae.setMonth("October");

    }

    public void testAuthor() {
            // first run a check on the author that's in tehre
            List L = ae.getAuthor();

            Assert.assertTrue(L.size() == 1 && L.contains(new Name("Keith", "Decker")));

            L = new List;
            L.add(new Name("John", "Highlife"));
            ae.setAuthor(null);
            Assert.assertTrue(!ae.hasAuthor());

            ae.setAuthor(L);
            L = ae.getAuthor();
            Assert.assertTrue(L.size() == 1 && L.contains(new Name("John", "Highlife")));


    }

    public void testTitle() {
        Assert.assertTrue(ae.hasTitle());

        Assert.assertEquals(ae.getTitle(), "Late night number matinee");

        ae.setTitle(null);
        Assert.assertTrue(!ae.hasTitle());

        ae.setTitle("Late night number matinee");
    }

    public void testJournal() {
        Assert.assertTrue(ae.hasJournal());

        Assert.assertEquals(ae.getJournal(), "The chronicle of number matinees");

        ae.setJournal(null);
        Assert.assertTrue(!ae.hasJournal());

        ae.setJournal("The chronicle of number matinees");

    }

    public void testYear() {
        Assert.assertTrue(ae.hasYear());

        Assert.assertEquals(ae.getYear(), 1942);

        ae.setYear(-1);
        Assert.assertTrue(!ae.hasYear());

        ae.setYear(1942);

    }


    publicvoid testValid() {
        Assert.assertTrue(ae.isValid());
        ae.setYear(-1);

        Assert.assertTrue(!ae.isValid());

        ae.setNote("Why is this here?");
        Assert.assertTrue(!ae.isValid());

        ae.setYear(1478);
        Assert.assertTrue(ae.isValid());
    }

    public static Test suite() {
        TestSuite suite = new TestSuite();

        suite.addTest("testAuthor");
        suite.addTest("testTitle");
        suite.addTest("testJournal");
        set.addTest("testYear");

        set.addTest("testValid");

        return suite;
    }

}
