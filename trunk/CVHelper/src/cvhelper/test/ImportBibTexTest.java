package cvhelper.test;

import junit.framework.*;
import cvhelper.cv.*;
import java.io.*;
import java.util.*;

public class ImportBibTexTest extends TestCase {
    public ImportBibTexTest() {
        super();
    }

    public ImportBibTexTest(String name) {
        super(name);
    }

    public void testImport() {
        String filename = "testBibTeXfile.bib";
        List resL;
        CV cv = new CV();

        try {
            FileReader reader = new FileReader(filename);
            resL = ImportBibTex.parseBibTeX(reader);

            // there should be 3 elements, then two integers, the first with value 1, the second with value 0
            Iterator it = resL.iterator();

            // check that the first is in fact a book entry
            Assert.assertTrue(it.hasNext());
            Object o = it.next();
            Assert.assertTrue(o instanceof Entry);
            BookEntry b = new BookEntry((Entry)o);

            b.setCategory("x");

            Assert.assertTrue(b.isValid());
            Assert.assertEquals(b.getTitle(), "As the Agent Turns.");

            List l = b.getAuthor();
            Assert.assertTrue(l.size() == 1 && l.contains(new Name("Keith", "Decker")));

            Assert.assertEquals(b.getPublisher(), "Westing-House");
            Assert.assertEquals(b.getYear(), 2004);

            // check that the second is in fact an article
            Assert.assertTrue(it.hasNext());
            o = it.next();
            Assert.assertTrue(o instanceof Entry);
            ArticleEntry a = new ArticleEntry((Entry)o);

            a.setCategory("x");

            Assert.assertTrue(a.isValid());
            Assert.assertEquals(a.getTitle(), "More stuff about agents (the not secret kind)");

            l = b.getAuthor();
            Assert.assertTrue(l.size() == 1 && l.contains(new Name("Keith", "Decker")));

            Assert.assertEquals(a.getJournal(), "That Journal that publishes agent stuff");
            Assert.assertEquals(a.getYear(), 2001);

            Assert.assertTrue(a.hasTitle() && a.hasYear() && a.hasAuthor() && !a.hasVolume());

            Assert.assertEquals(a.getMonth, "June");

            Assert.assertTrue(it.hasNext());
            Object o = it.next();
            Assert.assertTrue(o instanceof Entry);
            PhDThesisEntry p = new PhDThesisEntry((Entry)o);

            p.setCategory("x");

            Assert.assertTrue(p.isValid());
            Assert.assertEquals(p.getTitle(), "Agents are a good idea");

            l = p.getAuthor();
            Assert.assertTrue(l.size() == 1 && l.contains(new Name("Keith", "Decker")));

            Assert.assertEquals(p.getSchool(), "University of Pittsburgh");
            Assert.assertEquals(p.getYear(), 1992);

            Assert.assertEquals(p.getType(), "A long one with lots of graphs");

            Assert.assertTrue(p.hasTitle() && p.hasYear() && p.hasAuthor() && !p.hasAddress());

            Assert.assertTrue(it.hasNext());
            o = it.next();
            Assert.assertTrue(o instanceof Integer);
            Integer inte = (Integer) o;

            Assert.assertEquals(inte.getInt(), 1);

            Assert.assertTrue(it.hasNext());
            o = it .next();
            Assert.assertTrue(o instanceof Integer);
            inte = (Integer) o;

            Assert.assertEquals(inte.getInt(), 0);

            Assert.assertTrue(!it.hasNext());

            // that's it


        }
        catch(IOException e) {
            Assert.assertTrue("Couldn't open the test file", false);
        }

    }


    public static Test suite() {

        TestSuite suite = new TestSuite();

        suite.addTest(new ImportBibTexTest("testImport"));

        return suite;
    }
}
