package cvhelper.cv;

import java.io.*;
import java.util.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class TemplateGen {
    public TemplateGen() {
    }

    public static void printHeader(Writer w, String className) {
        try {
            w.write(
                    "package cvhelper.cv;\n\nimport bibtex.dom.*;\n\n/**\n * " + className + "\n * Extension of Entry class, implements special handling rules of " + className + " BibTeX entry\n *\n */\npublic class " +
                    className + " extends Entry { \n\n\t/**\n\t * " + className + "\n\t * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.\n\t *\n\t */\n\tpublic " + className +
                    "(BibtexEntry entry);\n\n\t/**\n\t * " + className + "\n\t * Basic constructor, a kind of \"copy\" constructor between BibTeX types, converts the underlying element to this type\n\t *\n\t */\n\tpublic class " + className +
                    "(Entry entry);\n\n");

        }
        catch(IOException e) {
            System.err.println("Could not writer header");
        }
    }

    public static void printReqSetter(Writer W, String className, String att, String type) {

        try {
            W.write("\t/**\n\t * set" + att + "\n\t * Sets the value of the " + att + " BibTeX field in the BibTeX Entry\n\t * Note that the field " + att + " is a required field of this BibTeX Type\n\t *\n\t */");

            if (type == "List") {

                W.write("\n\tpublic void set" + att + "(List L);");
            }
            else if(type == "String") {
                W.write("\n\tpublic void set" + att + "(String " + att + ");");
            }
            else if(type == "int") {
                W.write("\n\tpublic void set" + att + "(int " + att + ");");
            }
        }
        catch(IOException e) {
            System.err.println("Could write attribute line for " + att);
        }

    }

    public static void printSetter(Writer W, String className, String att, String type) {

        try {
            W.write("\t/**\n\t * set" + att + "\n\t * Sets the value of the " + att + " BibTeX field in the BibTeX Entry\n\t *\n\t */");

            if (type == "List") {

                W.write("\n\tpublic void set" + att + "(List L);");
            }
            else if(type == "String") {
                W.write("\n\tpublic void set" + att + "(String " + att + ");");
            }
            else if(type == "int") {
                W.write("\n\tpublic void set" + att + "(int " + att + ");");
            }
        }
        catch(IOException e) {
            System.err.println("Could write attribute line for " + att);
        }

    }


    public static void printGetter(Writer W, String className, String att, String type) {
        try {

            W.write("\t/**\n\t * get" + att + "\n\t * Gets the current value of the " + att + " BibTeX field in the BibTeX Entry\n\t *\n\t */");


            if(type == "List") {
                W.write("\n\tpublic List get" + att + "();");
            }
            else if(type == "String") {
                W.write("\n\tpublic String get" + att + "();");
            }
            else if(type == "int") {
                W.write("\n\tpublic int get" + att + "();");
            }
        }
        catch(IOException e) {
            System.err.println("Could not write getter for " + att);
        }
    }

    public static void printReqGetter(Writer W, String className, String att, String type) {
        try {

            W.write("\t/**\n\t * get" + att + "\n\t * Gets the current value of the " + att + " BibTeX field in the BibTeX Entry\n\t * Note the " + att + " field is a required field of this type of BibTeX entry\n\t *\n\t */");


            if(type == "List") {
                W.write("\n\tpublic List get" + att + "();");
            }
            else if(type == "String") {
                W.write("\n\tpublic String get" + att + "();");
            }
            else if(type == "int") {
                W.write("\n\tpublic int get" + att + "();");
            }
        }
        catch(IOException e) {
            System.err.println("Could not write getter for " + att);
        }
    }


    public static void printHas(Writer W, String className, String att) {
        try {
            W.write("\t/**\n\t * has" + att + "\n\t * Checks whether the " + att + " BibTeX field in the BibTeX Entry has been initialized\n\t *\n\t */");


            W.write("\n\tpublic boolean has" + att + "();");
        }
        catch(IOException e) {
            System.err.println("Could not generate has for " + att);
        }
    }

    public static void printValid(Writer W, String className, List L) {
        Iterator it;


        try {
            W.write("\t/**\n\t * isValid\n\t * Checks whether all the requried fields of this entry have been defined using the has{Field} methods.\n\t *\n\t */");
            W.write("\n\tpublic boolean isValid() {\n\t\treturn super.isValid()");

            it = L.iterator();
            while(it.hasNext()) {
                String att = (String) it.next();
                W.write(" && has" + att + "()");
            }
            W.write(";\n\t}\n");
        }
        catch(IOException e) {
            System.err.println("Could not write isValid section.");
        }
    }

    public static void printFooter(Writer W) {
        try { W.write("}"); } catch(IOException e) { System.err.println("Could not print footer");}
    }

    public static void staticGen(String path, String className, List ReqNames, List ReqTypes, List OptNames, List OptTypes) {

        try {
            FileWriter f = new FileWriter(path + className + ".java");

            printHeader(f, className);

            Iterator nameIt = ReqNames.iterator();
            Iterator typeIt = ReqTypes.iterator();

            while(nameIt.hasNext() && typeIt.hasNext()) {
                String name = (String) nameIt.next();
                String type = (String) typeIt.next();

                printReqSetter(f, className, name, type);
                f.write("\n\n");

                printReqGetter(f, className, name, type);
                f.write("\n\n");

                printHas(f, className, name);
                f.write("\n\n");
            }

            nameIt = OptNames.iterator();
            typeIt = OptTypes.iterator();

            while(nameIt.hasNext() && typeIt.hasNext()) {
                String name = (String) nameIt.next();
                String type = (String) typeIt.next();

                printSetter(f, className, name, type);
                f.write("\n\n");

                printGetter(f, className, name, type);
                f.write("\n\n");

                printHas(f, className, name);
                f.write("\n\n");

            }

            printValid(f, className, ReqNames);
            f.write("\n\n");

            printFooter(f);

            f.close();

        }
        catch(IOException e) {
            System.err.println("Could not open file for output");
        }
    }

    public static void main(String[] argv) {
        List ReqNames = new LinkedList();
        List ReqTypes = new LinkedList();

        List OptNames = new LinkedList();
        List OptTypes = new LinkedList();

        ReqNames.add("Title");
        ReqTypes.add("String");

        ReqNames.add("Author");
        ReqTypes.add("List");

        ReqNames.add("Note");
        ReqTypes.add("String");

        OptNames.add("Month");
        OptTypes.add("String");

        OptNames.add("Year");
        OptTypes.add("int");

        staticGen("/Users/rich/Desktop/", "UnPublishedEntry", ReqNames, ReqTypes, OptNames, OptTypes);
    }

}
