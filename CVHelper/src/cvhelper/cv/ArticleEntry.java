package cvhelper.cv;

import bibtex.dom.*;

/**
 * ArticleEntry
 * Extension of Entry class, implements special handling rules of ArticleEntry BibTeX entry
 *
 */
public class ArticleEntry extends Entry {

	/**
	 * ArticleEntry
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	ArticleEntry(BibtexEntry entry);

	/**
	 * ArticleEntry
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class ArticleEntry(Entry entry);

	/**
	 * setAuthor
	 * Sets the value of the Author BibTeX field in the BibTeX Entry
	 * Note that the field Author is a required field of this BibTeX Type
	 *
	 */
	public void setAuthor(List L);

	/**
	 * getAuthor
	 * Gets the current value of the Author BibTeX field in the BibTeX Entry
	 * Note the Author field is a required field of this type of BibTeX entry
	 *
	 */
	public List getAuthor();

	/**
	 * hasAuthor
	 * Checks whether the Author BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasAuthor();

	/**
	 * setTitle
	 * Sets the value of the Title BibTeX field in the BibTeX Entry
	 * Note that the field Title is a required field of this BibTeX Type
	 *
	 */
	public void setTitle(String Title);

	/**
	 * getTitle
	 * Gets the current value of the Title BibTeX field in the BibTeX Entry
	 * Note the Title field is a required field of this type of BibTeX entry
	 *
	 */
	public String getTitle();

	/**
	 * hasTitle
	 * Checks whether the Title BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasTitle();

	/**
	 * setPublication
	 * Sets the value of the Publication BibTeX field in the BibTeX Entry
	 * Note that the field Publication is a required field of this BibTeX Type
	 *
	 */
	public void setPublication(String Publication);

	/**
	 * getPublication
	 * Gets the current value of the Publication BibTeX field in the BibTeX Entry
	 * Note the Publication field is a required field of this type of BibTeX entry
	 *
	 */
	public String getPublication();

	/**
	 * hasPublication
	 * Checks whether the Publication BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasPublication();

	/**
	 * setYear
	 * Sets the value of the Year BibTeX field in the BibTeX Entry
	 * Note that the field Year is a required field of this BibTeX Type
	 *
	 */
	public void setYear(int Year);

	/**
	 * getYear
	 * Gets the current value of the Year BibTeX field in the BibTeX Entry
	 * Note the Year field is a required field of this type of BibTeX entry
	 *
	 */
	public int getYear();

	/**
	 * hasYear
	 * Checks whether the Year BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasYear();

	/**
	 * setVolume
	 * Sets the value of the Volume BibTeX field in the BibTeX Entry
	 *
	 */
	public void setVolume(int Volume);

	/**
	 * getVolume
	 * Gets the current value of the Volume BibTeX field in the BibTeX Entry
	 *
	 */
	public int getVolume();

	/**
	 * hasVolume
	 * Checks whether the Volume BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasVolume();

	/**
	 * setNumber
	 * Sets the value of the Number BibTeX field in the BibTeX Entry
	 *
	 */
	public void setNumber(int Number);

	/**
	 * getNumber
	 * Gets the current value of the Number BibTeX field in the BibTeX Entry
	 *
	 */
	public int getNumber();

	/**
	 * hasNumber
	 * Checks whether the Number BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasNumber();

	/**
	 * setPages
	 * Sets the value of the Pages BibTeX field in the BibTeX Entry
	 *
	 */
	public void setPages(int Pages);

	/**
	 * getPages
	 * Gets the current value of the Pages BibTeX field in the BibTeX Entry
	 *
	 */
	public int getPages();

	/**
	 * hasPages
	 * Checks whether the Pages BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasPages();

	/**
	 * setMonth
	 * Sets the value of the Month BibTeX field in the BibTeX Entry
	 *
	 */
	public void setMonth(String Month);

	/**
	 * getMonth
	 * Gets the current value of the Month BibTeX field in the BibTeX Entry
	 *
	 */
	public String getMonth();

	/**
	 * hasMonth
	 * Checks whether the Month BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasMonth();

	/**
	 * setNote
	 * Sets the value of the Note BibTeX field in the BibTeX Entry
	 *
	 */
	public void setNote(String Note);

	/**
	 * getNote
	 * Gets the current value of the Note BibTeX field in the BibTeX Entry
	 *
	 */
	public String getNote();

	/**
	 * hasNote
	 * Checks whether the Note BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasNote();

	/**
	 * isValid
	 * Checks whether all the requried fields of this entry have been defined using the has{Field} methods.
	 *
	 */
	public boolean isValid() {
		return super.isValid() && hasAuthor() && hasTitle() && hasPublication() && hasYear();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}


}
