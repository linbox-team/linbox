package cvhelper.cv;

import bibtex.dom.*;

/**
 * InCollection
 * Extension of Entry class, implements special handling rules of InCollection BibTeX entry
 *
 */
public class InCollectionEntry extends Entry {

	/**
	 * InCollection
	 * Basic constructor, creates a new Entry using the underlying BibTeXEntry element.
	 *
	 */
	InCollection(BibtexEntry entry);

	/**
	 * InCollection
	 * Basic constructor, a kind of "copy" constructor between BibTeX types, converts the underlying element to this type
	 *
	 */
	public class InCollection(Entry entry);

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
	 * setBookTitle
	 * Sets the value of the BookTitle BibTeX field in the BibTeX Entry
	 * Note that the field BookTitle is a required field of this BibTeX Type
	 *
	 */
	public void setBookTitle(String BookTitle);

	/**
	 * getBookTitle
	 * Gets the current value of the BookTitle BibTeX field in the BibTeX Entry
	 * Note the BookTitle field is a required field of this type of BibTeX entry
	 *
	 */
	public String getBookTitle();

	/**
	 * hasBookTitle
	 * Checks whether the BookTitle BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasBookTitle();

	/**
	 * setPublisher
	 * Sets the value of the Publisher BibTeX field in the BibTeX Entry
	 * Note that the field Publisher is a required field of this BibTeX Type
	 *
	 */
	public void setPublisher(String Publisher);

	/**
	 * getPublisher
	 * Gets the current value of the Publisher BibTeX field in the BibTeX Entry
	 * Note the Publisher field is a required field of this type of BibTeX entry
	 *
	 */
	public String getPublisher();

	/**
	 * hasPublisher
	 * Checks whether the Publisher BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasPublisher();

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
	 * setEditor
	 * Sets the value of the Editor BibTeX field in the BibTeX Entry
	 *
	 */
	public void setEditor(List L);

	/**
	 * getEditor
	 * Gets the current value of the Editor BibTeX field in the BibTeX Entry
	 *
	 */
	public List getEditor();

	/**
	 * hasEditor
	 * Checks whether the Editor BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasEditor();

	/**
	 * setVolumeOrNumber
	 * Sets the value of the VolumeOrNumber BibTeX field in the BibTeX Entry
	 *
	 */
	public void setVolumeOrNumber(int VolumeOrNumber);

	/**
	 * getVolumeOrNumber
	 * Gets the current value of the VolumeOrNumber BibTeX field in the BibTeX Entry
	 *
	 */
	public int getVolumeOrNumber();

	/**
	 * hasVolumeOrNumber
	 * Checks whether the VolumeOrNumber BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasVolumeOrNumber();

	/**
	 * setSeries
	 * Sets the value of the Series BibTeX field in the BibTeX Entry
	 *
	 */
	public void setSeries(String Series);

	/**
	 * getSeries
	 * Gets the current value of the Series BibTeX field in the BibTeX Entry
	 *
	 */
	public String getSeries();

	/**
	 * hasSeries
	 * Checks whether the Series BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasSeries();

	/**
	 * setType
	 * Sets the value of the Type BibTeX field in the BibTeX Entry
	 *
	 */
	public void setType(String Type);

	/**
	 * getType
	 * Gets the current value of the Type BibTeX field in the BibTeX Entry
	 *
	 */
	public String getType();

	/**
	 * hasType
	 * Checks whether the Type BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasType();

	/**
	 * setChapter
	 * Sets the value of the Chapter BibTeX field in the BibTeX Entry
	 *
	 */
	public void setChapter(int Chapter);

	/**
	 * getChapter
	 * Gets the current value of the Chapter BibTeX field in the BibTeX Entry
	 *
	 */
	public int getChapter();

	/**
	 * hasChapter
	 * Checks whether the Chapter BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasChapter();

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
	 * setAddress
	 * Sets the value of the Address BibTeX field in the BibTeX Entry
	 *
	 */
	public void setAddress(String Address);

	/**
	 * getAddress
	 * Gets the current value of the Address BibTeX field in the BibTeX Entry
	 *
	 */
	public String getAddress();

	/**
	 * hasAddress
	 * Checks whether the Address BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasAddress();

	/**
	 * setEdition
	 * Sets the value of the Edition BibTeX field in the BibTeX Entry
	 *
	 */
	public void setEdition(int Edition);

	/**
	 * getEdition
	 * Gets the current value of the Edition BibTeX field in the BibTeX Entry
	 *
	 */
	public int getEdition();

	/**
	 * hasEdition
	 * Checks whether the Edition BibTeX field in the BibTeX Entry has been initialized
	 *
	 */
	public boolean hasEdition();

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
		return super.isValid() && hasTitle() && hasAuthor() && hasBookTitle() && hasPublisher() && hasYear();
	}

        // This method is a simple extension to Entry's clone method
        public Object clone() { return null;}

}
