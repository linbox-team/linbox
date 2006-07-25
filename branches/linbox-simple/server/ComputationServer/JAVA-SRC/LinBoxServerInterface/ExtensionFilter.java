package LinBoxServerInterface;

import java.io.FilenameFilter;

/** A simple FilenameFilter which accepts only those files with the given
  * extenxion (e.g. xml, sms)
  */
public class ExtensionFilter implements FilenameFilter {
	
	private String extension;
	
	/** Specified by FilenameFilter interface */
	public boolean accept( java.io.File dir, String name ) {
		return name.endsWith( extension );
	}
	
	/** Constructor.
	  * @param extension The extension of files that should be accepted.
	  */
	public ExtensionFilter( String extension ) {
		this.extension = extension;
	}
}
