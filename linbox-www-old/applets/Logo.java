

import java.applet.*;	// for the class Applet that is extended by Logo
import java.awt.*;	// for classes Font, FontMetrics, Graphics 

public class Logo extends Applet{

    private final String lin = "LIN",
                         box = "BOX",
                         linbox = "LINBOX";

    public void init(){
        //font = new Font( "SansSerif", Font.BOLD, 48 );
    }

    public void paint( Graphics g ){

        Dimension appletDim = getSize(); // Get size of applet
        
        // Fill blue border
        g.setColor( Color.blue );
        g.fillRect( 0,0,appletDim.width,appletDim.height );
        int squareDim = appletDim.width < appletDim.height ? appletDim.width
                                                           : appletDim.height;
        // Create white square
        squareDim = 9*squareDim/10;	// Use of blue border of 5%
        int sqXOrig = (appletDim.width - squareDim)/2;
        int sqYOrig = (appletDim.height - squareDim)/2;
  
        g.setColor( Color.white );
        g.fill3DRect( sqXOrig,sqYOrig,squareDim,squareDim,true );

        // Create black box
        int bbDim = squareDim/2;
        int bbXOrig = sqXOrig + squareDim/4;
        int bbYOrig = sqYOrig + squareDim/4;
        g.setColor( Color.black );
        g.fill3DRect( bbXOrig,bbYOrig,bbDim,bbDim,false );
        // Create dot for multiply sumbol
        int dotXOrig = bbXOrig + bbDim + 4;
        int dotYOrig = bbYOrig + bbDim/2 - 2;
        g.fillOval( dotXOrig,dotYOrig,5,5 );
        // Create vector
        int vXTop = dotXOrig + 12;
        int vYTop = bbYOrig;
        int vXBot = vXTop;
        int vYBot = bbYOrig + bbDim;
        g.drawLine( vXTop,vYTop,vXBot,vYBot );
        g.drawLine( vXBot - 3*bbDim/32,vYBot - bbDim/8,vXBot,vYBot );
        g.drawLine( vXBot + 3*bbDim/32,vYBot - bbDim/8,vXBot,vYBot );

        // Put in LINBOX
        int fs = squareDim/8;
        Font font = new Font( "SansSerif", Font.BOLD, fs );
        g.setColor( Color.red );
        FontMetrics fm = getFontMetrics( font );
        int strXOrig = bbXOrig + (bbDim - fm.stringWidth( linbox ))/2;
        g.setFont( font );
        g.drawString( linbox,strXOrig,bbYOrig-2 );
        g.drawString( linbox,strXOrig,bbYOrig+bbDim+fm.getAscent()-fm.getDescent()+1 );
    }
}
