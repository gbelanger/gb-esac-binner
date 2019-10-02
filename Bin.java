package gb.esac.binner;

import org.apache.log4j.Logger;

/**

The class <code>Bin</code> extends <code>AbstractBin</code> which implements
<code>IBin</code>. 

Two additional methods are provided here for operations that allow to 
<code>split(double where)</code> somewhere along the length of the bin between its 
two bin edges, and to <code>join(Bin bin)</code> two bins together.

IMPORTANT: the ability to split and join bins must be universal to all bins. 
Therefore, all types of bins must extend the <code>Bin</code> class.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version August 2018

 **/

public class Bin extends AbstractBin {

    private static Logger logger  = Logger.getLogger(Bin.class);

    //  Constructors are private
    private Bin() {
	super();
    }

    public Bin(Bin bin) throws BinningException {
	super(bin);
    }

    public Bin(double leftEdge, double rightEdge) throws BinningException {
	super(leftEdge, rightEdge);
    }

    public Bin[] split(double whereToSplit) throws BinningException {
	if ( !this.contains(whereToSplit) ) {
	    throw new BinningException("Cannot split: Bin ["+this.getLeftEdge()+", "+this.getRightEdge()+"] does not contain value ("+whereToSplit+")");
	}
	Bin leftBin = new Bin(this.getLeftEdge(), whereToSplit);
	Bin rightBin = new Bin(whereToSplit, this.getRightEdge());
	return new Bin[] {leftBin, rightBin};
    }

    public Bin joinWith(IBin bin) throws BinningException {
	Bin newBin = null;
	// Adjacent bins with this one to the left
	if (this.getRightEdge() == bin.getLeftEdge()) {
	    newBin = new Bin(this.getLeftEdge(), bin.getRightEdge());
	}
	// Adjacent bins with this one to the right	
	else if (this.getLeftEdge() == bin.getRightEdge()) {
	    newBin = new Bin(bin.getLeftEdge(), this.getRightEdge());
	}
	// This contains the other
    	else if (this.contains(bin)) {
	    newBin = new Bin(this.getLeftEdge(), this.getRightEdge());
	}
	// Other contains this
	else if (bin.contains(this)) {
	    newBin = new Bin(bin.getLeftEdge(), bin.getRightEdge());
    	}
	// They overlap (this also covers perfect overlap)
    	else if (this.overlaps(bin)) {
    	    double newLeftEdge = Math.min(this.getLeftEdge(), bin.getLeftEdge());
    	    double newRightEdge = Math.max(this.getRightEdge(), bin.getRightEdge());
	    newBin = new Bin(newLeftEdge, newRightEdge);
    	}
	// Bins are disjointed
	else {
	    throw new BinningException("Cannot join: Bins are disjointed. "+
				       "Their bin edges are: ["+this.getLeftEdge()+", "+this.getRightEdge()+"] and ["+
				       bin.getLeftEdge()+", "+bin.getRightEdge()+"]");
	}
    	return newBin;
    }


}
