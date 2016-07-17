package gb.esac.binner;

import org.apache.log4j.Logger;


public class Bin extends AbstractBin {

    private static Logger logger  = Logger.getLogger(Bin.class);

    //  Constructors

    private Bin() {

    }

    private Bin(Bin bin) {
	setEdges(bin.getLeftEdge(), bin.getRightEdge());
    }

    public Bin(double leftEdge, double rightEdge) {
	setEdges(leftEdge, rightEdge);
    }

    public Bin[] split(double whereToSplit) throws BinningException {
	if ( !this.covers(whereToSplit) ) {
	    throw new BinningException("Bin does not cover this value: Cannot split");
	}
	Bin leftBin = new Bin(this.leftEdge, whereToSplit);
	Bin rightBin = new Bin(whereToSplit+Math.ulp(whereToSplit), this.rightEdge);
	return new Bin[] {leftBin, rightBin};
    }

    public Bin join(IBin bin) {

	Bin newBin = null;
	if ( this.covers(bin) ) {
	    newBin = new Bin(this);
	}
	if ( bin.covers(this) ) {
	    newBin = new Bin((Bin)bin);
	}
	if ( this.overlaps(bin) ) {
	    double newLeftEdge = Math.min(this.leftEdge, bin.getLeftEdge());
	    double newRightEdge = Math.max(this.rightEdge, bin.getRightEdge());
	    newBin = new Bin(newLeftEdge, newRightEdge);
	}
	return newBin;
    }


}