package gb.esac.binner;

import org.apache.log4j.Logger;


public class DataBin extends AbstractDataBin {

    private static Logger logger  = Logger.getLogger(DataBin.class);
    private static DataBinSplitter dataBinSplitter = new DataBinSplitter();

    // Constructors

    private DataBin() {

    }

    public DataBin(DataBin dataBin) {
	setEdges(dataBin.getEdges());
	setIntensity(dataBin.getIntensity());
	setError(dataBin.getError());
    }

    public DataBin(double leftEdge, double rightEdge, double intensity) {
	setEdges(leftEdge, rightEdge);
	setIntensity(intensity);
    }

    public DataBin(double leftEdge, double rightEdge, double intensity, double error) {
	setEdges(leftEdge, rightEdge);
	setIntensity(intensity);
	setError(error);
    }


    //  Abstract methods in AbstractBin that need to be implemented according to the child class

    public DataBin[] split(double whereToSplit, DataBin previousDataBin, DataBin nextDataBin) throws BinningException {
	return dataBinSplitter.split(whereToSplit, this, previousDataBin, nextDataBin);
    }

    public DataBin[] splitFirstBin(double whereToSplit, DataBin nextDataBin) throws BinningException {
	return dataBinSplitter.splitFirstBin(whereToSplit, this, nextDataBin);
    }

    public DataBin[] splitLastBin(double whereToSplit, DataBin previousDataBin) throws BinningException {
	return dataBinSplitter.splitLastBin(whereToSplit, this, previousDataBin);
    }


//     public DataBin join(IBin bin) {

// 	DataBin newBin = null;
// 	if ( this.covers(bin) ) {

// 	}
// 	if ( bin.covers(this) ) {

// 	}
// 	if ( this.overlaps(bin) ) {
// 	    double newLeftEdge = Math.min(this.leftEdge, bin.getLeftEdge());
// 	    double newRightEdge = Math.max(this.rightEdge, bin.getRightEdge());

// 	}
// 	return newBin;
//     }

}