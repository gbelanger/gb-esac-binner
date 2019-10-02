package gb.esac.binner;

import org.apache.log4j.Logger;

/**

   The class <code>DensityBin</code> extend <code>AbstractIntensityBin</code> 
   and represents an actual density bin. This means that the quantity is 
   something that needs to be multiplied with the bin width in order to 
   represent the absolute quanity. This implies that combining density bins
   must be done in particular way to take this into account. This is why
   there is a <code>IntensityBinSplitter</code> and a <code>IntensityBinCombiner</code>.

   @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
   @created March 2013
   @version August 2018

**/


public class DensityBin extends AbstractIntensityBin {

    private static Logger logger  = Logger.getLogger(DensityBin.class);
    private static IntensityBinSplitter intensityBinSplitter = new IntensityBinSplitter();

    // Constructors

    private DensityBin() {
	super();
    }

    //  with another DensityBin
    public DensityBin(DensityBin densityBin) throws BinningException {
	super(densityBin);
    }

    //  with IIntensity and IBin
    public DensityBin(IIntensity iDensity, IBin iBin) throws BinningException {
	super(iDensity,iBin);
    }

    //  with value
    public DensityBin(double leftEdge, double rightEdge, double density) throws BinningException {
	super(leftEdge, rightEdge, density);
    }
    public DensityBin(double leftEdge, double rightEdge, double density, String units) throws BinningException {
	super(leftEdge, rightEdge, density, units);
    }
    public DensityBin(double leftEdge, double rightEdge, double density, String units, String description) throws BinningException {	
	super(leftEdge, rightEdge, density, units, description);	
    }

    // with value and error
    public DensityBin(double leftEdge, double rightEdge, double density, double error) throws BinningException {
	super(leftEdge, rightEdge, density, error);
    }
    public DensityBin(double leftEdge, double rightEdge, double density, double error, String units) throws BinningException {
	super(leftEdge, rightEdge, density, error, units);
    }
    public DensityBin(double leftEdge, double rightEdge, double density, double error, String units, String description) throws BinningException {
	super(leftEdge, rightEdge, density, error, units, description);
    }

    public Density getIntensity() {
	return new Density(this.getValue(),this.getError(),this.getUnits(),this.getDescription());
    }
    
    //  Methods for splitting
    public DensityBin[] split(double whereToSplit, DensityBin previousBin, DensityBin nextBin) throws BinningException {
	return IntensityBinSplitter.split(whereToSplit, this, previousBin, nextBin);
    }

    public DensityBin[] split(double whereToSplit, DensityBin previousBin, DensityBin nextBin, boolean addNoise) throws BinningException {
	return IntensityBinSplitter.split(whereToSplit, this, previousBin, nextBin, addNoise);
    }

    public DensityBin[] splitFirstBin(double whereToSplit, DensityBin nextBin) throws BinningException {
	return IntensityBinSplitter.splitFirstBin(whereToSplit, this, nextBin);
    }

    public DensityBin[] splitLastBin(double whereToSplit, DensityBin previousBin) throws BinningException {
	return IntensityBinSplitter.splitLastBin(whereToSplit, this, previousBin);
    }

    //  Method for joining
    public DensityBin joinWith(DensityBin densityBin) throws BinningException {
	return IntensityBinCombiner.join(this, densityBin);
    }


}
