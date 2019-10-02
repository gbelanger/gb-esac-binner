package gb.esac.binner;

import org.apache.log4j.Logger;

/**

   The class <code>AbsoluteQuantityBin</code> extend <code>AbstractIntensityBin</code> 
   and represents an actual intensity bin. 

   @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
   @created March 2013
   @version August 2018

**/


public class AbsoluteQuantityBin extends AbstractIntensityBin {

    private static Logger logger  = Logger.getLogger(AbsoluteQuantityBin.class);

    // Constructors

    private AbsoluteQuantityBin() {
	super();
    }

    //  with another AbsoluteQuantityBin    
    public AbsoluteQuantityBin(AbsoluteQuantityBin absoluteQuantityBin) throws BinningException {
	super(absoluteQuantityBin);
    }

    //  with IIntensity and IBin
    public AbsoluteQuantityBin(IIntensity iIntensity, IBin iBin) throws BinningException {
	super(iIntensity, iBin);
    }

    //  with value
    public AbsoluteQuantityBin(double leftEdge, double rightEdge, double absoluteQuantity) throws BinningException {
	super(leftEdge, rightEdge, absoluteQuantity);
    }
    public AbsoluteQuantityBin(double leftEdge, double rightEdge, double absoluteQuantity, String units) throws BinningException {
	super(leftEdge, rightEdge, absoluteQuantity, units);
    }
    public AbsoluteQuantityBin(double leftEdge, double rightEdge, double absoluteQuantity, String units, String description) throws BinningException {
	super(leftEdge, rightEdge, absoluteQuantity, units, description);
    }


    //  Methods for splitting absoluteQuantity bins
    public AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin previousBin, AbsoluteQuantityBin nextBin) throws BinningException {
	return IntensityBinSplitter.split(whereToSplit, this, previousBin, nextBin);
    }
    
    public AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin previousBin, AbsoluteQuantityBin nextBin, boolean addNoise) throws BinningException {
	return IntensityBinSplitter.split(whereToSplit, this, previousBin, nextBin, addNoise);
    }

    public AbsoluteQuantityBin[] splitFirstBin(double whereToSplit, AbsoluteQuantityBin nextBin) throws BinningException {
	return IntensityBinSplitter.splitFirstBin(whereToSplit, this, nextBin);
    }

    public AbsoluteQuantityBin[] splitLastBin(double whereToSplit, AbsoluteQuantityBin previousBin) throws BinningException {
	return IntensityBinSplitter.splitLastBin(whereToSplit, this, previousBin);
    }


    //  Method for joining density bins
    public AbsoluteQuantityBin joinWith(AbsoluteQuantityBin absoluteQuantityBin) throws BinningException {
	return IntensityBinCombiner.join(this, absoluteQuantityBin);
    }

    
}
