package gb.esac.binner;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;
import org.apache.log4j.Logger;

/**

The class <code>IntensityBinSplitter</code> is the class that does the work of splitting 
<code>IntensityBin</code> objects. It has the methods needed to do this, which all call 
the main one <code>split(doule, IntensityBin, slopeMinue, slopePlus)</code> that passes 
the necessary arguments for the job. The other splitting methods all call this one 
after having computed the arguments. 

All methods are public except for <code>computeSlope</code> which is used to determine 
the trend in the data. In addition, there is the <code>addPoissonNoise</code> method
that can be used or not based on the boolean switch in the arguments of the main
<code>split</code> method.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.

 @created March 2013
 @modified August 2018
 @modified January 2020

 @version January 2020

 **/

final class IntensityBinSplitter {

    private static Logger logger  = Logger.getLogger(IntensityBinSplitter.class);

    // PUBLIC

    // For DensityBin
    public static DensityBin[] splitFirstBin(double whereToSplit, DensityBin thisBin, DensityBin nextBin) throws BinningException {
		double slopePlus = computeSlope(nextBin.getValue(), thisBin.getValue(), nextBin.getCentre(), thisBin.getCentre());
		return split(whereToSplit, thisBin, slopePlus, slopePlus);
    }

    public static DensityBin[] splitLastBin(double whereToSplit, DensityBin thisBin, DensityBin previousBin) throws BinningException {
		double slopeMinus = computeSlope(thisBin.getValue(), previousBin.getValue(), thisBin.getCentre(), previousBin.getCentre());
		return split(whereToSplit, thisBin, slopeMinus, slopeMinus);
    }

    public static DensityBin[] split(double whereToSplit, DensityBin thisBin, DensityBin previousBin, DensityBin nextBin) throws BinningException {
		boolean addNoise = false; //  Default is false
		return split(whereToSplit, thisBin, previousBin, nextBin, addNoise);
    }

    public static DensityBin[] split(double whereToSplit, DensityBin thisBin, DensityBin previousBin, DensityBin nextBin, boolean addNoise) throws BinningException {
		double slopeMinus = computeSlope(thisBin.getValue(), previousBin.getValue(), thisBin.getCentre(), previousBin.getCentre());
		double slopePlus = computeSlope(nextBin.getValue(), thisBin.getValue(), nextBin.getCentre(), thisBin.getCentre());
		return split(whereToSplit, thisBin, slopeMinus, slopePlus, addNoise);
    }


    //  For AbsoluteQuantityBin
    public static AbsoluteQuantityBin[] splitFirstBin(double whereToSplit, AbsoluteQuantityBin thisBin, AbsoluteQuantityBin nextBin) throws BinningException {
		double slopePlus = computeSlope(nextBin.getValue(), thisBin.getValue(), nextBin.getCentre(), thisBin.getCentre());
		return split(whereToSplit, thisBin, slopePlus, slopePlus);
    }

    public static AbsoluteQuantityBin[] splitLastBin(double whereToSplit, AbsoluteQuantityBin thisBin, AbsoluteQuantityBin previousBin) throws BinningException {
		double slopeMinus = computeSlope(thisBin.getValue(), previousBin.getValue(), thisBin.getCentre(), previousBin.getCentre());
		return split(whereToSplit, thisBin, slopeMinus, slopeMinus);
    }

    public static AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin thisBin, AbsoluteQuantityBin previousBin, AbsoluteQuantityBin nextBin) throws BinningException {
		boolean addNoise = false; //  Default is false
		return split(whereToSplit, thisBin, previousBin, nextBin, addNoise);
    }

    public static AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin thisBin, AbsoluteQuantityBin previousBin, AbsoluteQuantityBin nextBin, boolean addNoise) throws BinningException {
		double slopeMinus = computeSlope(thisBin.getValue(), previousBin.getValue(), thisBin.getCentre(), previousBin.getCentre());
		double slopePlus = computeSlope(nextBin.getValue(), thisBin.getValue(), nextBin.getCentre(), thisBin.getCentre());
		return split(whereToSplit, thisBin, slopeMinus, slopePlus, addNoise);
    }

    
    //  PRIVATE
    
    //  For DensityBin
    
    private static DensityBin[] split(double whereToSplit, DensityBin thisBin, double slopeMinus, double slopePlus) throws BinningException {
		boolean addNoise = false; // Default is false
		return split(whereToSplit, thisBin, slopeMinus, slopePlus, addNoise);
    }
    
    private static DensityBin[] split(double whereToSplit, DensityBin thisBin, double slopeMinus, double slopePlus, boolean addNoise) throws BinningException {
		if (!thisBin.contains(whereToSplit)) {
		    throw new BinningException("Cannot split: Bin ["+thisBin.getLeftEdge()+", "+thisBin.getRightEdge()+"] does not contain "+whereToSplit+".");
		}
		else {
		    logger.info("Splitting bin ["+thisBin.getLeftEdge()+", "+thisBin.getRightEdge()+"] at "+whereToSplit);
		}
		//  Distribute intensity according to the trend in the data
		double deltaMinus = whereToSplit - thisBin.getLeftEdge();
		double deltaPlus = thisBin.getRightEdge() - whereToSplit;
		double binCentreMinus = thisBin.getCentre() - 0.5*(thisBin.getWidth() - deltaMinus);
		double binCentrePlus = thisBin.getCentre() + 0.5*(thisBin.getWidth() - deltaPlus);
		double intensityMinus = thisBin.getValue() + slopeMinus*(binCentreMinus - thisBin.getCentre());
		double intensityPlus = thisBin.getValue() + slopePlus*(binCentrePlus - thisBin.getCentre());

		//  Add Poisson noise
		if (addNoise) {
		    logger.info("  Adding Poisson noise...");
		    double nTot = thisBin.getValue()*thisBin.getWidth();
		    double[] intensityMinusAndPlus = new double[] {intensityMinus, intensityPlus};
		    double[] deltaMinusAndPlus = new double[] {deltaMinus, deltaPlus};
		    double[] noisyIntensities = addPoissonNoise(nTot, intensityMinusAndPlus, deltaMinusAndPlus);
		    intensityMinus = noisyIntensities[0];
		    intensityPlus = noisyIntensities[1];
		    logger.info("  New noisy intensities are: "+intensityMinus+" and "+intensityPlus+".");	    
		}

		//  Distribute the variance according to size of parts
		double errorMinus = Double.NaN;
		double errorPlus = Double.NaN;
		if (thisBin.errorIsSet()) {
		    errorMinus = Math.sqrt(thisBin.getVariance()*thisBin.getWidth()/deltaMinus);
		    errorPlus = Math.sqrt(thisBin.getVariance()*thisBin.getWidth()/deltaPlus);

		}
		
		//  Contruct and return the two new bins
		DensityBin leftBin = new DensityBin(thisBin.getLeftEdge(), whereToSplit, intensityMinus, errorMinus);
		DensityBin rightBin = new DensityBin(whereToSplit, thisBin.getRightEdge(), intensityPlus, errorPlus);
		//logger.info("Splitting complete");
		return new DensityBin[] {leftBin, rightBin};
    }


    
    //  For AbstoluteQuantityBin
    
    private static AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin thisBin, double slopeMinus, double slopePlus) throws BinningException {
		boolean addNoise = false; // Default is false
		return split(whereToSplit, thisBin, slopeMinus, slopePlus, addNoise);
    }
    
    private static AbsoluteQuantityBin[] split(double whereToSplit, AbsoluteQuantityBin thisBin, double slopeMinus, double slopePlus, boolean addNoise) throws BinningException {
		if (!thisBin.contains(whereToSplit)) {
		    throw new BinningException("Cannot split: Bin ["+thisBin.getLeftEdge()+", "+thisBin.getRightEdge()+"] does not contain "+whereToSplit+".");
		}
		else {
		    logger.info("Splitting bin ["+thisBin.getLeftEdge()+", "+thisBin.getRightEdge()+"] at "+whereToSplit);
		}
		//  Distribute intensity according to the trend in the data
		double deltaMinus = whereToSplit - thisBin.getLeftEdge();
		double deltaPlus = thisBin.getRightEdge() - whereToSplit;
		double binCentreMinus = thisBin.getCentre() - 0.5*(thisBin.getWidth() - deltaMinus);
		double binCentrePlus = thisBin.getCentre() + 0.5*(thisBin.getWidth() - deltaPlus);
		double intensityMinus = thisBin.getValue() + slopeMinus*(binCentreMinus - thisBin.getCentre());
		double intensityPlus = thisBin.getValue() + slopePlus*(binCentrePlus - thisBin.getCentre());
		
		//  Add Poisson noise
		if (addNoise) {
		    logger.info("  Adding Poisson noise...");
		    double nTot = thisBin.getValue()*thisBin.getWidth();
		    double[] intensityMinusAndPlus = new double[] {intensityMinus, intensityPlus};
		    double[] deltaMinusAndPlus = new double[] {deltaMinus, deltaPlus};
		    double[] noisyIntensities = addPoissonNoise(nTot, intensityMinusAndPlus, deltaMinusAndPlus);
		    intensityMinus = noisyIntensities[0];
		    intensityPlus = noisyIntensities[1];
		    logger.info("  New noisy intensities are: "+intensityMinus+" and "+intensityPlus+".");	    
		}

		//  Contruct and return the two new bins
		AbsoluteQuantityBin leftBin = new AbsoluteQuantityBin(thisBin.getLeftEdge(), whereToSplit, intensityMinus);
		AbsoluteQuantityBin rightBin = new AbsoluteQuantityBin(whereToSplit, thisBin.getRightEdge(), intensityPlus);
		//logger.info("Splitting complete");
		return new AbsoluteQuantityBin[] {leftBin, rightBin};
    }

	    
    //  General 

    private static double computeSlope(double y1, double y2, double x1, double x2) {
		return (y1-y2)/(x1-x2);
    }
	    
    private static double[] addPoissonNoise(double nTot, double[] rateMinusAndPlus, double[] deltaMinusAndPlus) {
		MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
		Poisson poisson = new Poisson(1, engine);
		// Rates
		double rateMinus = rateMinusAndPlus[0];
		double ratePlus = rateMinusAndPlus[1];
		// Widths
		double deltaMinus = deltaMinusAndPlus[0];
		double deltaPlus = deltaMinusAndPlus[1];
		// Total events
		nTot = rateMinus*deltaMinus + ratePlus*deltaPlus;
		//System.out.println("nTot = "+nTot);
		double newRateMinus = 0;
		double newRatePlus = 0;

		if (deltaMinus >= deltaPlus) {
		    double nu = rateMinus*deltaMinus;
		    //System.out.println("nu = "+nu);
		    int n = poisson.nextInt(nu);
		    //System.out.println("n = "+n);
		    double nTotMinusN = nTot - n;
		    while (nTotMinusN < 0) {
			n = poisson.nextInt(nu);
			//System.out.println("n = "+n);
			nTotMinusN = nTot - n;
			//System.out.println("nTotMinusN = "+nTotMinusN);
		    }
		    newRateMinus = n/deltaMinus;
		    newRatePlus = nTotMinusN/deltaPlus;
		    //System.out.println();
		}
		else {
		    double nu = ratePlus*deltaPlus;
		    //System.out.println("nu = "+nu);
		    int n = poisson.nextInt(nu);
		    //System.out.println("n = "+n);
		    double nTotMinusN = nTot - n;
		    while ( nTotMinusN < 0 ) {
			n = poisson.nextInt(nu);
			//System.out.println("n = "+n);
			nTotMinusN = nTot - n;
			//System.out.println("nTotMinusN = "+nTotMinusN);
		    }
		    newRatePlus = n/deltaPlus;
		    newRateMinus = nTotMinusN/deltaMinus;
		    //System.out.println();
		}
		return new double[] {newRateMinus, newRatePlus};
    }
   

}
