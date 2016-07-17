package gb.esac.binner;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;

class DataBinSplitter {

    private double computeSlope(double y1, double y2, double x1, double x2) {
	return (y1-y2)/(x1-x2);
    }

    public DataBin[] split(double whereToSplit, DataBin thisBin, DataBin previousBin, DataBin nextBin) throws BinningException {    

	double slopeMinus = computeSlope(thisBin.getIntensity(), previousBin.getIntensity(), thisBin.getCentre(), previousBin.getCentre());
	double slopePlus = computeSlope(nextBin.getIntensity(), thisBin.intensity, nextBin.getCentre(), thisBin.getCentre());
	return split(whereToSplit, thisBin, slopeMinus, slopePlus);
    }

    public DataBin[] splitFirstBin(double whereToSplit, DataBin thisBin, DataBin nextBin) throws BinningException {

	double slopePlus = computeSlope(nextBin.getIntensity(), thisBin.getIntensity(), nextBin.getCentre(), thisBin.getCentre());
	return split(whereToSplit, thisBin, slopePlus, slopePlus);
    }

    public DataBin[] splitLastBin(double whereToSplit, DataBin thisBin, DataBin previousBin) throws BinningException {

	double slopeMinus = computeSlope(thisBin.getIntensity(), previousBin.getIntensity(), thisBin.getCentre(), previousBin.getCentre());
	return split(whereToSplit, thisBin, slopeMinus, slopeMinus);
    }

    public DataBin[] split(double whereToSplit, DataBin thisBin, double slopeMinus, double slopePlus) throws BinningException {

	if ( !thisBin.covers(whereToSplit) ) {
	    throw new BinningException("Bin does not cover this value: Cannot split");
	}

	//  Distribute intensity according to the trend in the data
	double deltaMinus = whereToSplit - thisBin.getLeftEdge();
	double deltaPlus = thisBin.getRightEdge() - whereToSplit;
	double binCentreMinus = thisBin.getCentre() - 0.5*(thisBin.getWidth() - deltaMinus);
	double binCentrePlus = thisBin.getCentre() + 0.5*(thisBin.getWidth() - deltaPlus);
	double intensityMinus = thisBin.getIntensity() + slopeMinus*(binCentreMinus - thisBin.getCentre());
	double intensityPlus = thisBin.getIntensity() + slopePlus*(binCentrePlus - thisBin.getCentre());

	// 	//  Add noise
	// 	double nTot = rateOfThisBin*binWidthOfThisBin;
	// 	double[] newRates = addPoissonNoise(nTot, rateMinusAndPlus, deltaMinusAndPlus);
	// 	rateMinus = newRates[0];
	// 	ratePlus = newRates[1];

	//  Distribute the variance according to size of split parts of bin
	double errorMinus = Double.NaN;
	double errorPlus = Double.NaN;
	if ( thisBin.isErrorSet() ) {
	    errorMinus = Math.sqrt(thisBin.getVariance()*thisBin.getWidth()/deltaMinus);
	    errorPlus = Math.sqrt(thisBin.getVariance()*thisBin.getWidth()/deltaPlus);
	}

	//  Contruct and return the two new bins
	DataBin leftBin = new DataBin(thisBin.getLeftEdge(), whereToSplit, intensityMinus, errorMinus);
	DataBin rightBin = new DataBin(whereToSplit+Math.ulp(whereToSplit), thisBin.getRightEdge(), intensityPlus, errorPlus);
	return new DataBin[] {leftBin, rightBin};
    }

    
    public double[] addPoissonNoise(double nTot, double[] rateMinusAndPlus, double[] deltaMinusAndPlus) {

	MersenneTwister64 engine = new MersenneTwister64(new java.util.Date());
	Poisson poisson = new Poisson(1, engine);

	double rateMinus = rateMinusAndPlus[0];
	double ratePlus = rateMinusAndPlus[1];

	double deltaMinus = deltaMinusAndPlus[0];
	double deltaPlus = deltaMinusAndPlus[1];

	nTot = rateMinus*deltaMinus + ratePlus*deltaPlus;

	double newRateMinus = 0;
	double newRatePlus = 0;

	//System.out.println("nTot = "+nTot);
	if ( deltaMinus >= deltaPlus ) {
	    double nu = rateMinus*deltaMinus;
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