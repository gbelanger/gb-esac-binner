package gb.esac.binner;

import java.util.Date;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;


public final class BinResampler {

    public static double[][] resample(DataBin[] oldBins, IBin[] newBins) {

	//  The rate in each new bin is estimated using the effective exposure (sum of old bins within 
	// new bin). The dead time between old bins, therefore only affects the error and not the rate.

	int nOldBins = oldBins.length;
	int nNewBins = newBins.length;

	//   Initialize variables
	int k = 0;   //  k is the index for the old bins
	double leftEdgeOfOldBin = oldBins[k].getLeftEdge();
	double rightEdgeOfOldBin = oldBins[k].getRightEdge();
	double tstop = oldBins[oldBins.length-1].getRightEdge();
	
	double[] rebinnedRates = new double[nNewBins];
	double[] rebinnedErrors = new double[nNewBins];

	double weight = 0;
	double weightedSumOfRates = 0;
	double sumOfWeights = 0;

	//  Loop through the new bins to resample
	//  i is the index for the new bins
	//  k is the index for the old bins

	int i = 0;  //  i is the index for the new bins
	double rightEdgeOfNewBin = newBins[i].getRightEdge();
	while ( i < nNewBins ) {

	    while ( k < nOldBins-1 && newBins[i].covers(oldBins[k]) ) {
		//   Enter here when the newbin contains at least 1 old bin
		
		weight = 1/oldBins[k].getVariance();
		weightedSumOfRates += weight*oldBins[k].getIntensity();
		sumOfWeights += weight;
		// logger.debug(weight+"\t"+weightedSumOfRates+"\t"+sumOfWeights);

		//   Move to the next old bin and define its edges
		k++;

// 		if ( k < nOldBins ) {
// 		    leftEdgeOfOldBin = oldBins[k].getLeftEdge();
// 		    rightEdgeOfOldBin = oldBins[k].getRightEdge();
// 		}
	    }
	    //   We get out of the while loop when the next old bin is not fully contained within the new bin

	    boolean weAreAtStartOrBeforeBin = newBins[i].getRightEdge() <= oldBins[k].getLeftEdge();
	    if ( weAreAtStartOrBeforeBin ) {

		//   We are either right at or before a bin:

		//   If there is a gap in the old bins, and therefore, the new bin ends before or at the start 
		//   of the next old bin, write out the final rate for the new bin and reset counts to 0

		if ( weightedSumOfRates != 0 ) {

 		    rebinnedRates[i] = weightedSumOfRates/sumOfWeights;
		    rebinnedErrors[i] = 1/Math.sqrt(sumOfWeights);

		    //  Reset values for the next new bin
		    weightedSumOfRates = 0;
		    sumOfWeights = 0;
		}
		else {
		    //  Enter here if we are before first bin
 		    rebinnedRates[i] = Double.NaN;
 		    rebinnedErrors[i] = Double.NaN;
		}

		//   Move to the next new bin
		i++;
		weAreAtStartOrBeforeBin = newBins[i].getRightEdge() <= oldBins[k].getLeftEdge();
		while ( weAreAtStartOrBeforeBin ) {

		    //logger.debug("HEY: weAreAtStartOrBeforeBin");
		    rebinnedRates[i] = Double.NaN;
		    rebinnedErrors[i] = Double.NaN;

		    i++;
		    weAreAtStartOrBeforeBin = newBins[i].getRightEdge() <= oldBins[k].getLeftEdge();
		}
	    }
	    else {
		//   We are within an old bin:

		//   If the new bin ends inside the next old bin, 
		//   then add the counts corresponding to the fraction of the old bin, 
		//   and write out the final rate for the new bin. 
		
		double rateOfPreviousBin=0, rateOfThisBin=0, rateOfNextBin=0;
		double centreOfPreviousBin=0, centreOfThisBin=0, centreOfNextBin=0;
		double varianceOfThisBin = 0;
		double slope=0, slopeMinus=0, slopePlus=0;
		
		double deltaMinus = 0, deltaPlus = 0;
		double binWidthOfThisBin = 0;
		double binCentreMinus = 0, binCentrePlus = 0;
		double rateMinus = 0, ratePlus = 0;
		double varianceMinus = 0, variancePlus = 0;

		boolean newBinIsContainedWithinThisOldBin = oldBins[k].covers(newBins[i]);
		boolean weAreWithinTheLastOldBin = (k == nOldBins-1);

		if ( weAreWithinTheLastOldBin ) {

		    rightEdgeOfNewBin = Math.min(rightEdgeOfNewBin, tstop);

		    rateOfPreviousBin = rates[k-1];
		    rateOfThisBin = rates[k];
		    varianceOfThisBin = Math.pow(errors[k], 2);

		    centreOfPreviousBin = oldBinCentres[k-1];
		    centreOfThisBin = oldBinCentres[k];
		    binWidthOfThisBin = oldBinWidths[k];

		    while ( i < nNewBins && newBinIsContainedWithinThisOldBin ) {

			//  The code enters here if the new bin is smaller than the old bin

			DataBin[] splitBins = oldBins[k].splitLastBin(rightEdgeOfNewBin, oldBins[k-1]);
			rateMinus = splitBins[0].getIntensity();
			ratePlus = splitBins[1].getIntensity();
			varianceMinus = splitBins[0].getVariance(); 
			variancePlus = splitBins[1].getVariance();

			//  Define the rebinned rate and error for this bin
			weight = 1/varianceMinus;
			weightedSumOfRates += weight*rateMinus;
			sumOfWeights += weight;

			rebinnedRates[i] = weightedSumOfRates/sumOfWeights;
			rebinnedErrors[i] = 1/Math.sqrt(sumOfWeights);
			

			//  Re-initialize values for the next new bin
			weightedSumOfRates = 0;
			sumOfWeights = 0;

			rateOfPreviousBin = rateMinus;
			rateOfThisBin = ratePlus;
			varianceOfThisBin = variancePlus;
			binWidthOfThisBin = deltaPlus;
			centreOfPreviousBin = binCentreMinus;
			centreOfThisBin = binCentrePlus;
			leftEdgeOfOldBin = rightEdgeOfNewBin;

			//  Move to next new bin, but make sure not to go past tstop
			i++;
			if ( i < nNewBins ) {
			    rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
			    newBinIsContainedWithinThisOldBin = oldBins[k].covers(newBins[i]);
			}
		    }
		}
		else {

		    //  We are within an old bin before the last bin

		    if ( k == 0 ) {
			//  We are in the first bin
			splitBins = oldBins[k].splitFirstBin(rightEdgeOfNewBin, oldBins[k+1]);
		    }
		    else {
			//  We are in a bin after first and before last
			splitBins = oldBins[k].split(rightEdgeOfNewBin, oldBins[k-1], oldBins[k+1]);
		    }
		    rateMinus = splitBins[0].getIntensity();
		    ratePlus = splitBins[1].getIntensity();
		    varianceMinus = splitBins[0].getVariance();
		    variancePlus = splitBins[1].getVariance();

		    //   Add the first part of the old bin to the new bin that ends here
		    weight = 1/varianceMinus;
		    sumOfWeights += weight;
		    weightedSumOfRates += weight*rateMinus;
		    rebinnedRates[i] = weightedSumOfRates/sumOfWeights;
		    rebinnedErrors[i] = 1/Math.sqrt(sumOfWeights);

		    //   Re-initialize to take into account the second piece of the old bin
		    weight = 1/variancePlus;
		    sumOfWeights = weight;
		    weightedSumOfRates = weight*ratePlus;

		    rateOfPreviousBin = rateMinus;
		    rateOfThisBin = ratePlus;
		    varianceOfThisBin = variancePlus;
		    binWidthOfThisBin = splitBins[1].getWidth();
		    centreOfPreviousBin = splitBins[0].getCentre();
		    centreOfThisBin = splitBins[1].getCentre();

		    leftEdgeOfOldBin = rightEdgeOfNewBin;

		    //  Move to next new bin, but make sure not to go past tstop
		    i++;
		    try {
			rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
		    }
		    catch ( ArrayIndexOutOfBoundsException e ) {}
		    newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;

		    while ( i < nNewBins && newBinIsContainedWithinThisOldBin ) {

			weight=0;
			weightedSumOfRates=0;
			sumOfWeights=0;

			//  Work out the rates for each part of the bin
			deltaMinusAndPlus = calcDeltaMinusAndPlus(leftEdgeOfOldBin, rightEdgeOfOldBin, rightEdgeOfNewBin); 
			deltaMinus = deltaMinusAndPlus[0];
			deltaPlus = deltaMinusAndPlus[1];
			
			binCentreMinusAndPlus = calcBinCentreMinusAndPlus(centreOfThisBin, binWidthOfThisBin, deltaMinusAndPlus);
			binCentreMinus = binCentreMinusAndPlus[0];
			binCentrePlus = binCentreMinusAndPlus[1];
			
			rateAndCentreOfThisBin = new double[] {rateOfThisBin, centreOfThisBin};
			double[] rateAndCentreOfPreviousBin = new double[] {rateOfPreviousBin, centreOfPreviousBin};
			double[] rateAndCentreOfNextBin = new double[] {rateOfNextBin, centreOfNextBin};
			slopeMinusAndPlus = calcSlopeMinusAndPlus(rateAndCentreOfThisBin, rateAndCentreOfPreviousBin, rateAndCentreOfNextBin);
			rateMinusAndPlus = calcRateMinusAndPlus(rateAndCentreOfThisBin, binCentreMinusAndPlus, slopeMinusAndPlus);
			rateMinus = rateMinusAndPlus[0];
			ratePlus = rateMinusAndPlus[1];
			
			// Add Poisson noise to rates
// 			nTot = rateOfThisBin*binWidthOfThisBin;
// 			newRates = addPoissonNoise(nTot, rateMinusAndPlus, deltaMinusAndPlus);
// 			rateMinus = newRates[0];
// 			ratePlus = newRates[1];

			//  Distribute the variance of the original bin
			varianceMinus = varianceOfThisBin*binWidthOfThisBin/deltaMinus;
			variancePlus = varianceOfThisBin*binWidthOfThisBin/deltaPlus;

			//   Add the first part of the old bin
			weight = 1/varianceMinus;
			weightedSumOfRates += weight*rateMinus;
			sumOfWeights += weight;
			rebinnedRates[i] = weightedSumOfRates/sumOfWeights;
			rebinnedErrors[i] = 1/Math.sqrt(sumOfWeights);
			
			//  Re-initialize values for the next new bin
			weight = 1/variancePlus;
			sumOfWeights = weight;
			weightedSumOfRates = weight*ratePlus;

			rateOfPreviousBin = rateMinus;
			rateOfThisBin = ratePlus;
			varianceOfThisBin = variancePlus;
			binWidthOfThisBin = deltaPlus;
			centreOfPreviousBin = binCentreMinus;
			centreOfThisBin = binCentrePlus;
			leftEdgeOfOldBin = rightEdgeOfNewBin;

			//  Move to next new bin, but make sure not to go past tstop
			i++;
			try {
			    rightEdgeOfNewBin = Math.min(newBinEdges[2*i+1], tstop);
			}
			catch ( ArrayIndexOutOfBoundsException e ) {}
			newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;
		    }

		    //   Move to the next old bin and define its edges
		    k++;
		    if ( k < nOldBins ) {
			leftEdgeOfOldBin = oldBinEdges[2*k];
			rightEdgeOfOldBin = oldBinEdges[2*k+1];
		    }
		    newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;

		}
	    }

	}

	return new double[][] {rebinnedRates, rebinnedErrors, newBinEdges};
    }

    public static double[] resample(double[] counts, double[] oldBinEdges, double[] newBinEdges) throws BinningException {

	double[] errors = new double[counts.length];
	for ( int i=0; i < counts.length; i++ ) {
	    errors[i] = 1.0;
	}

	double[][] resampledCounts = resample(counts, errors, oldBinEdges, newBinEdges);
	return resampledCounts[0];
    }
    
    public static double[][] resample(double[] rates, double[] errors, double[] oldBinEdges, double newBinWidth) throws BinningException {
    
	double xmin = oldBinEdges[0];
	double xmax = oldBinEdges[oldBinEdges.length-1];
	double[] newBinEdges = BinningUtils.getBinEdges(xmin, xmax, newBinWidth);
	
	return resample(rates, errors, oldBinEdges, newBinEdges);
    }

    private static double[] calcDeltaMinusAndPlus(double leftEdgeOfOldBin, double rightEdgeOfOldBin, double rightEdgeOfNewBin) {
	
	double deltaMinus = rightEdgeOfNewBin - leftEdgeOfOldBin;
	double deltaPlus = rightEdgeOfOldBin - rightEdgeOfNewBin;
	return new double[] {deltaMinus, deltaPlus};
	
    }

    private static double[] calcBinCentreMinusAndPlus(double centreOfThisBin, double binWidthOfThisBin, double[] deltaMinusAndPlus) {

	double deltaMinus = deltaMinusAndPlus[0];
	double deltaPlus = deltaMinusAndPlus[1];
	double binCentreMinus = centreOfThisBin - 0.5*(binWidthOfThisBin - deltaMinus);
	double binCentrePlus = centreOfThisBin + 0.5*(binWidthOfThisBin - deltaPlus);
	return new double[] {binCentreMinus, binCentrePlus};
    }

    private static double[] calcSlopeMinusAndPlus(double[] rateAndCentreOfThisBin, double[] rateAndCentreOfPreviousBin, double[] rateAndCentreOfNextBin) {

	double rateOfThisBin = rateAndCentreOfThisBin[0];
	double centreOfThisBin = rateAndCentreOfThisBin[1];
	double rateOfPreviousBin = rateAndCentreOfPreviousBin[0];
	double centreOfPreviousBin = rateAndCentreOfPreviousBin[1];
	double rateOfNextBin = rateAndCentreOfNextBin[0];
	double centreOfNextBin = rateAndCentreOfNextBin[1];
	double slopeMinus = (rateOfThisBin - rateOfPreviousBin)/(centreOfThisBin - centreOfPreviousBin);
	double slopePlus = (rateOfNextBin - rateOfThisBin)/(centreOfNextBin - centreOfThisBin);
	return new double[] {slopeMinus, slopePlus};
    }

    private static double[] calcRateMinusAndPlus(double[] rateAndCentreOfThisBin, double[] binCentreMinusAndPlus, double[] slopeMinusAndPlus) {

	double rateOfThisBin = rateAndCentreOfThisBin[0];
	double centreOfThisBin = rateAndCentreOfThisBin[1];
	double binCentreMinus = binCentreMinusAndPlus[0];
	double binCentrePlus = binCentreMinusAndPlus[1];
	double slopeMinus = slopeMinusAndPlus[0];
	double slopePlus = slopeMinusAndPlus[1];
	double rateMinus = rateOfThisBin + slopeMinus*(binCentreMinus - centreOfThisBin);
	double ratePlus = rateOfThisBin + slopePlus*(binCentrePlus - centreOfThisBin);
	return new double[] {rateMinus, ratePlus};
    }

    private static double[] addPoissonNoise(double nTot, double[] rateMinusAndPlus, double[] deltaMinusAndPlus) {

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