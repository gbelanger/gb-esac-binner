package gb.esac.binner;

import java.util.Date;

import cern.jet.random.Poisson;
import cern.jet.random.engine.MersenneTwister64;


/**

   The final class <code>Resampler</code> defines the methods to resample data.
   In contrast to rebinning (regrouping), resampling implies splitting bins.

   @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
   @created August 2010
   @version March 2013

**/


public final class Resampler {

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

    public static double[][] resample(double[] rates, double[] errors, double[] oldBinEdges, double[] newBinEdges) {

	//  The rate in each new bin is estimated using the effective exposure (sum of old bins within 
	// new bin). The dead time between old bins, therefore only affects the error and not the rate.

	int nOldBins = rates.length;
	int nNewBins = (int) newBinEdges.length/2;
	double[] oldBinCentres = new double[nOldBins];
	double[] oldBinWidths = new double[nOldBins];
	for ( int i=0; i < nOldBins; i++ ) {
	    oldBinCentres[i] = 0.5*(oldBinEdges[2*i] + oldBinEdges[2*i+1]);
	    oldBinWidths[i] = oldBinEdges[2*i+1] - oldBinEdges[2*i];
	}

	//   Initialize variables
	int k = 0;   //  k is the index for the old bins
	double leftEdgeOfOldBin = oldBinEdges[2*k];
	double rightEdgeOfOldBin = oldBinEdges[2*k+1];
	double tstop = oldBinEdges[oldBinEdges.length-1];
	//logger.debug("leftOld = "+leftEdgeOfOldBin+"\t rightOld = "+rightEdgeOfOldBin);
	
	double[] rebinnedRates = new double[nNewBins];
	double[] rebinnedErrors = new double[nNewBins];

	double weight = 0;
	double weightedSumOfRates = 0;
	double sumOfWeights = 0;


	//  Loop through the new bins to resample
	//  i is the index for the new bins
	//  k is the index for the old bins

	int i = 0;
	double rightEdgeOfNewBin = newBinEdges[2*i+1];

	//logger.debug("i="+i+" rightNew = "+rightEdgeOfNewBin);

	while ( i < nNewBins ) {

	    boolean oldBinIsContainedInNewBin = rightEdgeOfOldBin <= (rightEdgeOfNewBin + Math.ulp(0)); 

	    //logger.debug("oldBinIsContainedInNewBin = "+oldBinIsContainedInNewBin);

	    while ( k < nOldBins-1 && oldBinIsContainedInNewBin ) {

		//   Enter here when the newbin contains at least 1 old bin
		
		weight = Math.pow(errors[k], -2);
		weightedSumOfRates += weight*rates[k];
		sumOfWeights += weight;

		//logger.debug(weight+"\t"+weightedSumOfRates+"\t"+sumOfWeights);

		//   Move to the next old bin and define its edges
		k++;
		if ( k < nOldBins ) {

		    leftEdgeOfOldBin = oldBinEdges[2*k];
		    rightEdgeOfOldBin = oldBinEdges[2*k+1];

		    //logger.debug("Move to next old bin");
		    //logger.debug("leftOld = "+leftEdgeOfOldBin+"\t rightOld = "+rightEdgeOfOldBin);
		}

		oldBinIsContainedInNewBin = rightEdgeOfOldBin <= (rightEdgeOfNewBin + Math.ulp(0));
		//logger.debug("oldBinIsContainedInNewBin = "+oldBinIsContainedInNewBin);

	    }
	    //   We get out of the while loop when the next old bin is not fully contained within the new bin

	    boolean weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= (leftEdgeOfOldBin + Math.ulp(0));
	    //logger.debug(rightEdgeOfNewBin+" <= "+leftEdgeOfOldBin+" = "+weAreAtStartOrBeforeBin);

	    if ( weAreAtStartOrBeforeBin ) {

		//   We are either right at or before a bin
		//   If there is a gap in the old bins, and therefore, the new bin ends before or at the start 
		//   of the next old bin, write out the final rate for the new bin and reset counts to 0

		if ( weightedSumOfRates != 0 ) {

 		    rebinnedRates[i] = weightedSumOfRates/sumOfWeights;
		    rebinnedErrors[i] = 1/Math.sqrt(sumOfWeights);

		    //logger.debug("RebinnedRates["+i+"] = "+rebinnedRates[i]+" +/- "+rebinnedErrors[i]);

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
		//logger.debug("Move the next new bin");

		rightEdgeOfNewBin = newBinEdges[2*i+1] + Math.ulp(0);
		//logger.debug("i="+i+" rightNew = "+rightEdgeOfNewBin);

		weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= (leftEdgeOfOldBin + Math.ulp(0));
		//logger.debug("weAreAtStartOrBeforeBin = "+weAreAtStartOrBeforeBin);

		while ( weAreAtStartOrBeforeBin ) {

		    //logger.debug("HEY: weAreAtStartOrBeforeBin");
		    rebinnedRates[i] = Double.NaN;
		    rebinnedErrors[i] = Double.NaN;

		    //logger.debug("RebinnedRates["+i+"] = NaN");

		    //logger.debug("Move the next new bin");
		    i++;
		    rightEdgeOfNewBin = newBinEdges[2*i+1];

		    //logger.debug("i="+i+" rightNew = "+rightEdgeOfNewBin);
		    weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= leftEdgeOfOldBin;

		    //logger.debug("weAreAtStartOrBeforeBin = "+weAreAtStartOrBeforeBin);
		}
	    	
	    }
	    else {

		//   We are within an old bin
		//   If the new bin ends inside the next old bin, add the counts corresponding to the 
		//   fraction of the old bin, and write out the final rate for the new bin. 
		
		double rateOfPreviousBin=0, rateOfThisBin=0, rateOfNextBin=0;
		double centreOfPreviousBin=0, centreOfThisBin=0, centreOfNextBin=0;
		double varianceOfThisBin = 0;
		double slope=0, slopeMinus=0, slopePlus=0;
		
		double deltaMinus = 0, deltaPlus = 0;
		double binWidthOfThisBin = 0;
		double binCentreMinus = 0, binCentrePlus = 0;
		double rateMinus = 0, ratePlus = 0;
		double varianceMinus = 0, variancePlus = 0;

		boolean weAreWithinTheLastOldBin = k == nOldBins-1;
		boolean newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;

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

			double[] deltaMinusAndPlus = calcDeltaMinusAndPlus(leftEdgeOfOldBin, rightEdgeOfOldBin, rightEdgeOfNewBin); 
			deltaMinus = deltaMinusAndPlus[0];
			deltaPlus = deltaMinusAndPlus[1];

			double[] binCentreMinusAndPlus = calcBinCentreMinusAndPlus(centreOfThisBin, binWidthOfThisBin, deltaMinusAndPlus);
			binCentreMinus = binCentreMinusAndPlus[0];
			binCentrePlus = binCentreMinusAndPlus[1];

			slope = (rateOfThisBin - rateOfPreviousBin)/(centreOfThisBin - centreOfPreviousBin);
			double[] slopeMinusAndPlus = new double[] {slope, slope};
			double[] rateAndCentreOfThisBin = new double[] {rateOfThisBin, centreOfThisBin};
			double[] rateMinusAndPlus = calcRateMinusAndPlus(rateAndCentreOfThisBin, binCentreMinusAndPlus, slopeMinusAndPlus);
			rateMinus = rateMinusAndPlus[0];
			ratePlus = rateMinusAndPlus[1];

			// Add Poisson noise to rates
// 			double nTot = rateOfThisBin*binWidthOfThisBin;
// 			double[] newRates = addPoissonNoise(nTot, rateMinusAndPlus, deltaMinusAndPlus);
// 			rateMinus = newRates[0];
// 			ratePlus = newRates[1];

			//  Distribute the variance of the original bin
			varianceMinus = varianceOfThisBin*binWidthOfThisBin/deltaMinus;
			variancePlus = varianceOfThisBin*binWidthOfThisBin/deltaPlus;

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
			try {
			    rightEdgeOfNewBin = Math.min(newBinEdges[2*i+1], tstop);
			    newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;
			}
			catch ( ArrayIndexOutOfBoundsException e ) {
			    newBinIsContainedWithinThisOldBin = rightEdgeOfNewBin <= rightEdgeOfOldBin;
			}
			    
		    }
		}

		else {

		    try {
			//  We are within an old bin before the last bin

			rateOfPreviousBin = rates[k-1];
			rateOfThisBin = rates[k];
			rateOfNextBin = rates[k+1];
		    
			centreOfPreviousBin = oldBinCentres[k-1];
			centreOfThisBin = oldBinCentres[k];
			centreOfNextBin = oldBinCentres[k+1];

			double[] rateAndCentreOfThisBin = new double[] {rateOfThisBin, centreOfThisBin};
			double[] rateAndCentreOfPreviousBin = new double[] {rateOfPreviousBin, centreOfPreviousBin};
			double[] rateAndCentreOfNextBin = new double[] {rateOfNextBin, centreOfNextBin};

			double[] slopeMinusAndPlus = calcSlopeMinusAndPlus(rateAndCentreOfThisBin, rateAndCentreOfPreviousBin, rateAndCentreOfNextBin);
			slopeMinus = slopeMinusAndPlus[0];
			slopePlus = slopeMinusAndPlus[1];
		    }

		    catch ( ArrayIndexOutOfBoundsException e ) {

			//  We are in the first bin

			rateOfThisBin = rates[k];
			rateOfNextBin = rates[k+1];
		    
			centreOfThisBin = oldBinCentres[k];
			centreOfNextBin = oldBinCentres[k+1];

			slope = (rateOfNextBin - rateOfThisBin)/(centreOfNextBin - centreOfThisBin);
			slopeMinus = slope;
			slopePlus = slope;
		    }


		    //  Work out the rates for each part of the bin
		    binWidthOfThisBin = oldBinWidths[k];

		    double[] deltaMinusAndPlus = calcDeltaMinusAndPlus(leftEdgeOfOldBin, rightEdgeOfOldBin, rightEdgeOfNewBin); 
		    deltaMinus = deltaMinusAndPlus[0];
		    deltaPlus = deltaMinusAndPlus[1];
		    
		    double[] binCentreMinusAndPlus = calcBinCentreMinusAndPlus(centreOfThisBin, binWidthOfThisBin, deltaMinusAndPlus);
		    binCentreMinus = binCentreMinusAndPlus[0];
		    binCentrePlus = binCentreMinusAndPlus[1];

		    double[] rateAndCentreOfThisBin = new double[] {rateOfThisBin, centreOfThisBin};
		    double[] slopeMinusAndPlus = new double[] {slopeMinus, slopePlus};
		    double[] rateMinusAndPlus = calcRateMinusAndPlus(rateAndCentreOfThisBin, binCentreMinusAndPlus, slopeMinusAndPlus);
		    rateMinus = rateMinusAndPlus[0];
		    ratePlus = rateMinusAndPlus[1];
		    
		    // Add Poisson noise to rates
// 		    double nTot = rateOfThisBin*binWidthOfThisBin;
// 		    double[] newRates = addPoissonNoise(nTot, rateMinusAndPlus, deltaMinusAndPlus);
// 		    rateMinus = newRates[0];
// 		    ratePlus = newRates[1];

		    //  Distribute the variance of the original bin
		    varianceOfThisBin = Math.pow(errors[k], 2);
		    varianceMinus = varianceOfThisBin*binWidthOfThisBin/deltaMinus;
		    variancePlus = varianceOfThisBin*binWidthOfThisBin/deltaPlus;

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
