package gb.esac.binner;


import org.apache.log4j.Logger;

/**

The class <code>IntensityBinResampler</code> is the new version of <code>Resampler</code>
that makes use of the new set of classes developped to handle intensity bins.


 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created March 2013
 @version September 2018

 **/


public final class IntensityBinResampler {

    private static Logger logger  = Logger.getLogger(IntensityBinResampler.class);

    public static DensityBin[] resample(DensityBin[] oldBins, IBin[] newBins) throws BinningException {

	// The rate in each new bin is estimated using the effective exposure (sum of old bins within new bin).
	// The dead time between old bins, therefore only affects the error; not the rate.

	logger.info("Resampling intensity bins:");
	int nOldBins = oldBins.length;
	int nNewBins = newBins.length;
	double tstop = oldBins[oldBins.length-1].getRightEdge();
	DensityBin[] newDensityBins = new DensityBin[nNewBins]; // These are what we are defining
	int i = 0; // index for new bins
	int k = 0; // index for old bins
	logger.info("  nOldBins = "+nOldBins);
	logger.info("  nNewBins = "+nNewBins);
	double rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
	DensityBin currentNewDensityBin = new DensityBin(new Density(0,0),newBins[i]);
	//  Loop through the new bins and resample
	while (i < nNewBins) {
	    
	    while (newBins[i].contains(oldBins[k]) && k < nOldBins-1) {
		// Enter when the newbin contains at least 1 old bin
		logger.info("  newBins["+i+"] contains oldBins["+k+"]");
		currentNewDensityBin = currentNewDensityBin.joinWith(oldBins[k]);
		k++;
	    } // Exit loop when next old bin is not fully contained within new bin
	    
	    boolean weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= oldBins[k].getLeftEdge();
	    if (weAreAtStartOrBeforeBin) {
		logger.info("  We are before or at start of the next old bin: defining new intensity.");
		// Close the new intensity bin
		newDensityBins[i] = new DensityBin(currentNewDensityBin);
		// Move to next new bin and fill with NaNs until we reach next old bin
		i++;
		weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= oldBins[k].getLeftEdge();
		while (weAreAtStartOrBeforeBin) {
		    logger.info("  We are before next old bin: defining new intensity as NaN");
		    newDensityBins[i] = new DensityBin(new Density(Double.NaN,Double.NaN), newBins[i]);
		    i++;
		    rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
		    weAreAtStartOrBeforeBin = rightEdgeOfNewBin <= oldBins[k].getLeftEdge();
		} // Exit loop when we have entered the next old bin
	    }

	    else { // Enter when we are within an old bin
		DensityBin[] splitBins = new DensityBin[2];
	        rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop); // Make sure right edge of new bin is less than tstop
		boolean newBinIsContainedWithinThisOldBin = oldBins[k].contains(rightEdgeOfNewBin);
		boolean weAreWithinTheLastOldBin = (k == nOldBins-1);

		if (weAreWithinTheLastOldBin) {
		    logger.info("  We are within the _last_ old bin");
		    DensityBin thisOldBin = new DensityBin(oldBins[k]);
		    DensityBin previousOldBin = new DensityBin(oldBins[k-1]);

		    while (newBinIsContainedWithinThisOldBin && i < nNewBins) { //  Enter when the new bin is smaller than the old bin
			logger.info("  New bin is contained within this old bin: need to split and join first piece to current new intensity.");

			// split
			splitBins = thisOldBin.splitLastBin(rightEdgeOfNewBin, previousOldBin);

			// join two pieces to close this new bin
			newDensityBins[i] = currentNewDensityBin.joinWith(splitBins[0]);
			previousOldBin = new DensityBin(splitBins[0]);
			thisOldBin = new DensityBin(splitBins[1]);

			//  Move to next new bin
			i++;
			if (i < nNewBins) {
			    rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
			    logger.info("  Initialise next new bin with zero intensity");
			    currentNewDensityBin = new DensityBin(newBins[i].getLeftEdge(), rightEdgeOfNewBin, 0,0);
			    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
			}
			
		    } // Exit when new bin is no longer contained within the last old bin
		} 

		else { // We are within an old bin _before_ the last

		    if (k == 0) { // We are in the first bin
			logger.info("  We are within the _first_ old bin, and need to split it.");
			splitBins = oldBins[k].splitFirstBin(rightEdgeOfNewBin, oldBins[k+1]);
		    }

		    else { // We are in a bin after first and before last
			logger.info("  We are in a bin after the first and before the last, and need to split it.");
			splitBins = oldBins[k].split(rightEdgeOfNewBin, oldBins[k-1], oldBins[k+1]);
		    }

		    // Close new intensity bin with first part; add second part to currentNewBin; redefine bin variables
		    newDensityBins[i] = currentNewDensityBin.joinWith(splitBins[0]);
		    currentNewDensityBin = new DensityBin(splitBins[1]);
		    DensityBin thisOldBin = new DensityBin(splitBins[1]);
		    DensityBin previousOldBin = new DensityBin(oldBins[k-1]);
		    DensityBin nextOldBin = new DensityBin(oldBins[k+1]);

		    //  Move to next new bin
		    i++;
		    if (i < nNewBins) {
			rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
			newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
		    }
		    
		    while (newBinIsContainedWithinThisOldBin && i < nNewBins) {
			logger.info("  New bin is contained within this old bin: ");
			logger.info("    we need to split it, join first piece to current intensity to define new intensity,");
			logger.info("    and initialise current intensity with second piece.");
			splitBins = thisOldBin.split(rightEdgeOfNewBin, previousOldBin, nextOldBin);
			newDensityBins[i] = currentNewDensityBin.joinWith(splitBins[0]);
			currentNewDensityBin = new DensityBin(splitBins[1]);
			thisOldBin = new DensityBin(splitBins[1]);
			//  Move to next new bin
			i++;
			if (i < nNewBins) {
			    rightEdgeOfNewBin = Math.min(newBins[i].getRightEdge(), tstop);
			    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
			}
		    } // Exit when new bin is no longer contained within the old bin

		    
		    //   Move to the next old bin and define its edges
		    k++;
		    if ( k < nOldBins ) {
			logger.info("  Moving to the next old bin (k="+k+")");
			thisOldBin = new DensityBin(oldBins[k]);
		    }
		    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);

		} // End of else block - When we are within an old bin after first and before last 
		
	    } // End of else block - When we are within an old bin
	    
	} // End of while (i < nNewBins)
	
	return newDensityBins;
    }


    
    
    // public static double[] resample(double[] counts, double[] oldBinEdges, double[] newBinEdges) throws BinningException {
    // 	double[] errors = new double[counts.length];
    // 	for ( int i=0; i < counts.length; i++ ) {
    // 	    errors[i] = 1.0;
    // 	}
    // 	double[][] resampledCounts = resample(counts, errors, oldBinEdges, newBinEdges);
    // 	return resampledCounts[0];
    // }
    
    // public static double[][] resample(double[] rates, double[] errors, double[] oldBinEdges, double newBinWidth) throws BinningException {
    // 	double xmin = oldBinEdges[0];
    // 	double xmax = oldBinEdges[oldBinEdges.length-1];
    // 	double[] newBinEdges = BinningUtils.getBinEdges(xmin, xmax, newBinWidth);
    // 	return resample(rates, errors, oldBinEdges, newBinEdges);
    // }

  

}
