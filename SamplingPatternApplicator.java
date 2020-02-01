package gb.esac.binner;

import java.util.ArrayList;
import org.apache.log4j.Logger;

import gb.esac.binner.AbsoluteQuantityBin;
import gb.esac.binner.DensityBin;

/**

The class <code>IntensityBinResampler</code> is the new version of <code>Resampler</code>
that makes use of the new set of classes developped to handle intensity bins.

 @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
 @created January 2020 (based on IntensityBinResampler)
 @version January 2020

**/


public final class SamplingPatternApplicator {

    private static Logger logger  = Logger.getLogger(SamplingPatternApplicator.class);

    public static DensityBin[] applySamplingPattern(DensityBin[] intensityBins, AbsoluteQuantityBin[] samplingBins) throws BinningException {
    	ArrayList<DensityBin> newIntensityBinsList = new ArrayList<DensityBin>();
    	int i = 0; // intensity bins
    	int k = 0; // sampling bins
		// Loop through sampling bins one by one
		while (k < samplingBins.length && i < intensityBins.length) {
    		if (samplingBins[k].getRightEdge() == intensityBins[i].getRightEdge()) {
				double leftEdge = intensityBins[i].getLeftEdge();
    			double rightEdge = intensityBins[i].getRightEdge();
    			double value = intensityBins[i].getValue() * samplingBins[k].getValue();
    			double error = intensityBins[i].getError() * samplingBins[k].getValue();
    			newIntensityBinsList.add(new DensityBin(leftEdge,rightEdge,value,error));
    			i++; // move to next intensity bin
	    		k++; // move to next sampling bin
    		}
    		else if (samplingBins[k].getRightEdge() > intensityBins[i].getRightEdge()) {
				while (i < intensityBins.length && (samplingBins[k].getRightEdge() > intensityBins[i].getRightEdge())) {
					double leftEdge = intensityBins[i].getLeftEdge();
	    			double rightEdge = intensityBins[i].getRightEdge();
	    			double value = intensityBins[i].getValue() * samplingBins[k].getValue();
	    			double error = intensityBins[i].getError() * samplingBins[k].getValue();
	    			newIntensityBinsList.add(new DensityBin(leftEdge,rightEdge,value,error));
	    			i++; // move to next intensity bin
	    		} // get out when current sampling bin does not contain the next intensity bin
    		}
    		else {
    			DensityBin[] splitIntensityBin = new DensityBin[2];
	   			if (i == 0) {
	    			splitIntensityBin = intensityBins[i].splitFirstBin(samplingBins[k].getRightEdge(),intensityBins[i+1]);
    			}
    			else if (i == intensityBins.length-1) {
	    			splitIntensityBin = intensityBins[i].splitLastBin(samplingBins[k].getRightEdge(),intensityBins[i-1]);
	    		}
	    		else {
					splitIntensityBin = intensityBins[i].split(samplingBins[k].getRightEdge(),intensityBins[i-1],intensityBins[i+1]);
	    		}
    			double leftEdge = splitIntensityBin[0].getLeftEdge();
    			double rightEdge = splitIntensityBin[0].getRightEdge();
    			double value = splitIntensityBin[0].getValue() * samplingBins[k].getValue();
    			double error = splitIntensityBin[0].getError() * samplingBins[k].getValue();
  				//System.out.println("splitBins[0] = "+leftEdge+" "+rightEdge+"  "+value+" "+error);
    			newIntensityBinsList.add(new DensityBin(leftEdge,rightEdge,value,error));
    			leftEdge = splitIntensityBin[1].getLeftEdge();
    			rightEdge = splitIntensityBin[1].getRightEdge();
    			value = splitIntensityBin[1].getValue() * samplingBins[k+1].getValue();
    			error = splitIntensityBin[1].getError() * samplingBins[k+1].getValue();
    			//System.out.println("splitBins[1] = "+leftEdge+" "+rightEdge+"  "+value+" "+error);
    			newIntensityBinsList.add(new DensityBin(leftEdge,rightEdge,value,error));
    			i++; // move to next intensity bin
	    		k++; // move to next sampling bin
    		}
    	}
    	newIntensityBinsList.trimToSize();
    	int nObservedBins = newIntensityBinsList.size();
    	DensityBin[] observedIntensityBins = new DensityBin[nObservedBins];
    	for (int j = 0; j < nObservedBins; j++ ) {
    		observedIntensityBins[j] = (DensityBin) newIntensityBinsList.get(j);
    	}
    	return observedIntensityBins;
    }




  //   public static DensityBin[] applySamplingPattern(DensityBin[] intensityBins, IBin[] samplingBins) throws BinningException {

		// // The rate in each new bin is estimated using the effective exposure (sum of old bins within new bin).
		// // The dead time between old bins, therefore only affects the error; not the rate.

		// logger.info("Applying sampling pattern:");
		// int nIntensityBins = intensityBins.length;
		// int nSamplingBins = samplingBins.length;
		// double tstop = intensityBins[intensityBins.length-1].getRightEdge();
		// ArrayList<DensityBins> newIntensityBinsList = new ArrayList<DensityBins>();
		// //DensityBin[] newDensityBins = new DensityBin[nSamplingBins]; // These are what we are defining
		// int i = 0; // index for sampling bins
		// int k = 0; // index for intensity bins
		// logger.info("  nIntensityBins = "+nIntensityBins);
		// logger.info("  nSamplingBins = "+nSamplingBins);
		// double rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), intensityBins[k].getRightEdge());
		// double newIntensityValue = samplingBins[i].getValue()*intensityBins[k].getValue();
		// double newErrorValue = samplingBins[i].getError()*intensityBins[k].getValue();
		// DensityBin currentNewIntensityBin = new DensityBin(intensityBins[k].getLeftEdge(), rightEdgeOfNewBin, newIntensityValue, newErrorValue);

		// //  Loop through the new bins and resample
		// while (i < nSamplingBins) {
		    
		//     while (samplingBins[i].contains(intensityBins[k]) && k < nIntensityBins-1) {
		// 		// Enter when the newbin contains at least 1 old bin
		// 		logger.info("  samplingBins["+i+"] contains intensityBins["+k+"]");
		// 		rightEdgeOfNewBin = intensityBins[k].getRightEdge();
		// 		newIntensityValue = samplingBins[i].getValue()*intensityBins[k].getValue(); // sampling value is 0 or 1
		// 		newErrorValue = samplingBins[i].getValue()*intensityBins[k].getValue();				
		// 		currentNewIntensityBin = new DensityBin(intensityBins[k].getLeftEdge(), rightEdgeOfNewBin, newIntensityValue, newErrorValue);
		// 		newIntensityBinsList.add(currentNewIntensityBin);
		// 		k++;
		//     } // Exit loop when next intensity bin is not fully contained within sampling bin
		    
		//     boolean weAreAtStartOrBeforeNextIntensityBin = rightEdgeOfNewBin <= intensityBins[k].getLeftEdge();
		//     if (weAreAtStartOrBeforeNextIntensityBin) {
		// 		logger.info("  We are before or at start of the next intensity bin: defining new intensity.");
		// 		// Close the new intensity bin
		// 		newIntensityBinsList.add(currentNewIntensityBin);
		// 		// Move to next sampling bin 
		// 		i++;
		// 		// And fill with NaNs until we reach next intensity bin
		// 		weAreAtStartOrBeforeNextIntensityBin = rightEdgeOfNewBin <= intensityBins[k].getLeftEdge();
		// 		while (weAreAtStartOrBeforeNextIntensityBin) {
		// 		    logger.info("  We are before next intensity bin: defining new intensity as NaN");
		// 		    newDensityBinsList.add(new DensityBin(intensityBins[k].getLeftEdge(), intensityBins[k].getRightEdge(), Double.NaN, Double.NaN));
		// 		    i++;
		// 		    rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), tstop);
		// 		    weAreAtStartOrBeforeNextIntensityBin = rightEdgeOfNewBin <= intensityBins[k].getLeftEdge();
		// 		} // Exit loop when we have entered the next old bin
		//     }
		//     else { // Enter when we are within an old bin
		// 		DensityBin[] splitBins = new DensityBin[2];
		//         rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), tstop); // Make sure right edge of new bin is less than tstop
		// 		boolean newBinIsContainedWithinThisOldBin = intensityBins[k].contains(rightEdgeOfNewBin);
		// 		boolean weAreWithinTheLastOldBin = (k == nIntensityBins-1);

		// 		if (weAreWithinTheLastOldBin) {
		// 		    logger.info("  We are within the _last_ old bin");
		// 		    DensityBin thisOldBin = new DensityBin(intensityBins[k]);
		// 		    DensityBin previousOldBin = new DensityBin(intensityBins[k-1]);

		// 		    while (newBinIsContainedWithinThisOldBin && i < nSamplingBins) { //  Enter when the new bin is smaller than the old bin
		// 			logger.info("  New bin is contained within this old bin: need to split and join first piece to current new intensity.");

		// 			// split
		// 			splitBins = thisOldBin.splitLastBin(rightEdgeOfNewBin, previousOldBin);

		// 			// join two pieces to close this new bin
		// 			newDensityBins[i] = currentNewIntensityBin.joinWith(splitBins[0]);
		// 			previousOldBin = new DensityBin(splitBins[0]);
		// 			thisOldBin = new DensityBin(splitBins[1]);

		// 			//  Move to next new bin
		// 			i++;
		// 			if (i < nSamplingBins) {
		// 			    rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), tstop);
		// 			    logger.info("  Initialise next new bin with zero intensity");
		// 			    currentNewIntensityBin = new DensityBin(samplingBins[i].getLeftEdge(), rightEdgeOfNewBin, 0,0);
		// 			    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
		// 			}
				
		// 	    } // Exit when new bin is no longer contained within the last old bin
		// 	} 

		// 	else { // We are within an old bin _before_ the last

		// 	    if (k == 0) { // We are in the first bin
		// 			logger.info("  We are within the _first_ old bin, and need to split it.");
		// 			splitBins = intensityBins[k].splitFirstBin(rightEdgeOfNewBin, intensityBins[k+1]);
		// 	    }

		// 	    else { // We are in a bin after first and before last
		// 			logger.info("  We are in a bin after the first and before the last, and need to split it.");
		// 			splitBins = intensityBins[k].split(rightEdgeOfNewBin, intensityBins[k-1], intensityBins[k+1]);
		// 	    }

		// 	    // Close new intensity bin with first part; add second part to currentNewBin; redefine bin variables
		// 	    newDensityBins[i] = currentNewIntensityBin.joinWith(splitBins[0]);
		// 	    currentNewIntensityBin = new DensityBin(splitBins[1]);
		// 	    DensityBin thisOldBin = new DensityBin(splitBins[1]);
		// 	    DensityBin previousOldBin = new DensityBin(intensityBins[k-1]);
		// 	    DensityBin nextOldBin = new DensityBin(intensityBins[k+1]);

		// 	    //  Move to next new bin
		// 	    i++;
		// 	    if (i < nSamplingBins) {
		// 			rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), tstop);
		// 			newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
		// 	    }
			    
		// 	    while (newBinIsContainedWithinThisOldBin && i < nSamplingBins) {
		// 			logger.info("  New bin is contained within this old bin: ");
		// 			logger.info("    we need to split it, join first piece to current intensity to define new intensity,");
		// 			logger.info("    and initialise current intensity with second piece.");
		// 			splitBins = thisOldBin.split(rightEdgeOfNewBin, previousOldBin, nextOldBin);
		// 			newDensityBins[i] = currentNewIntensityBin.joinWith(splitBins[0]);
		// 			currentNewIntensityBin = new DensityBin(splitBins[1]);
		// 			thisOldBin = new DensityBin(splitBins[1]);
		// 			//  Move to next new bin
		// 			i++;
		// 			if (i < nSamplingBins) {
		// 			    rightEdgeOfNewBin = Math.min(samplingBins[i].getRightEdge(), tstop);
		// 			    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);
		// 			}
		// 	    } // Exit when new bin is no longer contained within the old bin

			    
		// 	    //   Move to the next old bin and define its edges
		// 	    k++;
		// 	    if ( k < nIntensityBins ) {
		// 			logger.info("  Moving to the next old bin (k="+k+")");
		// 			thisOldBin = new DensityBin(intensityBins[k]);
		// 	    }
		// 	    newBinIsContainedWithinThisOldBin = thisOldBin.contains(rightEdgeOfNewBin);

		// 	} // End of else block - When we are within an old bin after first and before last 
			
		//     } // End of else block - When we are within an old bin
		    
		// } // End of while (i < nSamplingBins)
		
		// return newDensityBins;
  //   }  

}
