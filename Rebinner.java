package gb.esac.binner;

import cern.colt.list.DoubleArrayList;
import gb.esac.tools.Converter;
import java.util.ArrayList;
import org.apache.log4j.Logger;


public class Rebinner {

    private static Logger logger = Logger.getLogger(Rebinner.class);

    public static double[][] rebinToMinSignif(double[] rates, double[] errors, double[] binEdges, double minSignif) {

	DoubleArrayList listOfNewBinEdges = new DoubleArrayList();
	listOfNewBinEdges.add(binEdges[0]);
	DoubleArrayList listOfNewRates = new DoubleArrayList();
	DoubleArrayList listOfNewErrors = new DoubleArrayList();

	int i=0;
	int nBins = rates.length;
	while ( i < nBins ) {

	    double weight = Math.pow(errors[i], -2);
	    double weightedSumOfRates = weight*rates[i];
	    double sumOfWeights = weight;;
	    double signif = weightedSumOfRates/Math.sqrt(sumOfWeights);

	    while ( i < nBins && signif < minSignif ) {
		
		i++;
		try {
		    weight = Math.pow(errors[i], -2);
		    weightedSumOfRates += weight*rates[i];
		    sumOfWeights += weight;
		    signif = weightedSumOfRates/Math.sqrt(sumOfWeights);
		}
		catch ( ArrayIndexOutOfBoundsException e ) {}
		
	    }
	    double newRate = weightedSumOfRates/sumOfWeights;
	    double newError = 1/Math.sqrt(sumOfWeights);
	    listOfNewRates.add(newRate);
	    listOfNewErrors.add(newError);
	    try {
		listOfNewBinEdges.add(binEdges[2*i+1]);
	    }
	    catch ( ArrayIndexOutOfBoundsException e ) {
		listOfNewBinEdges.add(binEdges[binEdges.length-1]);
	    }
	    i++;
	    if ( i < nBins ) {
		listOfNewBinEdges.add(binEdges[2*i]);
	    }
	}
 	listOfNewBinEdges.trimToSize();
 	listOfNewRates.trimToSize();
 	listOfNewErrors.trimToSize();


	//  Remove last bin
	listOfNewBinEdges.remove(listOfNewBinEdges.size()-1);
	listOfNewBinEdges.remove(listOfNewBinEdges.size()-2);
	listOfNewRates.remove(listOfNewRates.size()-1);
	listOfNewErrors.setSize(listOfNewErrors.size()-1);
	

	//  Define new arrays
	double[] newBinEdges = listOfNewBinEdges.elements();
	double[] newRates = listOfNewRates.elements();
	double[] newErrors = listOfNewErrors.elements();

	return new double[][] {newRates, newErrors, newBinEdges};
    }


    public static double[][] rebinAsPapadakis(double[] freq, double[] pow, int groupSize) {

	//  Algorithm is from Papadakis and Lawrence 1993 MNRAS, 261, pp 612-624
	//  Assuming input freqs are with constant linear binning 

	int ntot = freq.length;
	double oldBinWidth = freq[1] - freq[0];
	double xmin = freq[0] - 0.5*oldBinWidth;
	ArrayList avePow = new ArrayList();
	ArrayList aveFreq = new ArrayList();
	int i = 0;
	while ( i <= (ntot-groupSize) ) {
	    int j = 0;
	    double sumOfLogPow = 0;
	    double prodOfFreq = 1;
	    while ( j < groupSize && i < ntot ) {
		sumOfLogPow += Math.log10(pow[i]) + 0.253;
		prodOfFreq *= freq[i];
		i++;
		j++;
	    }
	    double expo = 1.0/ (double) j;
	    avePow.add(Math.pow(10, sumOfLogPow*expo));
	    aveFreq.add(Math.pow(prodOfFreq, expo));
	}
	
	int nNewBins = avePow.size();
	double papadakisPow = 0;
	double papadakisFreq = 0;
	double[] adjacentBinEdges = new double[nNewBins+1];
	adjacentBinEdges[0] = xmin;
	double halfBinWidth = 0;
	double newBinWidth = 0;
	double[][] papadakis = new double[3][nNewBins];
	for ( i = 0; i < nNewBins; i++ ) {
	    papadakisPow = ((Double) avePow.get(i)).doubleValue();
	    papadakisFreq = ((Double) aveFreq.get(i)).doubleValue();
	    halfBinWidth = papadakisFreq - adjacentBinEdges[i];
	    adjacentBinEdges[i+1] = adjacentBinEdges[i] + 2*halfBinWidth;
	    newBinWidth = 2*halfBinWidth;
	    papadakis[0][i] = papadakisFreq;
	    papadakis[1][i] = papadakisPow;
	    papadakis[2][i] = newBinWidth;
	}
	
	return papadakis;
    }



//     public static double[][] rebinRates(double[] rates, double[] errors, double[] oldBinEdges, double newBinWidth) throws BinningException {

// 	//   Construct adjacent newBinEdges that correspond to newBinTime
// 	double xmin = oldBinEdges[0];
// 	double xmax = oldBinEdges[oldBinEdges.length-1];
// 	double[] newBinEdges = BinningUtils.getLinearBinEdges(xmin, xmax, newBinWidth);

// 	return rebinRates(rates, errors, oldBinEdges, newBinEdges);
//     }


//     public static double[][] rebinRates(double[] rates, double[] errors, double[] oldBinEdges, double[] newBinEdges) throws BinningException {


// 	//    The rate in each new bin is estimated using the effective exposure (sum of old bins within 
// 	//    new bin). The dead time between old bins, therefore only affects the error and not the rate.
	
// 	//  Make sure that nNewBins <= nOldBins
// 	int nOldBins = rates.length;
// 	int nNewBins = (int) newBinEdges.length/2;
// 	boolean cannotRebin = nNewBins > nOldBins;
// 	if ( cannotRebin ) {
// 	    throw new BinningException("\nnNewBins ("+nNewBins+") must be <= nOldBins ("+nOldBins+")");
// 	}


// 	//   Initialize variables
// 	double tstop = oldBinEdges[oldBinEdges.length-1];
// 	double[] rebinnedRates = new double[nNewBins];
// 	double[] rebinnedErrors = new double[nNewBins];

// 	double counts = 0;
// 	double errorCounts = 0;
// 	double exposure = 0;
// 	double effNewBinTime = 0;

// 	int k = 0;   //  k is the index for the old bins
// 	double leftEdge = oldBinEdges[2*k];
// 	double rightEdge = oldBinEdges[2*k+1];

//  	for ( int i=0; i < nNewBins; i++ ) {

// 	    double rightEdgeOfNewBin = newBinEdges[2*i+1];

// 	    //   Sum the counts of the old bins while within the new bin

// 	    while ( k < nOldBins-1 && rightEdge <= rightEdgeOfNewBin ) {

// 		exposure = (rightEdge - leftEdge);
// 		effNewBinTime += exposure;
// 		counts += exposure*rates[k];
// 		errorCounts += Math.pow(exposure*errors[k], 2);

// 		//   Move to the next old bin and define its edges
// 		k++;

// 		if ( k < nOldBins ) {
// 		    leftEdge = oldBinEdges[2*k];
// 		    rightEdge = oldBinEdges[2*k+1];
// 		}

// 	    }
// 	    //   At this point, the next old bin is not fully contained within the new bin
		 

// 	    //   If there is a gap in the old bins, and therefore, the new bin ends before or at the start 
// 	    //   of the next old bin, write out the final rate for the new bin and reset counts to 0

// 	    if ( rightEdgeOfNewBin <= leftEdge ) {

// 		rebinnedRates[i] = counts/effNewBinTime;
// 		rebinnedErrors[i] = Math.sqrt(errorCounts)/effNewBinTime;

// 		////logger.debug("effNewBinTime="+effNewBinTime);
// 		////logger.debug("r="+rebinnedRates[0][i]+"\t\t e="+rebinnedRates[1][i]);

// 		effNewBinTime = 0;
// 		exposure = 0;
// 		counts = 0;
// 		errorCounts = 0;
// 	    }
	    
	    
// 	    //   If the new bin ends inside the next old bin, add the counts corresponding to the 
// 	    //   fraction of the old bin, and write out the final rate for the new bin. 

// 	    //   Here we reset the counts to the other fraction of the old bin. 

// 	    else {

// 		if ( k == nOldBins-1 ) {

// 		    rightEdgeOfNewBin = Math.min(rightEdgeOfNewBin, tstop);
		
// 		    //   Add last bit of counts from the first part of the old bin

// 		    exposure = (rightEdgeOfNewBin - leftEdge);
// 		    counts += rates[k]*exposure;
// 		    errorCounts += Math.pow(errors[k]*exposure, 2);
// 		    effNewBinTime += exposure;
// 		    rebinnedRates[i] = counts/effNewBinTime;
// 		    rebinnedErrors[i] = Math.sqrt(errorCounts)/effNewBinTime;

// 		    ////logger.debug("effNewBinTime="+effNewBinTime);
// 		    ////logger.debug("r="+rebinnedRates[0][i]+"\t\t e="+rebinnedRates[1][i]);
// 		}

// 		else {

// 		    //   Add last bit of counts from the first part of the old bin

// 		    exposure = (rightEdgeOfNewBin - leftEdge);
// 		    counts += rates[k]*exposure;
// 		    errorCounts += Math.pow(errors[k]*exposure, 2);
// 		    effNewBinTime += exposure;
// 		    rebinnedRates[i] = counts/effNewBinTime;
// 		    rebinnedErrors[i] = Math.sqrt(errorCounts)/effNewBinTime;

// 		    ////logger.debug("effNewBinTime="+effNewBinTime);
// 		    ////logger.debug("r="+rebinnedRates[0][i]+"\t\t e="+rebinnedRates[1][i]);


// 		    //   Reset to take into account the second piece of the old bin

// 		    exposure = rightEdgeOfNewBin - leftEdge;
// 		    counts = exposure*rates[k];
// 		    errorCounts = Math.pow(errors[k]*exposure, 2);
// 		    effNewBinTime = exposure;


// 		    //   Move to the next old bin and define its edges

// 		    k++;
// 		    if ( k < nOldBins ) {
// 			leftEdge = oldBinEdges[2*k];
// 			rightEdge = oldBinEdges[2*k+1];
// 		    }

// 		}
// 	    }

// 	}

// 	return new double[][] {rebinnedRates, rebinnedErrors, newBinEdges};

//     }



     public static double[] rebinRates(double[] rates, double[] oldBinEdges, double newBinTime, boolean directSum) throws BinningException {


// 	if ( directSum == false ) {

// 	    return Resampler.resampleRates(rates, oldBinEdges, newBinTime);
// 	}
// 	else {

	    //  Define arrays
	    int nbins = oldBinEdges.length/2;
	    double duration = oldBinEdges[oldBinEdges.length-1] - oldBinEdges[0];
	    int nNewBins = (int) Math.round(duration/newBinTime);
	    newBinTime = duration/nNewBins;
	    double[][] rebinnedRates = new double[2][nNewBins];   //  Contains rates and errors


	    //   Rebin
	    double sumOfCounts = 0;
	    int k = 0;
	    double time = 0;
	    double rightEdgeOfNewBin = 0;
	    double rate = 0;
	    double error = 0;
	    double frac = 0;
	    double effExposure = 0;
	    double sumOfBins = 0;
	    double binSize = 0;

	    for ( int i=0; i < nbins; i++ ) {

		//System.out.println("k = "+k);
		
		time = oldBinEdges[2*i+1];
		rightEdgeOfNewBin = (k+1)*newBinTime;
		rate = rates[i];
		double diff = rightEdgeOfNewBin - time;
		
		if ( time < rightEdgeOfNewBin && diff > 1e-6 ) {

		    //  Calculate the size of the bin from edges 
		    binSize = (oldBinEdges[2*i+1] - oldBinEdges[2*i]);

		    //  Sum the counts  
		    sumOfCounts += rate*binSize;

		}
		else if ( diff > 0 && diff < 1e-6 ) {

		    //  Calculate the size of each bin from edges 
		    binSize = (oldBinEdges[2*i+1] - oldBinEdges[2*i]);

		    //  Sum the counts  
		    sumOfCounts += rate*binSize;

		    //  Calculate the rebinned rates  
		    rebinnedRates[0][k] = sumOfCounts/newBinTime;
		    rebinnedRates[1][k] = Math.sqrt(sumOfCounts)/newBinTime;

		    // Re-initialise sumOfRates and simOfBins  
		    sumOfCounts = 0;
		    
		    // Increment to the next new bin  
		    k++;
		}
		else  {

		    //  Add the last fraction of the rates  
		    frac = rightEdgeOfNewBin - oldBinEdges[2*i];
		    sumOfCounts += rate*frac;

		    //  Calculate the rebinned rates  
		    rebinnedRates[0][k] = sumOfCounts/newBinTime;
		    rebinnedRates[1][k] = Math.sqrt(sumOfCounts)/newBinTime;

		    
		    // Re-initialise sumOfRates and sumOfBins  
		    sumOfCounts = rate*(oldBinEdges[2*i+1] - rightEdgeOfNewBin);
		    sumOfBins = (oldBinEdges[2*i+1] - rightEdgeOfNewBin);
		    
		    // Increment to the next new bin  
		    k++;
		}
	    }
	    
	    return rebinnedRates[0];
	    
// 	}

     }



}