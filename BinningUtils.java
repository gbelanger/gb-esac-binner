package gb.esac.binner;

import cern.colt.list.DoubleArrayList;
import java.util.ArrayList;
import org.apache.log4j.Logger;


public final class BinningUtils {

    private static Logger logger = Logger.getLogger(BinningUtils.class);

    public static double[] getBinEdges(double xmin, double xmax, double binWidth) throws BinningException {
	if ( xmin >= xmax ) {
	    throw new BinningException("Cannot construct bin edges: xmax >= xmin");
	}
	double range = xmax - xmin;
	if ( range < binWidth ) {
	    throw new BinningException("Cannot constuct bin edges: (xmax-xmin) < binWidth.");
	}
	//  Construct bin edges
	DoubleArrayList edges = new DoubleArrayList();
	double leftEdge = xmin - 0.5*Math.ulp(xmin);  //  move back the first left edge by a ULP
	double rightEdge = leftEdge + binWidth;
	while ( rightEdge <= (xmax+binWidth/2) ) {
	    edges.add(leftEdge);
	    edges.add(rightEdge);
	    leftEdge = rightEdge;
	    rightEdge = leftEdge + binWidth;
	}
	edges.trimToSize();
	double lastEdge = edges.get(edges.size()-1);
	lastEdge += 0.5*Math.ulp(lastEdge);  //  Push forward the last right edge by a ULP
	edges.set(edges.size()-1, lastEdge);
	return edges.elements();
    }

    public static double[] getBinEdges(double xmin, double xmax, int nbins) throws BinningException {
	double binWidth = (xmax-xmin)/nbins;
	return getBinEdges(xmin, xmax, binWidth);
    }

    public static double[] getBinEdges(double binWidth, int nbins) throws BinningException {
	double xmin = 0;
	double xmax = binWidth*nbins;
	return getBinEdges(xmin, xmax, binWidth);
    }

    public static double[] getBinEdges(double xmin, double[] binWidths) throws BinningException {
	//  Construct bin edges
	DoubleArrayList edges = new DoubleArrayList();
	double leftEdge = xmin - 0.5*Math.ulp(xmin);  //  move back the first left edge by 0.5 ULP
	double rightEdge = leftEdge + binWidths[0];
	for ( int i=1; i < binWidths.length; i++ ) {
	    edges.add(leftEdge);
	    edges.add(rightEdge);
	    leftEdge = rightEdge;
	    rightEdge = leftEdge + binWidths[i];
	}
	edges.trimToSize();
	double lastEdge = edges.get(edges.size()-1);
	lastEdge += 0.5*Math.ulp(lastEdge);  //  Push forward the last right edge by 0.5 ULP
	edges.set(edges.size()-1, lastEdge);
	return edges.elements();
    }

    public static double[] getBinEdgesFromLeftAndRightEdges(double[] leftEdges, double[] rightEdges) throws BinningException {
	if ( leftEdges.length != rightEdges.length) {
	    throw new BinningException("Left and right bin edges array lengths not equal");
	}
	int nBinEdges = leftEdges.length + rightEdges.length;
	double[] binEdges = new double[nBinEdges];
	for ( int i=0 ; i < leftEdges.length; i++ ) {
	    binEdges[2*i] = leftEdges[i];
	    binEdges[2*i+1] = rightEdges[i];
	}
	return binEdges;
    }
    
    public static double[] getBinEdgesInLogSpace(double xmin, double xmax, int nbins) throws BinningException {
	if ( xmin <= 0 || xmax <= 0 ) {
	    throw new BinningException("Cannot contruct log bin edges: xmin <= 0 or xmax <= 0");
	}
	double logBinWidth = Math.log(xmax/xmin)/nbins;
	double[] singleBinEdgesInLogSpace = new double[nbins+1];
	for ( int i=0; i <= nbins; i++ ) {
	    singleBinEdgesInLogSpace[i] = Math.log(xmin) + logBinWidth*i;
	}
	double[] binEdges = new double[nbins*2];
	binEdges[0] = Math.exp(singleBinEdgesInLogSpace[0]);
	binEdges[1] = Math.exp(singleBinEdgesInLogSpace[1]);
	for ( int i=1; i < nbins; i++ ) {
	    binEdges[2*i] = Math.exp(singleBinEdgesInLogSpace[i]);
	    binEdges[2*i+1] = Math.exp(singleBinEdgesInLogSpace[i+1]);
	}
	return binEdges;
    }

    public static double[] getBinEdgesFromBinCentres(double[] binCentres) throws BinningException {
	int nBins = binCentres.length;
	if ( nBins < 2 ) {
	    throw new BinningException("Cannot construct bin edges: there is only one bin");
	}
	double[] halfBinWidths = new double[nBins];
	double[] binEdges = new double[2*nBins];
	halfBinWidths[0] = (binCentres[1] - binCentres[0])/2;
	binEdges[0] = binCentres[0] - halfBinWidths[0];
	binEdges[1] = binCentres[0] + halfBinWidths[0];
	for ( int i=1; i < nBins; i++ ) {
	    binEdges[2*i] = binEdges[2*(i-1)+1];
	    halfBinWidths[i] = binCentres[i] - binEdges[2*i];
	    binEdges[2*i+1] = binCentres[i] + halfBinWidths[i];
	    logger.debug(i+"\t"+binEdges[2*i]+"\t"+binEdges[2*i+1]+"\t"+halfBinWidths[i]);
	    if ( halfBinWidths[i] < 0 ) {
	    	logger.debug("i="+i+" "+halfBinWidths[i]);
	     	System.exit(-1);
	    }
	}
	return binEdges;
    }

    public static double[] getBinEdgesFromBinCentresAndHalfWidths(double[] binCentres, double[] halfBinWidths) throws BinningException {
	if ( binCentres.length != halfBinWidths.length ) {
	    throw new BinningException("Incompatible array lengths: binCentres.length != halfBinWidths.length");
	}
	int nBins = binCentres.length;
	double[] binEdges = new double[2*nBins];
	for ( int i=0; i < nBins; i++ ) {

	    binEdges[2*i] = binCentres[i] - halfBinWidths[i];
	    binEdges[2*i+1] = binCentres[i] + halfBinWidths[i];	    
	}
	return binEdges;
    }


    public static double[] getBinWidthsFromBinEdges(double[] binEdges) {
	int nBins = binEdges.length/2;
	double[] binWidths = new double[nBins];
	for ( int i=0; i < nBins; i++ ) {
	    binWidths[i] = binEdges[2*i+1] - binEdges[2*i];
	}
	return binWidths;
    }

    public static double[] getHalfBinWidthsFromBinEdges(double[] binEdges) {
	int nBins = binEdges.length/2;
	double[] halfBinWidths = new double[nBins];
	for ( int i=0; i < nBins; i++ ) {
	    halfBinWidths[i] = (binEdges[2*i+1] - binEdges[2*i])/2;
	}
	return halfBinWidths;
    }

    public static double[] getHalfBinWidthsFromBinCentres(double[] binCentres) throws BinningException {
	return getHalfBinWidthsFromBinEdges(getBinEdgesFromBinCentres(binCentres));
    }

    public static double[] getBinCentresFromBinEdges(double[] binEdges) {
	int nBins = binEdges.length/2;
	double[] binCentres = new double[nBins];
	for ( int i=0; i < nBins; i++ ) {
	    binCentres[i] = (binEdges[2*i+1] + binEdges[2*i])/2;
	}
	return binCentres;
    }

    public static double[] getBinCentres(double xmin, double xmax, double binWidth) throws BinningException {
	return getBinCentresFromBinEdges(getBinEdges(xmin, xmax, binWidth));
    }

    public static double[] getBinCentres(double xmin, double xmax, int nBins) throws BinningException {
	return getBinCentresFromBinEdges(getBinEdges(xmin, xmax, nBins));
    }

    public static void checkArrayLengthsAreEqual(double[] array1, double[] array2) throws BinningException {
	if ( array1.length != array2.length ) {
	    throw new BinningException("Number of elements is different: array1.length ("+array1.length+") != array2.length ("+array2.length+")");
	}
    }


}
