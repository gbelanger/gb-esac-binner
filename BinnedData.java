package gb.esac.timeseries;


import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import gb.esac.tools.BasicStats;
import gb.esac.tools.MinMax;
import java.io.IOException;
import java.util.Arrays;
import org.apache.log4j.Logger;
import gb.esac.tools.DataUtils;


/**
 *

The class <code>BinnedData</code> is an abstract class that represents and binned data.

 *
 * @author <a href="mailto:guilaume.belanger@esa.int">Guillaume Belanger</a>
 * @version 1.0 (June 2015, ESAC)
 */
public abstract class BinnedData {
    private static Logger logger  = Logger.getLogger(BinnedData.class);

    private double[] binEdges;
    private double[] leftBinEdges;
    private double[] rightBinEdges;
    private int nBins;
    private double[] binCentres;
    private double[] binWidths;
    private double[] halfBinWidths;
    private boolean binWidthIsConstant = false;
    private double minBinWidth;
    private double maxBinWidth;
    private double avgBinWidth;
    private double binCentreAtMinBinHeight;
    private double binCentreAtMaxBinHeight;

    private double[] gapEdges;
    private double[] gapLengths;
    private int nGaps;
    private int nNonNaNs;
    private double minGap;
    private double maxGap;
    private double meanGap;
    private double sumOfGaps;
    private boolean thereAreGaps = false; 

    private int nSamplingFunctionBins;
    private double[] samplingFunctionValues;
    private double[] samplingFunctionBinEdges;

    private double[] binHeights;
    private double minBinHeight;
    private double maxBinHeight;
    private double meanBinHeight;
    private double sumOfBinHeights;
    private double sumOfSquaredBinHeights;
    private double varianceInBinHeights;
    private double meanDeviationInBinHeights;
    private double skewnessInBinHeights;
    private double kurtosisInBinHeights;

    private void printInfo() {
	logger.info("Bin heights are defined");
	logger.info("  Sum of binHeights = "+this.sumOfBinHeights);
	logger.info("  Mean binHeight = "+this.meanBinHeight);
	logger.info("  Min binHeight = "+this.minBinHeight+" (binCentre="+this.binCentreAtMinBinHeight+")");
	logger.info("  Max binHeight = "+this.maxBinHeight+" (binCentre="+this.binCentreAtMaxBinHeight+")");	
	logger.info("  Variance in bin heights = "+this.varianceInBinHeights);
	logger.info("  Fractional Rms = "+this.fractionalRms);
    }

    private void setBinEdges(double xStart, double[] binEdges) {
	this.binEdges = new double[binEdges.length];
	this.leftBinEdges = new double[binEdges.length/2];
	this.rightBinEdges = new double[binEdges.length/2];
	this.nBins = binEdges.length/2;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binEdges[2*i] = binEdges[2*i];
	    this.binEdges[2*i+1] = binEdges[2*i+1];
	    this.leftBinEdges[i] = binEdges[2*i];
	    this.rightBinEdges[i] = binEdges[2*i+1];
	}
	this.xStart = xStart;
	this.xRange = this.binEdges[this.binEdges.length-1] - this.binEdges[0];
	this.xStop = this.xStart + this.xRange;
	this.xMid = (this.xStart + this.xStop)/2;
	logger.info("Series has "+this.nBins+" bins");
	logger.info("  XStart = "+this.xStart);
    	logger.info("  XStop = "+this.xStop);
	logger.info("  XRange = "+this.xRange);
	this.binCentres = new double[this.nBins];
	this.binWidths = new double[this.nBins];
	this.halfBinWidths = new double[this.nBins];
	double min = Double.MAX_VALUE;
	double max = -Double.MAX_VALUE;
	double avg = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binCentres[i] = (this.binEdges[2*i] + this.binEdges[2*i+1])/2;
	    this.binWidths[i] = this.binEdges[2*i+1] - this.binEdges[2*i];
	    this.halfBinWidths[i] = this.binWidths[i]/2.0;
	    min = Math.min(min, binWidths[i]);
	    max = Math.max(max, binWidths[i]);
	    avg += binWidths[i];
	}
	this.minBinWidth = min;
	this.maxBinWidth = max;
	this.avgBinWidth = avg/this.nBins;

	//  Check if bin width is constant, excluding the last bin
	double[] widths = getBinWidths();
	double[] w = new double[widths.length-1];
	for ( int i=0; i < widths.length-1; i++ ) {
	    w[i] = widths[i];
	}
	double var = BasicStats.getVariance(w);
	if ( var < 1e-10 || Double.isNaN(var) ) {
	    this.binWidthIsConstant = true;
	    logger.info("  Bin width is constant = "+this.binWidths[0]);
	}
	else {
	    this.binWidthIsConstant = false;
	    logger.warn("  Bin width is not constant");
	    logger.info("  Min bin width = "+this.minBinWidth);
	    logger.info("  Max bin width = "+this.maxBinWidth);
	    logger.info("  Average bin width = "+this.avgBinWidth);
	}

	// Define gap edges and sampling function
	this.gapEdges = new double[2*(this.nBins-1)];
	this.gapLengths = new double[this.nBins-1];
	DoubleArrayList samplingFuncValuesList = new DoubleArrayList();
	DoubleArrayList samplingFuncEdgesList = new DoubleArrayList();
	samplingFuncEdgesList.add(this.binEdges[0]);
	samplingFuncEdgesList.add(this.binEdges[1]);
	samplingFuncValuesList.add(1); // never starts with a gap because if there is one it is taken out
	double minGap = Double.MAX_VALUE;
	double maxGap = -Double.MAX_VALUE;
	int nGaps = 0;
	double sumOfGaps = 0;
	for ( int i=1; i < this.nBins; i++ ) {
 	    double gap = this.binEdges[2*i] - this.binEdges[2*i-1];
	    if ( gap > Math.ulp(0d) ) {
		nGaps++;
		sumOfGaps += gap;
		samplingFuncEdgesList.add(this.binEdges[2*i-1]);
		samplingFuncEdgesList.add(this.binEdges[2*i]);
		samplingFuncValuesList.add(0);
	    }
	    samplingFuncEdgesList.add(this.binEdges[2*i]);
	    samplingFuncEdgesList.add(this.binEdges[2*i+1]);
	    samplingFuncValuesList.add(1);
	    
	    minGap = Math.min(minGap, gap);
	    maxGap = Math.max(maxGap, gap);
 	    this.gapLengths[i-1] = gap;
	    this.gapEdges[2*(i-1)] = this.binEdges[2*i-1];
	    this.gapEdges[2*(i-1)+1] = this.binEdges[2*i];
	}
	this.nGaps = nGaps;
	this.sumOfGaps = sumOfGaps;
	this.meanGap = sumOfGaps/nGaps;
	this.maxGap = maxGap;
	this.minGap = minGap;
	if ( this.maxGap > Math.ulp(0d) ) {
	    this.thereAreGaps = true;
	    logger.warn("There are "+nGaps+" gaps in timeline");
	    logger.info("  Total gap time = "+sumOfGaps);
	    logger.info("  Gap fraction wrt xRange = "+(sumOfGaps/this.xRange));
	    logger.info("  Mean gap = "+meanGap);
	    logger.info("  Max gap = "+maxGap);
	}
	else {
	    this.thereAreGaps = false;
	    logger.info("No gaps in timeline");
	}
	samplingFuncValuesList.trimToSize();
	samplingFuncEdgesList.trimToSize();
	this.samplingFunctionValues = samplingFuncValuesList.elements();
	this.samplingFunctionBinEdges = samplingFuncEdgesList.elements();
	this.nSamplingFunctionBins = (this.samplingFunctionValues).length;
	logger.info("Sampling function is defined");
	logger.info("  nZeros = "+this.nGaps);
	logger.info("  nOnes = "+this.nBins);
    }

    private void setBinHeights(double[] binHeights) {
	this.binHeights = new double[this.nBins];
	this.rates = new double[this.nBins];
	double minBinHeight = Double.MAX_VALUE; 
	double maxBinHeight = -Double.MAX_VALUE; 
	double sumOfBinHeights = 0;
	double sumOfSquaredBinHeights = 0;
	double sumOfRates = 0;
	double sumOfSquaredRates = 0;
	double minRate = Double.MAX_VALUE;
	double maxRate = -Double.MAX_VALUE; 
	int nNaNs = 0;
	int nNonNaNs = 0;
	for ( int i=0; i < this.nBins; i++ ) {
	    this.binHeights[i] = binHeights[i];
	    this.rates[i] = this.binHeights[i]/this.binWidths[i];

	    if ( Double.isNaN(this.binHeights[i]) ) {
		//logger.warn("NaN encountered in binHeights: index "+i+". Excluding from calculations.");
		nNaNs++;
		this.thereAreGaps = true;
	    }
	    else {
		minBinHeight = Math.min(minBinHeight, this.binHeights[i]);
		maxBinHeight = Math.max(maxBinHeight, this.binHeights[i]);
		sumOfBinHeights += this.binHeights[i];
		sumOfSquaredBinHeights += this.binHeights[i]*this.binHeights[i];
		minRate = Math.min(minRate, this.rates[i]);
		maxRate = Math.max(maxRate, this.rates[i]);
		sumOfRates += this.rates[i];
		sumOfSquaredRates += this.rates[i]*this.rates[i];
		nNonNaNs++;
	    }
	}
	this.nNonNaNs = nNonNaNs;
	this.minRate = minRate;
	this.maxRate = maxRate;
	this.minBinHeight = minBinHeight;
	this.maxBinHeight = maxBinHeight;
	this.sumOfBinHeights = sumOfBinHeights;
	this.sumOfSquaredBinHeights = sumOfSquaredBinHeights;
	this.sumOfSquaredRates = sumOfSquaredRates;
	setStatsOnBinHeights();
    }

    private void setStatsOnBinHeights() {
	this.binCentreAtMinBinHeight = this.binCentres[DataUtils.getIndex(this.minBinHeight, this.binHeights)];
	this.binCentreAtMaxBinHeight = this.binCentres[DataUtils.getIndex(this.maxBinHeight, this.binHeights)];
	this.meanBinHeight = this.sumOfBinHeights/this.nNonNaNs;
 	this.varianceInBinHeights = Descriptive.sampleVariance(this.nNonNaNs, this.sumOfBinHeights, this.sumOfSquaredBinHeights);
	this.meanDeviationInBinHeights = Descriptive.meanDeviation(new DoubleArrayList(this.binHeights), this.meanBinHeight);
	this.skewnessInBinHeights = Descriptive.sampleSkew(new DoubleArrayList(this.binHeights), this.meanBinHeight, this.varianceInBinHeights);
	this.skewnessStandardError = Descriptive.sampleSkewStandardError(this.nBins);
	this.kurtosisInBinHeights = Descriptive.sampleKurtosis(new DoubleArrayList(this.binHeights), this.meanBinHeight, this.varianceInBinHeights);
	this.kurtosisStandardError =  Descriptive.sampleKurtosisStandardError(this.nBins);
    }

    private void setErrorsOnRates(double[] errors) {
	this.errorsOnRates = new double[this.nBins];
	for ( int i=0; i < this.nBins; i++ ) {
	    if ( Double.isNaN(errors[i]) ) {
		if ( ! Double.isNaN(this.rates[i]) ) {
		    logger.warn("There is a NaN value in errors whose corresponding rate is not NaN. Setting error from mean binHeights per bin.");
		    double uncertainty = Math.sqrt(this.meanBinHeight);
		    this.errorsOnRates[i] = uncertainty/this.binWidths[i];
		}
	    }
	    else {
		this.errorsOnRates[i] = errors[i];
	    }
	}
	this.errorsAreSet = true;
    }

    //  About Bins
    public int nBins() { return this.nBins; }
    public double xStart() { return this.xStart; }
    public double xStop() { return this.xStop; }
    public double xMid() { return this.xMid; }
    public double xRange() { return this.xRange; }
    public double[] getBinCentres() { return Arrays.copyOf(this.binCentres, this.binCentres.length); }
    public double[] getBinWidths() { return Arrays.copyOf(this.binWidths, this.binWidths.length); }
    public double[] getHalfBinWidths() { return Arrays.copyOf(this.halfBinWidths, this.halfBinWidths.length); }
    public double[] getBinEdges() { return Arrays.copyOf(this.binEdges, this.binEdges.length); }
    public double[] getLeftBinEdges() { return Arrays.copyOf(this.leftBinEdges, this.leftBinEdges.length); }
    public double[] getRightBinEdges() { return Arrays.copyOf(this.rightBinEdges, this.rightBinEdges.length); }
    public double binCentreAtMinBinHeight() { return this.binCentreAtMinBinHeight; }
    public double binCentreAtMaxBinHeight() { return this.binCentreAtMaxBinHeight; }
    public double minBinWidth() { return this.minBinWidth; }
    public double maxBinWidth() { return this.maxBinWidth; }
    public double avgBinWidth() { return this.avgBinWidth; }
    public double binWidth()  throws TimeSeriesException {
	if ( this.binWidthIsConstant )
	    return binWidths[0];
	else {
	    throw new TimeSeriesException("BinWidth is not constant. Use getBinWidths()");
	}
    }

    //  About Gaps
    public int nGaps() { return this.nGaps; }
    public double[] getGapEdges() { return Arrays.copyOf(this.gapEdges, this.gapEdges.length); }
    public double[] getGapLengths() { return Arrays.copyOf(this.gapLengths, this.gapLengths.length); }
    public double minGap() { return this.minGap; }
    public double maxGap() { return this.maxGap; }
    public double meanGap() { return this.meanGap; }
    public double sumOfGaps() { return this.sumOfGaps; }
    public int nSamplingFunctionBins() { return this.nSamplingFunctionBins; }
    public double[] getSamplingFunctionValues() { return Arrays.copyOf(this.samplingFunctionValues, this.samplingFunctionValues.length); }
    public double[] getSamplingFunctionBinEdges() { return Arrays.copyOf(this.samplingFunctionBinEdges, this.samplingFunctionBinEdges.length); }

    //  About Intensities
    public double[] getBinHeights() { return Arrays.copyOf(this.binHeights, this.binHeights.length); }
    public double sumOfBinHeights() { return this.sumOfBinHeights; }
    public double meanBinHeight() { return this.meanBinHeight; }
    public double minBinHeight() { return this.minBinHeight; }
    public double maxBinHeight() { return this.maxBinHeight; }
    public double varianceInBinHeights() { return this.varianceInBinHeights; }
    public double meanDeviationInBinHeights() { return this.meanDeviationInBinHeights; }
    public double kurtosisInBinHeights() { return this.kurtosisInBinHeights; }
    public double kurtosisStandardError() { return this.kurtosisStandardError; }
    public double skewnessInBinHeights() { return this.skewnessInBinHeights; }
    public double skewnessStandardError() { return this.skewnessStandardError; }
    public double fractionalRms() { return this.fractionalRms; }

    //  Boolean checkers
    public boolean binWidthIsConstant() { return this.binWidthIsConstant; }
    public boolean thereAreGaps() { return this.thereAreGaps; }
    public boolean errorsAreSet() { return this.errorsAreSet; }

    //  Write to ascii file
    public void writeBinHeightsAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeBinHeightsAsQDP(this, filename);
    }
    public void writeBinHeightsAsQDP(double[] function, String filename) throws IOException {
	TimeSeriesWriter.writeBinHeightsAsQDP(this, function, filename);
    }
    public void writeRatesAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, filename);
    }
    public void writeRatesAsQDP(double[] function, String filename) throws IOException {
	TimeSeriesWriter.writeRatesAsQDP(this, function, filename);
    }
    public void writeBinHeightsAndSamplingAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeBinHeightsAndSamplingAsQDP(this, filename);
    }
    public void writeRatesAndSamplingAsQDP(String filename) throws IOException {
	TimeSeriesWriter.writeRatesAndSamplingAsQDP(this, filename);
    }

}