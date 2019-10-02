package gb.esac.binner;

import gb.esac.io.AsciiDataFileWriter;
import gb.esac.timeseries.TimeSeries;
import gb.esac.timeseries.TimeSeriesFileException;
import gb.esac.timeseries.TimeSeriesMaker;
import gb.esac.timeseries.TimeSeriesResampler;
import gb.esac.tools.Stats;
import hep.aida.IHistogram1D;
import java.io.IOException;
import org.apache.log4j.Logger;

/**

   The class <code>ResamplerTest</code> is for testing <code>Resampler</code>.

   @author <a href="mailto: guilaume.belanger@esa.int">Guillaume Belanger</a>, ESA/ESAC.
   @created August 2010
   @version March 2013

**/

public class ResamplerTest {

    private static Logger logger = Logger.getLogger(ResamplerTest.class);

    public static void main(String[] args) throws Exception  {

	resampleRates();
    }

    public static void resampleRates() throws Exception {

	// public static double[][] resampleRates(double[] rates, double[] errors, double[] oldBinEdges, double[] newBinEdges) {	


	int nLoops = 1000;
	double[] integralsOfRates = new double[nLoops];
	double[] integralsOfResampledRates = new double[nLoops];
	double[] integralsOfRebinnedRates = new double[nLoops];
	for ( int k=0; k < nLoops; k++ ) {

	    int nOldBins = 50;
	    double tStart = 0;
	    double tStop = 100;
	    double oldBinWidth = (tStop-tStart)/nOldBins;
	    double[] oldBinEdges = BinningUtils.getBinEdges(tStart, tStop, nOldBins);
	    double[] rates = new double[nOldBins];
	    double[] errors = new double[nOldBins];
	    for ( int i=0; i < nOldBins; i++ ) {
		rates[i] = 10*Math.random();
		errors[i] = 1;
	    }


	    // System.out.println("DEV /XS");
	    // System.out.println("READ SERR 2");
	    // System.out.println("LINE STEP ON");
	    // System.out.println("MA 21 ON");
	    // System.out.println("CS 1.3");
	    // System.out.println("SKIP SINGLE");
	    double integral = 0;
	    for ( int i=0; i < nOldBins; i++ ) {
		double binCentre = (oldBinEdges[2*i] + oldBinEdges[2*i+1])/2;
		integral += rates[i]*oldBinWidth;
		// System.out.println(binCentre+"\t"+rates[i]+"\t"+errors[i]);
	    }
	    integralsOfRates[k] = integral;
	    // System.out.println("NO NO NO");


	    int nNewBins = 24;
	    //double[] newBinEdges = BinningUtils.getLinearBinEdges(tStart, tStop, nNewBins);
	    //double[][] ratesAndErrors = Resampler.resampleRates(rates, errors, oldBinEdges, newBinEdges);
	    double newBinWidth = (tStop-tStart)/nNewBins;
	    double[][] ratesAndErrors = Resampler.resample(rates, errors, oldBinEdges, newBinWidth);
	    double[] r = ratesAndErrors[0];
	    double[] e = ratesAndErrors[1];
	    double[] newBinEdges = ratesAndErrors[2];


	    integral=0;
	    for ( int i=0; i < nNewBins; i++ ) {
		double binCentre = (newBinEdges[2*i] + newBinEdges[2*i+1])/2;
		integral += r[i]*newBinWidth;
		// System.out.println(binCentre+"\t"+r[i]+"\t"+e[i]);
	    }
	    integralsOfResampledRates[k] = integral;
	    // System.out.println("NO NO NO");


// 	    double[][] rebRatesAndErrors = Rebinner.rebinRates(rates, errors, oldBinEdges, newBinWidth);
// 	    r = rebRatesAndErrors[0];
// 	    e = rebRatesAndErrors[1];
// 	    newBinEdges = rebRatesAndErrors[2];
// 	    integral=0;
// 	    for ( int i=0; i < nNewBins; i++ ) {
// 		double binCentre = (newBinEdges[2*i] + newBinEdges[2*i+1])/2;
// 		integral += r[i]*newBinWidth;
// 		// System.out.println(binCentre+"\t"+r[i]+"\t"+e[i]);
// 	    }
// 	    integralsOfRebinnedRates[k] = integral;


	// 	nNewBins = 30;
	// 	//newBinEdges = BinningUtils.getLinearBinEdges(tStart, tStop, nNewBins);
	// 	//ratesAndErrors = Resampler.resampleRates(rates, errors, oldBinEdges, newBinEdges);
	// 	newBinWidth = (tStop-tStart)/nNewBins;
	// 	ratesAndErrors = Resampler.resampleRates(rates, errors, oldBinEdges, newBinWidth);
	// 	r = ratesAndErrors[0];
	// 	e = ratesAndErrors[1];
	// 	newBinEdges = ratesAndErrors[2];
	// 	for ( int i=0; i < nNewBins; i++ ) {
	// 	    double binCentre = (newBinEdges[2*i] + newBinEdges[2*i+1])/2;
	// 	    System.out.println(binCentre+"\t"+r[i]+"\t"+e[i]);
	// 	}


	}

	
// 	int nHistoBins = 25;
// 	double xmin = Stats.getMin(new double[] {Stats.getMin(integralsOfRates), Stats.getMin(integralsOfResampledRates), Stats.getMin(integralsOfRebinnedRates)}) - 1e-10;
// 	double xmax = Stats.getMax(new double[] {Stats.getMax(integralsOfRates), Stats.getMax(integralsOfResampledRates), Stats.getMax(integralsOfRebinnedRates)}) + 1e-10;
// 	IHistogram1D histoOfRates = Binner.makeHisto(integralsOfRates, xmin, xmax, nHistoBins);
// 	IHistogram1D histoOfResampledRates = Binner.makeHisto(integralsOfResampledRates, xmin, xmax, nHistoBins);
// 	IHistogram1D histoOfRebinnedRates = Binner.makeHisto(integralsOfRebinnedRates, xmin, xmax, nHistoBins);

// 	AsciiDataFileWriter out = new AsciiDataFileWriter("histoRates.qdp");
// 	out.writeHisto(histoOfRates, "Integral of Rates");
// 	out = new AsciiDataFileWriter("histoResampledRates.qdp");
// 	out.writeHisto(histoOfResampledRates, "Integral of Resampled Rates");
// 	out = new AsciiDataFileWriter("histoRebinnedRates.qdp");
// 	out.writeHisto(histoOfRebinnedRates, "Integral of Rebinned Rates");



	TimeSeries ts = (TimeSeries) TimeSeriesMaker.makeTimeSeries("/Users/gbelanger/javaProgs/flare.fits");
	double[] counts = ts.getBinHeights();
	double[] oldBinEdges = ts.getBinEdges();
	double binWidth = 120.0;
	double[] newBinEdges = BinningUtils.getBinEdges(oldBinEdges[0], oldBinEdges[oldBinEdges.length-1], binWidth);

	TimeSeries ts2 = TimeSeriesResampler.resample(ts, binWidth);
	ts2.writeCountsAsQDP("ts2.qdp");

	double[] resampledCounts = Resampler.resample(counts, oldBinEdges, newBinEdges);
	TimeSeries ts3 = TimeSeriesMaker.makeTimeSeries(newBinEdges, resampledCounts);
	ts2.writeCountsAsQDP("ts3.qdp");
	

    }

}
