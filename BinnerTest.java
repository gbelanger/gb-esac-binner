package gb.esac.binner;

import gb.esac.eventlist.EventList;
import gb.esac.montecarlo.WhiteNoiseGenerator;
import java.io.IOException;
import java.util.Arrays;
import org.apache.log4j.Logger;

public class BinnerTest { 

    private static Logger logger = Logger.getLogger(BinnerTest.class);

    public static void main(String[] args) throws Exception  {
	double mean = 5;
	double duration = 40;
	double[] times = WhiteNoiseGenerator.generateArrivalTimes(mean, duration);
	EventList evlist = new EventList(times);
	evlist.writeTimesAsQDP("times.qdp");
	binDataMethodOneTest(times);
	binDataMethodTwoTest(times);
    }


    private static void binDataMethodOneTest(double[] data) throws IOException {
	
	// public static double[][] binData(double[] data, int nbins) {	
	
	int nBins = 5;
	double[][] binnedTimes = Binner.binData(data, nBins);
	double[] binEdges = binnedTimes[1];
	for ( int i=0; i < nBins; i++ ) {
	    double binCentre = (binEdges[2*i+1] + binEdges[2*i])/2;
	    System.out.println(binCentre+"\t"+binnedTimes[0][i]);
	}
    }


    private static void binDataMethodTwoTest(double[] data) {

	// public static double[][] binData(double[] data, double[] binEdges) {
	
	int nBins = 5;
	double xmin = data[0];
	double[] binEdges = new double[] {0, 8, 8, 16, 16, 24, 24, 32, 32, 40};
	for ( int i=0; i < binEdges.length; i++ ) {
	    binEdges[i] += xmin;
	}
	double[] binnedTimes = Binner.binData(data, binEdges);
	for ( int i=0; i < nBins; i++ ) {
	    double binCentre = (binnedTimes[2*i] + binnedTimes[2*i+1])/2;
	    System.out.println(binCentre+"\t"+binnedTimes[i]);
	}

    }

}
