package gb.esac.binner;

public class SamplingPatternApplicatorTester {

	public static void main(String[] args) throws Exception {

		//  Create simple intensity bins
		double[] leftEdges = new double[] {0, 1, 2, 3};//, 2};
		double[] rightEdges = new double[] {1, 2, 3, 4};//, 3};
		double[] intensities = new double[] {4, 5, 6, 7};//, 1};
		double[] errors = new double[] {1, 2, 3, 4};//, 1};
		int nOld = intensities.length;
		DensityBin[] oldIntensityBins = new DensityBin[nOld];
		for (int i = 0; i < nOld; i++) {
		    oldIntensityBins[i] = new DensityBin(leftEdges[i], rightEdges[i], intensities[i], errors[i]);
		    System.out.println(oldIntensityBins[i].getLeftEdge()+" "+oldIntensityBins[i].getRightEdge()+" "+oldIntensityBins[i].getValue());
		}

		//  Create sampling bins
		int nNew = 5;
		double[] binEdges = BinningUtils.getBinEdges(0, rightEdges[rightEdges.length-1], nNew);
		double[] values = new double[nNew];
		AbsoluteQuantityBin[] samplingBins = new AbsoluteQuantityBin[nNew];
		for (int i = 0; i < nNew; i++) {
			values[i] = 1;
		    samplingBins[i] = new AbsoluteQuantityBin(binEdges[2*i], binEdges[2*i+1], values[i]);
		}
		// for (AbsoluteQuantityBin sampling : samplingBins) {
		// 	System.out.println(sampling.getLeftEdge()+" "+sampling.getRightEdge()+" "+sampling.getValue());
		// }

		// Resample intensity bins into the new bins
		DensityBin[] newIntensityBins = SamplingPatternApplicator.applySamplingPattern(oldIntensityBins,samplingBins);
		System.out.println("Resulting observed bins:");
		for (DensityBin bin : newIntensityBins) {
			System.out.println("["+bin.getLeftEdge()+", "+bin.getRightEdge()+"] : "+bin.getValue()+" +/- "+bin.getError());
		}

	}

}
