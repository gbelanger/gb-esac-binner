package gb.esac.binner;

public class IntensityBinResamplerTester {

    public static void main(String[] args) throws BinningException {

	//  Create 3 simple intensity bins
	double[] leftEdges = new double[] {0, 1};//, 2};
	double[] rightEdges = new double[] {1, 2};//, 3};
	double[] intensities = new double[] {1, 1};//, 1};
	double[] errors = new double[] {1, 1};//, 1};
	int nOld = intensities.length;
	DensityBin[] oldIntensityBins = new DensityBin[nOld];
	for (int i=0; i < nOld; i++) {
	    oldIntensityBins[i] = new DensityBin(leftEdges[i], rightEdges[i], intensities[i], errors[i]);
	}

	//  Create new bins according to which will be resampled the intensity bins
	int nNew = 1;
	double[] binEdges = BinningUtils.getBinEdges(0, rightEdges[rightEdges.length-1], nNew);
	Bin[] newBins = new Bin[nNew];
	for (int i = 0; i < nNew; i++) {
	    newBins[i] = new Bin(binEdges[2*i], binEdges[2*i+1]);
	}

	// Resample intensity bins into the new bins
	DensityBin[] newIntensityBins = IntensityBinResampler.resample(oldIntensityBins, newBins);
    }
    
}
