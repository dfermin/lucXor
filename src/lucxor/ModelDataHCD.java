/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import lombok.extern.slf4j.Slf4j;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;

/**
 *
 * @author dfermin
 */
@Slf4j
public class ModelDataHCD {

	// --Commented out by Inspection (2020-05-18 00:58):private final static double NORMAL_CONSTANT = 1.0 / Math.sqrt(2.0 * Math.PI);

	private final int chargeState;
	int numPSM;
	private static final int ntick = 2000;

    private double bIntMean;
    private double bIntVar;
    private double yIntMean;
    private double yIntVar;
	private double negIntMean;
    private double negIntVar;
	private double posDistMean;
    private double posDistVar;


	private double[] bTickMarksInt;
    private double[] yTickMarksInt;
    private double[] negTickMarksInt;
	private double[] posTickMarksDist;
	private double[] f_int_b;
    private double[] f_int_y;
    private double[] f_int_neg;
	private double[] f_dist;

    private final double[] b_int;
    private final double[] y_int;
    private final double[] n_int;
	private final double[] pos_dist;

	// Default constructor
	public ModelDataHCD(int z, List<Peak> peaks) {

		chargeState = z;
		numPSM = 0;
        ArrayList<Peak> posPks = new ArrayList<>(peaks.size());
        ArrayList<Peak> negPks = new ArrayList<>(peaks.size());
		
		int b = 0;
		int y = 0;
		int n = 0;
		int p = 0;
		
		for(Peak pk : peaks) {
			if(pk.isMatched()) {
				p++;
				if(pk.getMatchedIonStr().startsWith("b")) { posPks.add(pk); b++; }
				if(pk.getMatchedIonStr().startsWith("y")) { posPks.add(pk); y++; }
			}
			else {
                negPks.add(pk);
                n++;
            }
		}
		
		// To speed things up, we separate the data into primitive arrays
		b_int = new double[ b ];
		y_int = new double[ y ];
		pos_dist = new double[ p ];
		
		
		b = 0;
		y = 0;
		for(Peak pk : posPks) {
			if(pk.getMatchedIonStr().startsWith("b")) {
				b_int[ b ] = pk.getNormIntensity();
				b++;
			}
			else if(pk.getMatchedIonStr().startsWith("y")) {
				y_int[ y ] = pk.getNormIntensity();
				y++;
			}
		}
		
		p = 0;
		for(Peak pk : posPks) {
			pos_dist[ p ] = pk.getDist();
			p++;
		}

        // We will limit the size of the negative distribution to speed things up
        int limitN = (b + y);
        if(limitN < Constants.MIN_NUM_NEG_PKS) limitN += Constants.MIN_NUM_NEG_PKS;

        if(limitN > n) limitN = n; // prevents segfault on insufficient data for modeling

        n_int = new double[ limitN ];
        Collections.shuffle(negPks);
		for(n = 0; n < limitN; n++) {
            Peak pk = negPks.get(n);
			n_int[ n ] = pk.getNormIntensity();
		}
		
		posPks.clear();
		negPks.clear();
	}

	
// --Commented out by Inspection START (2020-05-18 00:55):
//	// This function trims the extreme values out of the distributions
//	public void percentileTrim(char ionType, int dataType, double percentTrim) {
//		double[] d = null;
//		double[] f;
//		int N;
//		int n;
//		int a;
//		int b;
//
//		if(ionType == 'b') {
//			if(dataType == 0) d = b_int;
////			if(dataType == 1) d = b_dist;
//		}
//
//		if(ionType == 'y') {
//			if(dataType == 0) d = y_int;
////			if(dataType == 1) d = y_dist;
//		}
//
//		if(ionType == 'n') d = n_int;
//
//		N = Objects.requireNonNull(d).length;
//
//		n = (int) (((double) N) * (percentTrim * 0.5)); // we want 1/2  of the trim from both sides
//
//		// get the limits
//		a = n;
//		b = (N - 1 - n);
//
//		Arrays.sort(d);
//
//		f = Arrays.copyOfRange(d, a, b);
//
//		if(ionType == 'b') {
//			if(dataType == 0) b_int = f;
////			if(dataType == 1) b_dist = f;
//		}
//
//		if(ionType == 'y') {
//			if(dataType == 0) y_int = f;
////			if(dataType == 1) y_dist = f;
//		}
//
//		if(ionType == 'n') n_int = f;
//
//	}
// --Commented out by Inspection STOP (2020-05-18 00:55)

	
	
	
	// Compute the mean for pos and neg data. This function must be called
	// *before* the calcVar function
	public void calcMean() {
		double N;
		double sum;
		
		
		// B-ions
		N = (double) b_int.length;
		sum = 0;
		for(double b : b_int) sum += b;
		bIntMean = sum / N;
		
		
		// Y-ions
		N = (double) y_int.length;
		sum = 0;
		for(double y : y_int) sum += y;
		yIntMean = sum / N;

		// Positive peak distances
		posDistMean = getMode(pos_dist);
		
		
		// Noise peaks
		N = (double) n_int.length;
		sum = 0;
		for(double n : n_int) sum += n;
		negIntMean = sum / N;
		
		
		// For the negative distribution we assume a uniform distribution over a
		// 1 Dalton window. As such, the log of a uniform distribution between 
		// -1 to 1 is zero
	}
	
	
	// Compute the variance of the data
	public void calcVar() {
		double N, v;
		
		// B-ions
		v = 0;
		N = ((double) b_int.length) - 1.0;
		for(double b : b_int) {
			double x = b - bIntMean;
			v += Math.pow(x, 2.0);
		}
		bIntVar = (v / N);

		
		// Y-ions
		v = 0;
		N = ((double) y_int.length) - 1.0;
		for(double b : y_int) {
			double x = b - yIntMean;
			v += Math.pow(x, 2.0);
		}
		yIntVar = (v / N);
			
		
		
		// Positive Peak Distance
		v = 0;
		N = ((double) pos_dist.length) - 1.0;
		for(double n : pos_dist) {
			double x = n - posDistMean;
			v += Math.pow(x, 2.0);
		}
		posDistVar = (v / N);
		
		
		
		// Noise Peaks
		v = 0;
		N = ((double) n_int.length) - 1.0;
		for(double n : n_int) {
			double x = n - negIntMean;
			v += Math.pow(x, 2.0);
		}
		negIntVar = (v / N);
		
	}


	/**
	 * Function computes the estimated parameters for the non-parametric model
	 * for the intensity of the peaks assigned to this ModelData object
	 * @param ionType ionType
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	void estimateNPIntensity(char ionType) throws InterruptedException, ExecutionException {
		int N;
		double[] norm_ints;
		double[] tickMarksInt;
		double[] f_int;
		double minI, maxI;
		double variance, bw;
		double kernelResult;
		double t;
		
		if(ionType == 'b') {
			N = b_int.length;
			norm_ints = b_int;
			variance = bIntVar;
			log.info("+" + chargeState + "  Estimating NP Model for b-ion intensities");
		}
		else if(ionType == 'y') {
			N = y_int.length;
			norm_ints = y_int;
			variance = yIntVar;
			log.info("+" + chargeState + "  Estimating NP Model for y-ion intensities");
		}
		else {
			N = n_int.length;
			norm_ints = n_int;
			variance = negIntVar;
			log.info("+" + chargeState + "  Estimating NP Model for noise peak intensities");
		}

		// get the range of the normalized intensities
		Arrays.sort(norm_ints);
		minI = norm_ints[0];
		maxI = norm_ints[ (N-1) ];

		double padding = 0.1;
		if(minI < 0) t = minI + (MathFunctions.roundDouble(minI,4) * padding);
		else t = minI - (MathFunctions.roundDouble(minI,4) * padding);
		minI = MathFunctions.roundDouble(t,4);
		
		if(maxI < 0) t = maxI - (MathFunctions.roundDouble(maxI,4) * padding);
		else t = maxI + (MathFunctions.roundDouble(maxI,4) * padding);
		maxI = MathFunctions.roundDouble(t,4);
		
		
		tickMarksInt = new double[ ntick ];

		// construct the tick marks for the x-axis of the plot
		for(int ii = 0; ii < ntick; ii++) {
			double tic = minI + ((double) ii) * (maxI - minI) / ((double) (ntick - 1));
			tickMarksInt[ii] = tic;
		}

		// Compute estimate of positive intensity bandwidth
		double sigma = Math.sqrt( variance );
		bw = ( 1.06 * (sigma / Math.pow(((double)N), 0.2)) ) * 0.5;

		f_int = new double[ ntick ]; // this will hold the estimated height
		
		//how big a chunk of the norm_ints each task will process
		int block = norm_ints.length / Globals.numThreads;
		
		
		// Iterate over each tick mark
		for(int i = 0; i < ntick; i++) {
			double tic = tickMarksInt[i];
			
			// Threadpool of executor objects
			ExecutorService pool = Executors.newFixedThreadPool(Globals.numThreads);
		
			// Create a List to hold the tasks to be performed
			List< Future<Double> > taskList = new ArrayList<>(Globals.numThreads);
		
			
			// break up the data in norm_ints into thread chunks
			for(int cpu = 0; cpu < Globals.numThreads; cpu++) {
				int start = cpu * block;
				int end   = start + block;
				if(cpu == (Globals.numThreads-1) ) end = norm_ints.length;
				
				// get a subset of the data to submit to the worker thread
				double[] subAry = Arrays.copyOfRange(norm_ints, start, end);
				
				// create a task and submit it to the pool
				Callable<Double> C = new NormalDensityWorkerThread(subAry, tic, bw);
				Future<Double> task = pool.submit(C);
				taskList.add(task);

			}
			
			// Actually launch the tasks and store the results to kernelResults
			kernelResult = 0;
			for(Future<Double> future : taskList) {
				kernelResult += future.get();
			}
			kernelResult /= (((double)N) * bw);
			
			if(kernelResult <= Constants.TINY_NUM) kernelResult = Constants.TINY_NUM;
			
			f_int[ i ] = kernelResult;
			
			pool.shutdown(); // shutdown the threadpool
		}
				
		if(ionType == 'b') {
			bTickMarksInt = tickMarksInt;
			f_int_b = f_int;
        }
		else if(ionType == 'y') {
			yTickMarksInt = tickMarksInt;
			f_int_y = f_int;
        }
		else {
			negTickMarksInt = tickMarksInt;
			f_int_neg = f_int;
		}
		
	}

	
	// Function computes the estimated parameters for the non-parametric model
	// for the m/z distances for the positive peaks in this model object
	void estimateNPPosDist() throws InterruptedException, ExecutionException {
		
		int N;
		double min_dist, max_dist;
		double bw;
		double kernelResult;
		double t;
		double[] f_ary;

		log.info("+" + chargeState + "  Estimating NP Model for b and y ions m/z distances");

		N = pos_dist.length;
		
		// get the range of the m/z distance values
		Arrays.sort(pos_dist);
		min_dist = pos_dist[0];
		max_dist = pos_dist[ (N-1) ];
		
		// we want to pad the extremes of the distribution to smooth out the curve
		double padding = 0.1;
		if(min_dist < 0) t = min_dist + (MathFunctions.roundDouble(min_dist, 4) * padding);
		else t = min_dist - (MathFunctions.roundDouble(min_dist, 4) * padding);
		min_dist = MathFunctions.roundDouble(t, 4);
		
		if(max_dist < 0) t = max_dist - (MathFunctions.roundDouble(max_dist, 4) * padding);
		else t = max_dist + (MathFunctions.roundDouble(max_dist, 4) * padding);
		max_dist = MathFunctions.roundDouble(t, 4);
			
		posTickMarksDist = new double[ ntick ];

		for(int ii = 0; ii < ntick; ii++) {
			double tic = min_dist + ((double) ii) * (max_dist - min_dist) / ((double) (ntick - 1));
			posTickMarksDist[ii] = tic;
		}

		// Compute estimate of positive intensity bandwidth
		double sigma = Math.sqrt( posDistVar );
		bw = ( 1.06 * (sigma / Math.pow(((double)N), 0.2)) ) * 0.1;
		
		f_ary = new double[ ntick ];

		//how big a chunk of the norm_ints each task will process
		int block = pos_dist.length / Globals.numThreads;
		
		
		// Iterate over each tick mark
		for(int i = 0; i < ntick; i++) {
			double tic = posTickMarksDist[i];
			
			// Threadpool of executor objects
			ExecutorService pool = Executors.newFixedThreadPool(Globals.numThreads);

			// Create a list to hold the tasks to be performed
			List< Future<Double> > taskList = new ArrayList<>(Globals.numThreads);

			
			// break up the data in D into thread chunks
			for(int cpu = 0; cpu < Globals.numThreads; cpu++) {
				int start = cpu * block;
				int end   = start + block;
				if(cpu == (Globals.numThreads-1)) end = pos_dist.length;
				
				double[] subAry = Arrays.copyOfRange(pos_dist, start, end);
				
				// construct the tasks to submit to the threadpool
				Callable<Double> C = new NormalDensityWorkerThread(subAry, tic, bw);
				Future<Double> task = pool.submit(C);
				taskList.add(task); // put the "task" into the queue to be executed
			}
			
			// Now launch each task
			kernelResult = 0;
			for(Future<Double> task : taskList) {
				kernelResult += task.get(); // run the job and capture the result
			}
			kernelResult /= (((double)N) * bw);
			
			if(kernelResult <= Constants.TINY_NUM) kernelResult = Constants.TINY_NUM;
			
			f_ary[ i ] = kernelResult;
			
			pool.shutdown(); // brint down the threadpool
		}
		
		f_dist = f_ary;
    }
	
	
// --Commented out by Inspection START (2020-05-18 00:55):
//	public static double normalDensity(double curTickMark, double curScore, double h) {
//		double res;
//
//		double x = (curTickMark - curScore) / h;
//
//		res = NORMAL_CONSTANT * Math.exp( (-0.5*x*x) );
//		return res;
//	}
// --Commented out by Inspection STOP (2020-05-18 00:55)


    /**
     * Function computes the mode of the given list of values
     * @param ary
     * @return
     */
	private static double getMode(double[] ary) {
		
		double mode;
		int Nbins = ntick;
		double binWidth = 0.0001;
		final double LIMIT = 0.1;
		
		double[] v = new double[ntick];
		
		
		// iterate over values
		for(double L : ary) {
			for(int j = 0; j < (Nbins - 1); j++) {
				double a = -LIMIT + ( ((double) j) * binWidth);
				double b = a + binWidth;
				
				if( (L >= a) && (L < b) ) {
					v[j] += v[j] + 1.0;
					break;
				}
			}
		}
		
		// Now find the largest value in 'v'
		double maxValue = 0;
		int maxValIdx = 0;
		for(int i = 0; i < Nbins; i++) {
			if(v[i] > maxValue) {
				maxValue = v[i];
				maxValIdx = i;
			}
		}
		
		mode = -LIMIT + ( ((double) maxValIdx) * binWidth ) + (binWidth/2); 
		
		return mode;
	}

	
	// Function returns the log-transformed non-parametric density for the given
	// peak intensity using the specified ion type
	double getLogNPdensityInt(char ionType, double x) {
		double a, b, tmp1, tmp2, fx;
		int N = 0;
		double sum = 0;
		double[] tickMarksInt = null;
		double[] f_int = null;
		
		if(ionType == 'b') {
			tickMarksInt = bTickMarksInt;
			f_int = f_int_b;
		}
		else if(ionType == 'y') {
			tickMarksInt = yTickMarksInt;
			f_int = f_int_y;
		}
		else if(ionType == 'n') {
			tickMarksInt = negTickMarksInt;
			f_int = f_int_neg;
		}
		
		
		double start_tick = Objects.requireNonNull(tickMarksInt)[ 0 ];
		double end_tick = tickMarksInt[ (ntick-1) ];
		double start_val = Objects.requireNonNull(f_int)[0];
		double end_val = f_int[ (ntick-1) ];
		
		// Figure out which two bins encompass the value of 'x'
		for(int j = (ntick-1); j >= 1; j--) { // working backwards
			int i = j - 1;
			a = tickMarksInt[ i ];
			b = tickMarksInt[ j ];
			
			if( (x >= a) && (x < b) ) {
				// We have found the bin that contains 'x'
				// Now compute the area under the curve upto the point that includes 'x'
				tmp1 = (b - x) / (b - a);
				tmp2 = (x - a) / (b - a);
				
				fx = (tmp1 * f_int[ i ]) + (tmp2 * f_int[ j ]);
				sum = fx;
				break; // leave loop, you're done
			}
		}
		
		if(x <= start_tick) sum = start_val;
		else if(x >= end_tick) sum = end_val;
		
		return Math.log(sum);
	}
	
	
	// Function returns the log-transformed non-parametric density for the given
	// peak distance using the positive model
	double getLogNPdensityDistPos(double x) {
		double a, b, tmp1, tmp2, fx;
		double sum = 0;
		
		double start_tick = posTickMarksDist[ 0 ];
		double end_tick = posTickMarksDist[ (ntick-1) ];
		double start_val = f_dist[0];
		double end_val = f_dist[ (ntick-1) ];
		
		// Figure out which two bins encompass the value of 'x'
		for(int j = (ntick-1); j >= 1; j--) { // working backwards
			int i = j - 1;
			a = posTickMarksDist[ i ];
			b = posTickMarksDist[ j ];
			
			if( (x >= a) && (x <b) ) {
				// We have found the bin that contains 'x'
				// Now compute the area under the curve upto the point that includes 'x'
				tmp1 = (b - x) / (b - a);
				tmp2 = (x - a) / (b - a);
				
				fx = (tmp1 * f_dist[ i ]) + (tmp2 * f_dist[ j ]);
				sum = fx;
				break; // leave loop, you're done
			}
		}
		
		if(x <= start_tick) sum = start_val;
		else if(x >= end_tick) sum = end_val;
		
		return Math.log(sum);
	}
	
	
	
	void writeModelPks() throws IOException {	
		File debugF = new File("debug_model_pks_HCD.txt");
		FileWriter fw;
		BufferedWriter bw;
		String line;
		double  dist,  normI;
		
		if(!debugF.exists()) {
			fw = new FileWriter(debugF);
			bw = new BufferedWriter(fw);
			String hdr = "charge\tdataType\tvalue\n";
			bw.write(hdr);
		}
		else{
			fw = new FileWriter(debugF, true); // open for appending
			bw = new BufferedWriter(fw);
		}


		for (double v : b_int) {
			normI = MathFunctions.roundDouble(v, 4);
			line = chargeState + "\tyi\t" +
					normI + "\n";
			bw.write(line);
		}

		for (double v : y_int) {
			normI = MathFunctions.roundDouble(v, 4);
			line = chargeState + "\tyi\t" + normI + "\n";
			bw.write(line);
		}

		for (double v : n_int) {
			normI = MathFunctions.roundDouble(v, 4);
			line = chargeState + "\tni\t" + normI + "\n";
			bw.write(line);
		}

		for (double v : pos_dist) {
			dist = MathFunctions.roundDouble(v, 4);
			line = chargeState + "\td\t" + dist + "\n";
			bw.write(line);
		}
		
		bw.close();
	}
	
	
	public void printStats() {
		
		String line = "Z = " + chargeState + ":  " +
				      "b-ions Intensity: (mean, std): (" +
				MathFunctions.roundDouble(bIntMean, 4) + ", " +
				MathFunctions.roundDouble(Math.sqrt(bIntVar), 4) + ") N = " +
					   b_int.length + "\n" +
					  "Z = " + chargeState + ":  " +
				      "y-ions Intensity: (mean, std): (" +
				MathFunctions.roundDouble(yIntMean, 4) + ", " +
				MathFunctions.roundDouble(Math.sqrt(yIntVar), 4) + ") N = " +
					   y_int.length + "\n" +
				
					  "Z = " + chargeState + ":  " +
				      "Matched Peak Distance: (mean, std): (" +
				MathFunctions.roundDouble(posDistMean, 4) + ", " +
				MathFunctions.roundDouble(Math.sqrt(posDistVar), 4) + ") N = " +
					   pos_dist.length + "\n" +
					   
					   "Z = " + chargeState + ":  " +
					   "Noise peak Intensity: (mean, std): (" +
				MathFunctions.roundDouble(negIntMean, 4) + ", " +
				MathFunctions.roundDouble(Math.sqrt(negIntVar), 4) + ") N = " +
					   n_int.length + "\n";
		
		log.info(line);
	}
	
	
	// Function writes to disk the data for the non-parametric model specified.
	// This is a debugging function
	public void write_density_data(int dataType) throws IOException {
		
		String outName = null;
		FileWriter fw;
		BufferedWriter bw;
		
		switch(dataType) {
			case 1:
				outName = "debug_HCD_nonParam_intensity.txt";
				break;
			case 2:
				outName = "debug_HCD_nonParam_distance.txt";
				break;
			default:
				break;
		}
		
		File outF = new File(Objects.requireNonNull(outName));
		if(!outF.exists()) {
			fw = new FileWriter(outF);
			bw = new BufferedWriter(fw);
			
			if(dataType == 1) bw.write("dataType\tcharge\ttic\tintensity\n");
			if(dataType == 2) bw.write("dataType\tcharge\ttic\tdistance\n");
		}
		else {
			fw = new FileWriter(outF, true); // open for appending
			bw = new BufferedWriter(fw);
		}
		
		
		if(dataType == 1) {
			for(int i = 0; i < ntick; i++) {
				bw.write(
						"b\t" +
								chargeState + "\t" +
								bTickMarksInt[i] + "\t" +
								f_int_b[i] + "\n"
				);
				
				bw.write(
						"y\t" +
								chargeState + "\t" +
								yTickMarksInt[i] + "\t" +
								f_int_y[i] + "\n"
				);
				
				bw.write(
						"n\t" +
								chargeState + "\t" +
								negTickMarksInt[i] + "\t" +
								f_int_neg[i] + "\n"
				);
			}
		}
		if(dataType == 2) {  // distance
			for(int i = 0; i < ntick; i++) {
				bw.write(
						"p\t" +
								chargeState + "\t" +
								posTickMarksDist[i] + "\t" +
								f_dist[i] + "\n"
				);
			}
		}
		
		bw.close();
	}
	
	
}
