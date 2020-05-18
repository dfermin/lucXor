/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor.algorithm;

import gnu.trove.map.TMap;
import gnu.trove.map.hash.THashMap;
import lucxor.common.PSMList;
import lucxor.utils.Constants;
import lucxor.LucXorConfiguration;
import lucxor.utils.MathFunctions;
import lucxor.common.PSM;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.*;

/**
 *
 * @author dfermin
 */
public class FLR {

	private final List<PSM> targetPSMs;
	private final List<PSM> decoyPSMs;
    private TMap<Integer, double[]> minorMapG;
	private TMap<Integer, double[]> minorMapL;

    // f0 = decoy
    // f1 = real

	private double[] pos;
	private double[] neg;
	private double[] tickMarks;
	private double[] f0;
	private double[] f1;
	private double[] globalFDR;
	private double[] localFDR;
	private double maxDeltaScore;
	
	// variance of delta score for both distributions
	private double deltaScoreVar_pos;
	private double deltaScoreVar_neg;
	
	// mean of delta score for both distributions
	private double deltaScoreMu_pos;
	private double deltaScoreMu_neg;
	
	// bandwidth for histograms
	private double bw_real;
	private double bw_decoy;
	private int Nreal;
	private int Ndecoy;
	
	private static final int NMARKS = 10001; // we want N+1 bins for the FLR


	public FLR() {
		targetPSMs = new ArrayList<>();
		decoyPSMs = new ArrayList<>();
		maxDeltaScore = 0d;
		Nreal = 0;
		Ndecoy = 0;
	}


	
	// Using primitive arrays is a lot faster
	public void prepArrays() {
		pos = new double[ targetPSMs.size() ];
		neg = new double[ decoyPSMs.size() ];
		
		for(int i = 0; i < targetPSMs.size(); i++) {
			pos[i] = targetPSMs.get(i).getDeltaScore();
		}
		
		for(int i = 0; i < decoyPSMs.size(); i++) {
			neg[i] = decoyPSMs.get(i).getDeltaScore();
		}
		
		Nreal = pos.length;
		Ndecoy = neg.length;
	}
	
	
	public void initializeTickMarks() {
		
		maxDeltaScore *= 1.001; // need maxDeltaScore to be slightly larger for bins

		tickMarks = new double[NMARKS];
		for(int i = 0; i < NMARKS; i++) {
			double x = ((double) i) * maxDeltaScore / ((double) NMARKS);
			tickMarks[i] = x;
		}
		
		
		// initialize other arrayLists
		localFDR = new double[ Nreal ];
		globalFDR = new double[ Nreal ];

		calcDeltaScoreMean();
		calcDeltaScoreVar();
		
		getBandWidth(Constants.DECOY); // decoys
		getBandWidth(Constants.REAL); // real
		
		System.err.print(
			"FLR bandwidth (pos): " + MathFunctions.roundDouble(bw_real, 6) + "\n" +
			"FLR bandwidth (neg): " + MathFunctions.roundDouble(bw_decoy, 6) + "\n"
		);
		
	}

	
	private void getBandWidth(int dataType) {
		double sigma;
		double N;
		double result;
		double x;
		
		if(dataType == Constants.REAL) { // real
			sigma = Math.sqrt(deltaScoreVar_pos);
			N = (double) targetPSMs.size();
			x = Math.pow(N, 0.2);
			
			result = 1.06 * (sigma / x);
			bw_real = result;
			
		}
		
		if(dataType == Constants.DECOY) { // decoy
			sigma = Math.sqrt(deltaScoreVar_neg);
			N = (double) decoyPSMs.size();
			x = Math.pow(N, 0.2);
			
			result = 1.06 * (sigma / x);
			bw_decoy = result;
		}		
	}

	
	// Compute the mean of the delta score for each data type
	private void calcDeltaScoreMean() {
		double sum = 0d;
		double N;
		
		// for the forwards
		for(double d : pos) sum += d;
		N = (double) pos.length;
		deltaScoreMu_pos = sum / N;
		
		
		// for the decoys
		sum = 0d;
		for(double d : neg) sum += d;
		N = (double) neg.length;
		deltaScoreMu_neg = sum / N;
	}

	
	// Compute the variance of the delta scores
	private void calcDeltaScoreVar() {
		
		double v = 0;
		double N;
		double x;
		
		// for the forwards
		for(double d : pos) {
			x = d - deltaScoreMu_pos;
			v += Math.pow( x, 2.0 ); 
		}
		N = (double) pos.length - 1;
		deltaScoreVar_pos = v / N;
		
		
		// for the decoys
		v = 0d;
		for(double d : neg) {
			x = d - deltaScoreMu_neg;
			v += Math.pow( x, 2.0 ); 
		}
		N = (double) neg.length - 1;
		deltaScoreVar_neg = v / N;
	}


	/**
	 * Function evaluates each deltaScore at each tick mark storing the result in
	 * either f0 (decoys) or f1 (real)
	 * @param dataType
	 * @throws InterruptedException
	 * @throws ExecutionException
	 */
	public void evalTickMarks(int dataType) throws InterruptedException, ExecutionException {
		
		double kernelResult;
		double[] dataAry = null;
		double[] res = null;
		double bw = 0;
		double N = 0;
		
		if(dataType == Constants.DECOY) { // decoys
			dataAry = neg;
			bw = bw_decoy;
			N = (double) Ndecoy;
			res = new double[ NMARKS ];
			f0 = new double[ NMARKS ];
		}
		
		if(dataType == Constants.REAL) { // real
			dataAry = pos;
			bw = bw_real;
			N = (double) Nreal;
			res = new double[ NMARKS ];
			f1 = new double[ NMARKS ];
		}
		
		// how big a chunk of the forward deltaScores to process per cpu
		int block = dataAry.length / LucXorConfiguration.getNumThreads();
		
		for(int i = 0; i < NMARKS; i++) {
			double tic = tickMarks[i];
			
			// Threadpool of executor objects
			ExecutorService pool = Executors.newFixedThreadPool(LucXorConfiguration.getNumThreads());

			// Create a list to hold the tasks to be performed
			List< Future<Double> > taskList = new ArrayList<>(LucXorConfiguration.getNumThreads());

			
			// break up the data in pos into thread chunks
			for(int cpu = 0; cpu < LucXorConfiguration.getNumThreads(); cpu++) {
				int start = cpu * block;
				int end   = start + block;
				if(cpu == (LucXorConfiguration.getNumThreads()-1)) end = dataAry.length;
				
				double[] subAry = Arrays.copyOfRange(dataAry, start, end);
				
				// construct the tasks to submit to the thread pool
				Callable<Double> C = new NormalDensityWorkerThread(subAry, tic, bw);
				Future<Double> task = pool.submit(C);
				taskList.add(task); // put the "task" to be executed into the queue
			}
			
			// Launch each task
			kernelResult = 0;
			for(Future<Double> task : taskList) {
				kernelResult += task.get(); // actually start the job & return result
			}
			kernelResult /= (N * bw);
			
			if(kernelResult > Constants.TINY_NUM) res[i] = kernelResult;
			else res[i] = Constants.TINY_NUM;
			
			pool.shutdown(); // bring down the threadpool
		}
		
		if(dataType == Constants.DECOY) f0 = res;
		
		if(dataType == Constants.REAL) f1 = res;

	}

	
	// Function computes the local and global FDRs (FLR here)
	public void calcBothFDRs() {
		double AUC_rev_0; // Area-Under Curve from end of tick marks working backwards (f0 data)
		double AUC_rev_1; // Area-Under Curve from end of tick marks working backwards (f1 data)
		double ratio;
		
		double Nreal2 = (double) Nreal;
		double Ndecoy2 = (double) Ndecoy;

		int i = 0;
		for(double tmp_score : pos) {
			
			if(tmp_score < 0.1) tmp_score = 0.1;
			
			
			// GlobalFLR
			AUC_rev_0 = getGlobalAUC(tmp_score, Constants.DECOY);
			AUC_rev_1 = getGlobalAUC(tmp_score, Constants.REAL);
			
			ratio = (Ndecoy2/Nreal2) * (AUC_rev_0 / AUC_rev_1);
			globalFDR[i] = ratio;
			
			// localFDR
			AUC_rev_0 = getLocalAUC(tmp_score, Constants.DECOY);
			AUC_rev_1 = getLocalAUC(tmp_score, Constants.REAL);
			
			ratio = (Ndecoy2/Nreal2) * (AUC_rev_0 / AUC_rev_1);
			localFDR[i] = ratio;

			i++; // increment array index
			
		}
	}
	
	
	// Function computes the area under the curve (AUC) for the bin that contains the passed score
	private double getLocalAUC(double x, int whichF) {
		double result = 0d;
		
		double a, b;
		double sum = 0d;
		double tmp1, tmp2;
		double fx;
		
		double start_tick = tickMarks[0];
		double end_tick = tickMarks[ (NMARKS-1) ];
		double start_val, end_val;
		
		
		// For f0
		if(whichF == Constants.DECOY) { // decoys
			start_val = f0[0];
			end_val = f0[ (NMARKS-1) ];
			
			//determine which two bins encompass the value of 'x'
			for(int j = (NMARKS-1); j >= 1; j--) { // working backwards here
				int i = j -1;
				a = tickMarks[i];
				b = tickMarks[j];
				
				if(x >= a) {
					// We have reached the bin that contains 'x'
					// Now compute the area under the curve upto & including 'x'
					tmp1 = (b-x) / (b-a);
					tmp2 = (x-a) / (b-a);
					fx = (tmp1 * f0[i]) + (tmp2 * f0[j]);
					sum += fx;
					break; // leave loop since you are done.
				}
			}
			
			if(x <= start_tick) sum = start_val;
			else if( x >= end_tick) sum = end_val;
			
			result = sum;
		}
		
		
		// For f1
		if(whichF == Constants.REAL) {
			start_val = f1[0];
			end_val = f1[ (NMARKS-1) ];
			
			//determine which two bins encompass the value of 'x'
			for(int j = (NMARKS-1); j >= 1; j--) { // working backwards here
				int i = j -1;
				a = tickMarks[i];
				b = tickMarks[j];
				
				if(x >= a) {
					// We have reached the bin that contains 'x'
					// Now compute the area under the curve upto & including 'x'
					tmp1 = (b-x) / (b-a);
					tmp2 = (x-a) / (b-a);
					fx = (tmp1 * f1[i]) + (tmp2 * f1[j]);
					sum += fx;
					break; // leave loop since you are done.
				}
			}
			if(x <= start_tick) sum = start_val;
			else if( x >= end_tick) sum = end_val;
			
			result = sum;
		}
		
		return result;
	}
	
	
	// Function computes the area under the cuver for the bin that contains the
	// passed score. This function is identical to getLocalAUC except that it
	// computes the area before and after the bin containing x
	private double getGlobalAUC(double x, int whichF) {

		double a, b;
		double sum = 0d;
		double tmp1, tmp2;
		double fx;
		
		// For f0
		if(whichF == Constants.DECOY) { // decoy
			sum = 0d;
			//determine which two bins encompass the value of 'x'
			for(int j = (NMARKS-1); j >= 1; j--) { // working backwards here
				int i = j -1;
				a = tickMarks[i];
				b = tickMarks[j];
				
				if(x < a) sum += (b - a) * (0.5 * (f0[j] + f0[i]));
				else {
					// We have reached the bin that contains 'x'
					// Now compute the area under the curve upto & including 'x'
					tmp1 = (b-x) / (b-a);
					tmp2 = (x-a) / (b-a);
					fx = (tmp1 * f0[i]) + (tmp2 * f0[j]);
					sum += (b-x) * (0.5 * (fx + f0[j]));
					break; // leave loop since you are done.
				}
			}
		}
		
		
		// For f1
		if(whichF == Constants.REAL) { // real
			sum = 0d;
			//determine which two bins encompass the value of 'x'
			for(int j = (NMARKS-1); j >= 1; j--) { // working backwards here
				int i = j -1;
				a = tickMarks[i];
				b = tickMarks[j];
				
				if(x < a) sum += (b - a) * (0.5 * (f1[j] + f1[i]));
				else {
					// We have reached the bin that contains 'x'
					// Now compute the area under the curve upto & including 'x'
					tmp1 = (b-x) / (b-a);
					tmp2 = (x-a) / (b-a);
					fx = (tmp1 * f1[i]) + (tmp2 * f1[j]);
					sum += (b-x) * (0.5 * (fx + f1[j]));
					break; // leave loop since you are done.
				}
			}
		}
		
		return sum;
	}
	
	
	// Function prepares minorMaps for FDR data
	public void setMinorMaps() {
		
		ArrayList<Double> scoreList = new ArrayList<>();
		HashMap<Double, double[]> localMap = new HashMap<>();
        double[] FDRary;
		double ds; // deltaScore
		double FDR;  // FDR
		int i;


        // 0 = globalFDR, 1 = localFDR
        for(int iter = 0; iter < 2; iter++) {

            if(iter == 0) {
                FDRary = globalFDR;
                minorMapG = new THashMap<>();
            }
            else {
                FDRary = localFDR;
                minorMapL = new THashMap<>();
            }

            localMap.clear();
            scoreList.clear();

			for(i = 0; i < Nreal; i++) {
                ds = pos[i]; // delta score
                FDR = FDRary[i];
                scoreList.add(ds);

                double[] X = new double[2];
                X[0] = ds;
                X[1] = FDR;
                localMap.put(ds, X);
            }

            Collections.sort(scoreList); // sort delta score values from low to high

            i = 0;
            for(double L : scoreList) {
                double[] pair = localMap.get(L);

                // Store the pair (deltascore, FDR) into the minorMap
                // The map is indexed from i = 0 to i < N where at i = 0
                // the lowest delta score and it's associated FDR are stored
                if(iter == 0) minorMapG.put(i, pair);
                else minorMapL.put(i, pair);

                i++;
            }
        } // end loop over iter
    }


	/**
	 * Function performs minorization of the FDR data to smooth out the curve
	 * in the event the data is too sparse. In big data sets, this will only improve things
	 */
	public void performMinorization() {
        int i,j,N;
        int curStart, curEnd;
        double f_expect;
        double slope;
        boolean cont;

        double[] deqPtr;
        TMap<Integer, double[]> mapPtr;

        // 0 = global FDR, 1 = localFDR
        for(int iter = 0; iter < 2; iter++) {

            if(iter == 0) {
                deqPtr = globalFDR;
                mapPtr = minorMapG;
            }
            else {
                deqPtr = localFDR;
                mapPtr = minorMapL;
            }

            N = mapPtr.size();

            // initialize arrays for this iteration
            double[] x = new double[N];
            double[] f = new double[N];
            double[] fcopy = new double[N];
            double[] forig = new double[N];
            boolean[] isMinorPoint = new boolean[N];

            for(i = 0; i < N; i++) {
                double[] M = mapPtr.get(i);

                x[i] = M[0]; // delta score
                f[i] = M[1]; // FDR
                forig[i] = f[i];
                fcopy[i] = 0;
                isMinorPoint[i] = false;
            }

            // find the minimum value in f and record it's index too
            int minId = 0;
            double minVal = f[0];
            for(i = 1; i < N; i++) {
                if( f[i] < minVal ) {
                    minVal = f[i];
                    minId = i;
                }
            }

            slope = (0.0 - minVal) / ( maxDeltaScore * 1.1 - x[minId] );
            i = minId;
            while(i < N) {
                f[i] = minVal + slope * ( x[i] - x[minId] );
                i++;
            }



            // find the maximum value in f and record it's index too
            int maxId = 0;
            double maxVal = f[0];
            i = 1;
            while(x[i] < (x[(N-1)]/2.0)) {
                if(f[i] >= maxVal) {
                    maxVal = f[i];
                    maxId = i;
                }
                i++;
            }

            slope = maxVal / ( x[maxId] - x[(N-1)] );
            i = maxId - 1;
            while(i >= 0) {
                f[i] = maxVal - slope * ( x[maxId] - x[i] );
                i--;
            }

            for(i = 0; i < maxId; i++) isMinorPoint[i] = true;
            curStart = maxId;
            curEnd = maxId + 1;

            while( x[curStart] >= x[curEnd] ) curEnd++;


            while( curStart < (N-1) ) {

                slope = ( f[curEnd] - f[curStart] ) / ( x[curEnd] - x[curStart] );

                if(slope > 0.0) curEnd++;
                else {
                    slope = ( f[curEnd] - f[curStart] ) / ( x[curEnd] - x[curStart] );

                    cont = true;
                    i = curStart+1;
                    while( cont && (i < N) && (i < curEnd) ) {
                        f_expect = f[curStart] + slope * (x[i] - x[curStart]);
						cont = !(f_expect > f[i]);
                        i++;
                    }

                    if(cont) {
                        isMinorPoint[curEnd] = true;
                        curStart = curEnd;
                        curEnd = curStart + 1;
                        //while( (x[curStart] >= x[curEnd]) && (curEnd < N) ) curEnd++;
                    }
                    else curEnd++;
                }
            }

            isMinorPoint[ N-1 ] = true; // check if this value is bigger than the second last minor point
            for(i = 0; i < N; i++) fcopy[i] = f[i];

            curStart = 0;
            curEnd = curStart + 1;
            while(!isMinorPoint[curEnd]) curEnd++;

            while( (curStart < (N-1)) && (curEnd < N) ) {
                i = curStart + 1;
                slope = (fcopy[curEnd] - fcopy[curStart]) / (x[curEnd] - x[curStart]);

                while(i < curEnd) {
                    f_expect = fcopy[curStart] + slope * (x[i] - x[curStart]);
                    if(fcopy[i] > f_expect) f[i] = f_expect;
                    i++;
                }

                curStart = curEnd;
                curEnd = curStart + 1;
                if(curEnd >= N) curEnd = N-1;
                while( (!isMinorPoint[curEnd]) && (curEnd < N) ) curEnd++;
            }


            // now map back to globalFDR: f[] has the minorized FDR values in increasing order of x[], which holds delta score
            for(i = 0; i < N; i++) {
                for(j = 0; j < N; j++) {
                    if(pos[i] == x[j]) {
                        deqPtr[i] = f[j];
                        break;
                    }
                }
            }

            if(iter == 0) globalFDR = deqPtr;
            else localFDR = deqPtr;
        } // end loop over iter

    }

	
	
	// Function assigns FDR values to the PSMs
	public void assignFDRs(PSMList psmList) {
		
		int N = targetPSMs.size();
		
		for(int i = 0; i < N; i++) {
            double g_FDR = (globalFDR[i] > 1.0 ? 1.0 : globalFDR[i]);
            double l_FDR = (localFDR[i] > 1.0 ? 1.0 : localFDR[i]);

			targetPSMs.get(i).setGlobalFDR(g_FDR);
			targetPSMs.get(i).setLocalFDR(l_FDR);
		}

		for (PSM realPSM : targetPSMs) {
			for (PSM p : psmList) {
				if (p.getSpecId().equalsIgnoreCase(realPSM.getSpecId())) {
					p.setGlobalFDR(realPSM.getGlobalFDR());
					p.setLocalFDR(realPSM.getLocalFDR());
					break;
				}
			}
		}
		
		// fill in remaining PSMs
		for(PSM p : psmList) {
			if(p.isDecoy()) {
				p.setGlobalFDR(Double.NaN);
				p.setLocalFDR(Double.NaN);
			}
			else if(p.getDeltaScore() <= Constants.MIN_DELTA_SCORE) {
				p.setGlobalFDR(1.0);
				p.setLocalFDR(1.0);
			}
		}
	}

	public void addDecoy(PSM psm){
		this.decoyPSMs.add(psm);
	}

	public void addTarget(PSM psm){
		this.targetPSMs.add(psm);
	}

	public void setMaxDeltaScore(double maxDeltaScore) {
		this.maxDeltaScore = maxDeltaScore;
	}
}
