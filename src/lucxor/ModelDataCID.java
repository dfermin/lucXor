/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucxor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author dfermin
 */
public class ModelDataCID {


	int chargeState;
	int numPSM; // holds number of PSMs used for this charge state
	double mu_int_B, mu_int_Y, mu_int_U;
	double var_int_B, var_int_Y, var_int_U;
	double mu_dist_B, mu_dist_Y, mu_dist_U;
	double var_dist_B, var_dist_Y, var_dist_U;

	double[] b_intensity = null;
	double[] b_distance = null;
	double[] y_intensity = null;
	double[] y_distance = null;
	double[] u_intensity = null;
	double[] u_distance = null;

	
	// This is an adjustment factor for the standard deviation suggested by hwchoi
	final double CID_ADJUST = 16.0 / 25.0; 
	
	public ModelDataCID(int z, List<Peak> peaks) {
        chargeState = z;
        numPSM = 0;

        int Nb = 0, Ny = 0, Nu = 0;

        for (Peak pk : peaks) {
            if (pk.isMatched()) {
                if (pk.getMatchedIonStr().startsWith("b")) {
                    Nb++;
                }
                if (pk.getMatchedIonStr().startsWith("y")) {
                    Ny++;
                }
            } else {
                Nu++;
            } //unmatched_pks.add(pk);
        }


        b_intensity = new double[Nb];
        b_distance = new double[Nb];
        y_intensity = new double[Ny];
        y_distance = new double[Ny];

        Nb = 0;
        Ny = 0;

        for (Peak pk : peaks) {
            if (pk.isMatched()) {
                if (pk.getMatchedIonStr().startsWith("b")) {
                    b_intensity[Nb] = pk.getNormIntensity();
                    b_distance[Nb] = pk.getDist();
                    Nb++;
                }

                if (pk.getMatchedIonStr().startsWith("y")) {
                    y_intensity[Ny] = pk.getNormIntensity();
                    y_distance[Ny] = pk.getDist();
                    Ny++;
                }
            }
        }


        // We will limit the size of the negative distribution to speed things up
        int limitN = (Nb + Ny);
        if(limitN < Constants.MIN_NUM_NEG_PKS) limitN += Constants.MIN_NUM_NEG_PKS;

        if(limitN > Nu) limitN = Nu; // prevents segfault on insufficient data for modeling
        u_intensity = new double[ limitN ];
        u_distance = new double[ limitN ];

        ArrayList<Peak> negPks = new ArrayList<>();
        for (Peak pk : peaks) {
            if (!pk.isMatched()) negPks.add(pk);
        }
        Collections.shuffle(negPks);
        for(int i = 0; i < limitN; i++) {
            Peak pk = negPks.get(i);
            u_intensity[i] = pk.getNormIntensity();
            u_distance[i]  = pk.getDist();
        }
        negPks.clear();
        negPks = null;

    }
	

	public void calcMean() {
		
		double N, sum;
		
		// Intensity
		sum = 0;
		N = (double) b_intensity.length;
		for(double d : b_intensity) sum += d;
		mu_int_B = (sum / N);
		
		sum = 0;
		N = (double) y_intensity.length;
		for(double d : y_intensity) sum += d;
		mu_int_Y = (sum / N);
		
		sum = 0;
		N = (double) u_intensity.length;
		for(double d : u_intensity) sum += d;
		mu_int_U = (sum / N);
		


		// Distance
		sum = 0;
		N = (double) b_distance.length;
		for(double d : b_distance) sum += d;
		mu_dist_B = (sum / N);

        sum = 0;
        N = (double) y_distance.length;
        for(double d : y_distance) sum += d;
        mu_dist_Y = (sum / N);
		
		sum = 0;
		N = (double) u_distance.length;
		for(double d : u_distance) sum += d;
		mu_dist_U = 0;

	}

	public void calcVar() {
		double N, v;
		
		// Intensity
		v = 0;
		N = (double) b_intensity.length - 1.0;
		for(double d : b_intensity) {
			double x = d - mu_int_B;
			v += Math.pow(x, 2.0);
		}
		var_int_B = (v / N);
		
		v = 0;
		N = (double) y_intensity.length - 1.0;
		for(double d : y_intensity) {
			double x = d - mu_int_Y;
			v += Math.pow(x, 2.0);
		}
		var_int_Y = (v / N);
		
		v = 0;
		N = (double) u_intensity.length - 1.0;
		for(double d : u_intensity) {
			double x = d - mu_int_U;
			v += Math.pow(x, 2.0);
		}
		var_int_U = (v / N);

		
		// Distance
		v = 0;
		N = (double) b_distance.length - 1.0;
		for(double d : b_distance) {
			double x = d - mu_dist_B;
			v += Math.pow(x, 2.0);
		}
		var_dist_B = (v / N);
		
		v = 0;
		N = (double) y_distance.length - 1.0;
		for(double d : y_distance) {
			double x = d - mu_dist_Y;
			v += Math.pow(x, 2.0);
		}
		var_dist_Y = (v / N);
		
		v = 0;
		N = (double) u_distance.length - 1.0;
		for(double d : u_distance) {
			double x = d - mu_dist_U;
			v += Math.pow(x, 2.0);
		}
		var_dist_U = (v / N);
		
		
		var_dist_B *= CID_ADJUST;
		var_dist_Y *= CID_ADJUST;

		
	}
	
	
	public void printSummaryStats() {
		
		System.err.print("+" + chargeState + ": " + numPSM + " PSMs for modeling.\n");
		System.err.print("-----------------------------------------------------\n");
		
		System.err.print(
			"+" + chargeState + "\tb-ions Intensity (mu, sigma): (" +
        MathFunctions.roundDouble(mu_int_B, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_int_B), 4) + ") N = " +
			b_intensity.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\ty-ions Intensity (mu, sigma): (" +
        MathFunctions.roundDouble(mu_int_Y, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_int_Y), 4) + ") N = " +
			y_intensity.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tNoise Intensity (mu, sigma): (" +
        MathFunctions.roundDouble(mu_int_U, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_int_U), 4) + ") N = " +
			u_intensity.length + "\n\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tb-ions m/z Accuracy (mu, sigma): (" +
        MathFunctions.roundDouble(mu_dist_B, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_dist_B), 4) + ") N = " +
			b_distance.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\ty-ions m/z Accuracy (mu, sigma): (" +
        MathFunctions.roundDouble(mu_dist_Y, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_dist_Y), 4) + ") N = " +
			y_distance.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tNoise Distance (mu, sigma): (" +
        MathFunctions.roundDouble(mu_dist_U, 4) + ", " +
        MathFunctions.roundDouble(Math.sqrt(var_dist_U), 4) + ") N = " +
			u_distance.length + "\n\n"
		);
	}


    /***************
     * Function clears out arrays to make more memory available. You call this function AFTER you have all the
     * modeling parameters recorded.
     * @param
     */
    public void clearArrays() {

        b_distance = null;
        b_intensity = null;
        y_distance = null;
        y_intensity = null;
        u_distance = null;
        u_intensity = null;
    }


    /******************
     * Function writes the peaks used for modeling to disk
     */
    public void writeModelPks() throws IOException {
        File debugF = new File("debug_model_pks_CID.txt");
        FileWriter fw;
        BufferedWriter bw;
        String line;
        double dist, normI;

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


      for (double v : b_intensity) {
        normI = MathFunctions.roundDouble(v, 4);
        line = chargeState + "\tyi\t" + normI + "\n";
        bw.write(line);
      }

      for (double v4 : y_intensity) {
        normI = MathFunctions.roundDouble(v4, 4);
        line = chargeState + "\tyi\t" + normI + "\n";
        bw.write(line);
      }

      for (double v3 : u_intensity) {
        normI = MathFunctions.roundDouble(v3, 4);
        line = chargeState + "\tni\t" + normI + "\n";
        bw.write(line);
      }

      for (double v2 : b_distance) {
        dist = MathFunctions.roundDouble(v2, 4);
        line = chargeState + "\tbd\t" + dist + "\n";
        bw.write(line);
      }

      for (double v1 : y_distance) {
        dist = MathFunctions.roundDouble(v1, 4);
        line = chargeState + "\tyd\t" + dist + "\n";
        bw.write(line);
      }

      for (double v : u_distance) {
        dist = MathFunctions.roundDouble(v, 4);
        line = chargeState + "\tnd\t" + dist + "\n";
        bw.write(line);
      }

        bw.close();
    }
}
