/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lucXor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author dfermin
 */
public class ModelData_CID {
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
	
	public ModelData_CID(int z, ArrayList<PeakClass> peaks) {
        chargeState = z;
        numPSM = 0;

        int Nb = 0, Ny = 0, Nu = 0;

        for (PeakClass pk : peaks) {
            if (pk.matched) {
                if (pk.matchedIonStr.startsWith("b")) {
                    Nb++;
                }
                if (pk.matchedIonStr.startsWith("y")) {
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

        for (PeakClass pk : peaks) {
            if (pk.matched) {
                if (pk.matchedIonStr.startsWith("b")) {
                    b_intensity[Nb] = pk.norm_intensity;
                    b_distance[Nb] = pk.dist;
                    Nb++;
                }

                if (pk.matchedIonStr.startsWith("y")) {
                    y_intensity[Ny] = pk.norm_intensity;
                    y_distance[Ny] = pk.dist;
                    Ny++;
                }
            }
        }


        // We will limit the size of the negative distribution to speed things up
        int limitN = (Nb + Ny);
        if(limitN < constants.MIN_NUM_NEG_PKS) limitN += constants.MIN_NUM_NEG_PKS;

        if(limitN > Nu) limitN = Nu; // prevents segfault on insufficient data for modeling
        u_intensity = new double[ limitN ];
        u_distance = new double[ limitN ];

        ArrayList<PeakClass> negPks = new ArrayList<>();
        for (PeakClass pk : peaks) {
            if (!pk.matched) negPks.add(pk);
        }
        Collections.shuffle(negPks);
        for(int i = 0; i < limitN; i++) {
            PeakClass pk = negPks.get(i);
            u_intensity[i] = pk.norm_intensity;
            u_distance[i]  = pk.dist;
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
			globals.round_dbl(mu_int_B, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_int_B), 4) + ") N = " +
			b_intensity.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\ty-ions Intensity (mu, sigma): (" +
			globals.round_dbl(mu_int_Y, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_int_Y), 4) + ") N = " +
			y_intensity.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tNoise Intensity (mu, sigma): (" +
			globals.round_dbl(mu_int_U, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_int_U), 4) + ") N = " +
			u_intensity.length + "\n\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tb-ions m/z Accuracy (mu, sigma): (" +
			globals.round_dbl(mu_dist_B, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_dist_B), 4) + ") N = " +
			b_distance.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\ty-ions m/z Accuracy (mu, sigma): (" +
			globals.round_dbl(mu_dist_Y, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_dist_Y), 4) + ") N = " +
			y_distance.length + "\n"
		);
		
		System.err.print(
			"+" + chargeState + "\tNoise Distance (mu, sigma): (" +
			globals.round_dbl(mu_dist_U, 4) + ", " +
			globals.round_dbl(Math.sqrt(var_dist_U), 4) + ") N = " +
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
        FileWriter fw = null;
        BufferedWriter bw = null;
        String line;
        double mz, dist, relI, normI;

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


        for(int b = 0; b < b_intensity.length; b++) {
            normI = globals.round_dbl(b_intensity[b], 4);
            line = Integer.toString(chargeState) + "\tyi\t" +
                    Double.toString(normI) + "\n";
            bw.write(line);
        }

        for(int y = 0; y < y_intensity.length; y++) {
            normI = globals.round_dbl(y_intensity[y], 4);
            line = Integer.toString(chargeState) + "\tyi\t" +
                    Double.toString(normI) + "\n";
            bw.write(line);
        }

        for(int n = 0; n < u_intensity.length; n++) {
            normI = globals.round_dbl(u_intensity[n], 4);
            line = Integer.toString(chargeState) + "\tni\t" +
                    Double.toString(normI) + "\n";
            bw.write(line);
        }

        for(int b = 0; b < b_distance.length; b++) {
            dist = globals.round_dbl(b_distance[b], 4);
            line = Integer.toString(chargeState) + "\tbd\t" +
                    Double.toString(dist) + "\n";
            bw.write(line);
        }

        for(int y = 0; y < y_distance.length; y++) {
            dist = globals.round_dbl(y_distance[y], 4);
            line = Integer.toString(chargeState) + "\tyd\t" +
                    Double.toString(dist) + "\n";
            bw.write(line);
        }

        for(int n = 0; n < u_distance.length; n++) {
            dist = globals.round_dbl(u_distance[n], 4);
            line = Integer.toString(chargeState) + "\tnd\t" +
                    Double.toString(dist) + "\n";
            bw.write(line);
        }

        bw.close();
    }
}
