package lucxor;


import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by dfermin on 3/31/14.
 */
public class Spectrum {

    public int N; // number of peaks
	public int maxI_index; // the max. intensity peaks index in this.raw_intensity[]
    public double maxI; // max intensity;


    public double[] mz = null;
    public double[] raw_intensity = null;
	public double[] rel_intensity = null;
	public double[] norm_intensity = null; 


    /*
     * Default constructor we use for Spectrum
     * @param mz
     * @param intensities
     */
    public Spectrum(double[] mz, double[] intensities) {
        N = mz.length;
        maxI = 0;
        maxI_index = 0;

        // This should never happen but just to avoid a segfault we put this conditional in
        if(mz.length != intensities.length)
            N = Math.min(mz.length, intensities.length);

        TIntArrayList candPks = new TIntArrayList(N);
        for(int i = 0; i < N; i++) {
            if(intensities[i] == 0) continue; // skip zero intensity peaks
            candPks.add(i);
        }

        N = candPks.size(); // *new* peak list size

        this.mz = new double[ N ];
        this.raw_intensity = new double[ N ];
        this.norm_intensity = new double[ N ];
        this.rel_intensity = new double[ N ];

        // record passed values into arrays
        for(int i = 0; i < N; i++) {
            int curIdx = candPks.get(i);

            if(intensities[ curIdx ] > maxI) {
                maxI = intensities[ curIdx ];
                maxI_index = curIdx;
            }
            this.mz[i] = mz[ curIdx ];
            this.raw_intensity[i] = intensities[ curIdx ];

        }

        candPks.clear();
        candPks = null;

        // compute relative intensity and prep remaining arrays
        for(int i = 0; i < N; i++) {
            this.rel_intensity[i] = (this.raw_intensity[i] / maxI) * 100.0;
        }
    }


    public Spectrum() {} // empty constructor


    /*
     * Function calculates the relative intensity of the spectrum based upon the value of
     * maxI. You only need to call this function when the value of 'maxI' has changed.
     */
    public void calcRelativeIntensity() {
        for(int i = 0; i < N; i++) {
            this.rel_intensity[i] = (this.raw_intensity[i] / maxI) * 100.0;
        }
    }


    /*
     * Function computes the log(rel_intensity/median_intensity)
     */
    public void medianNormalizeSpectra() {
        int mid = N / 2;
        double medianI = 0d;

        // Need to sort the peak intensities from low to high
        ArrayList<Double> pksI = new ArrayList<>(N);
        for(int i = 0; i < N; i++) pksI.add(rel_intensity[i]);

        Collections.sort(pksI); // intensities sorted from low to high

        if(N % 2  == 0) { // even number of elements
            int a = mid - 1;

            double aI = pksI.get(a);
            double bI = pksI.get(mid);
            medianI = (aI + bI) / 2.0;
        }
        else { // odd number of elements
            medianI = pksI.get(mid);
        }

        for(int i = 0; i < N; i++) {
            double d = rel_intensity[i] / medianI;
            norm_intensity[i] = FastMath.log(d);
        }
    }

    /*
     * Get a particular peak value mz or intensity
     */
    public double getPeak(int idx, int dt) {
        double ret = 0;

        switch(dt) {
            case 1:
                ret = this.mz[idx];
                break;
            case 2:
                ret = this.raw_intensity[idx];
                break;
            case 3:
                ret = this.rel_intensity[idx];
                break;
            case 4:
                ret = this.norm_intensity[idx];
                break;
        }

        return ret;
    }

    /*
     * Returns true if the mz variable is empty or unset
     * @return
     */
    public boolean isEmpty() {
        if(null == this.mz) return true;
        return this.mz.length == 0;
    }


    /**********
     * This function returns the index for the element corresponding to the given mz value
     */
    public int findIndexByMZ(double passed_mz) {
        int ret = -1;

        for(int i = 0; i < N; i++) {
            if(this.mz[i] == passed_mz) {
                ret = i;
                break;
            }
        }
        return ret;
    }

    /**
     * Get the peak information
     * @param idx index of the Peak
     * @return new Peak
     */
    public Peak getPeakClassInstance(int idx) {
        Peak ret = new Peak();
        ret.setMz(this.mz[ idx ]);
        ret.setRawIntensity(this.raw_intensity[ idx ]);
        ret.setRelIntensity(this.rel_intensity[ idx ]);
        ret.setNormIntensity(this.norm_intensity[ idx ]);
        return ret;
    }
}
