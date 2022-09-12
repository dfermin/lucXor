package lucxor.common;


import gnu.trove.list.array.TIntArrayList;
import lombok.Builder;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.Collections;

/**
 * Created by dfermin on 3/31/14.
 */
@Builder
public class Spectrum {

    public int N; // number of peaks
	public int maxI_index; // the max. intensity peaks index in this.raw_intensity[]
    public double maxI; // max intensity;


    public double[] mz = null;
    public double[] raw_intensity = null;
	public double[] rel_intensity = null;
	private double[] norm_intensity = null;


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

        // compute relative intensity and prep remaining arrays
        for(int i = 0; i < N; i++) {
            this.rel_intensity[i] = (this.raw_intensity[i] / maxI) * 100.0;
        }
    }

    public Spectrum(int n, int maxI_index, double maxI, double[] mz, double[] raw_intensity, double[] rel_intensity, double[] norm_intensity) {
        N = n;
        this.maxI_index = maxI_index;
        this.maxI = maxI;
        this.mz = mz;
        this.raw_intensity = raw_intensity;
        this.rel_intensity = rel_intensity;
        this.norm_intensity = norm_intensity;
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
        double medianI;

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
