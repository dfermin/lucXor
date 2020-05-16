package lucxor;

import java.io.IOException;

/**
 * Created by dfermin on 2/17/14.
 */



public class ScoringWorkerThread implements Runnable {

    private PSM curPSM; // the PSM to be scored
    private int jobIdx; // index of this PSM in the Globals.psmList
    private int runMode_;

    // Default constructor for this class
    public ScoringWorkerThread(PSM externalPSM, int rn, int i) {
        this.curPSM = externalPSM;
        this.jobIdx = i;
        this.runMode_ = rn;
    }

    @Override
    public void run() {
        try {
            synchronized (this) {
                this.curPSM.generatePermutations(this.runMode_);
                this.curPSM.scorePermutations();

                if ((this.jobIdx % 100) == 0) System.err.print(this.jobIdx + " ");
                if ((this.jobIdx % 1000) == 0) System.err.print("\n");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
