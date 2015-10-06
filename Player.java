import java.util.Vector;
import java.util.Random;
import java.util.ArrayList;

class Player {

    public static ArrayList<HMM> roundHmms;
    public static HMM[] speciesHmm = new HMM[6];
    public static int mvmt = 0;
    public static boolean[] speciesSeen = new boolean[] {false, false, false, false, false, false};

    public Player() {
    }

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to get the best action.
         * This skeleton never shoots.
         */

        if(mvmt == 0) {
            if(pState.getRound() > 0)
                roundHmms.clear();

            ArrayList<HMM> hmms = new ArrayList<HMM>();

            // Init HMM for each bird :
            for(int i = 0 ; i < pState.getNumBirds() ; i++) {
                double[][] A_init = randStochMat(5,5);
                double[][] B_init = randStochMat(5,9);
                double[][] pi_init_2D = randStochMat(1,5);
                double[] pi_init = new double[5];
                for(int j = 0 ; j < 5 ; j++) {
                    pi_init[j] = pi_init_2D[0][j];
                }

                HMM h = new HMM(A_init, B_init, pi_init);
                hmms.add(h);
            }

            roundHmms = new ArrayList<HMM>(hmms);
        }

        if(mvmt > 99) {
            for(int i = 0 ; i < pState.getNumBirds() ; i++) {
                // Get the sequence of observations
                int[] bSeq = new int[mvmt + 1];
                for(int j = 0 ; j < bSeq.length ; j++) {
                    bSeq[j] = pState.getBird(i).getObservation(j);
                }
                roundHmms.get(i).estimateModel(bSeq);
                Player.smoothStoch(roundHmms.get(i).A);
                Player.smoothStoch(roundHmms.get(i).B);
            }
        }
        // This line chooses not to shoot.

        mvmt++;
        return cDontShoot;

        // This line would predict that bird 0 will move right and shoot at it.
        // return Action(0, MOVE_RIGHT);
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        int[] lGuess = new int[pState.getNumBirds()];

        if(pState.getRound() == 0) {
            for (int i = 0; i < pState.getNumBirds(); ++i) {
                lGuess[i] = Constants.SPECIES_PIGEON;
            }
        } else {
            for (int i = 0; i < pState.getNumBirds(); ++i) {

                int[] bSeq = new int[99];
                for(int j = 0 ; j < bSeq.length ; j++) {
                    bSeq[j] = pState.getBird(i).getObservation(j);
                }

                HMM birdHmm = Player.roundHmms.get(i);
                double proba = -1, maxProba = -1;
                int mostProbSp = -1;

                for(int k = 0 ; k < 6 ; k++) {
                    if(speciesSeen[k]) {
                        proba = (speciesHmm[k]).estimateProbabilityOfEmissionSequence(bSeq);
                        if(proba > maxProba) {
                            maxProba = proba;
                            mostProbSp = k;
                        }
                    }
                }

                lGuess[i] = mostProbSp;
            }
        }
        System.err.println();
        System.err.print("Round : " + pState.getRound() + ", guesses: ");
        for(int i = 0 ; i < lGuess.length ; i++) {
            System.err.print(lGuess[i] + " ");
        }
        System.err.println();
        return lGuess;
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {

        System.err.print("Reveal:             ");
        for(int i = 0 ; i < pState.getNumBirds() ; i++) {
            System.err.print(pSpecies[i] + " ");

            if(!speciesSeen[pSpecies[i]]) {
                speciesSeen[pSpecies[i]] = true;
                speciesHmm[pSpecies[i]] = roundHmms.get(i);
            } else {
                // Reestimate model
                int[] bSeq = new int[99];
                for(int j = 0 ; j < bSeq.length ; j++) {
                    bSeq[j] = pState.getBird(i).getObservation(j);
                }

                // Clearing hmm with smoothing function & reestimate model
                Player.smoothStoch(speciesHmm[pSpecies[i]].A);
                Player.smoothStoch(speciesHmm[pSpecies[i]].B);

                speciesHmm[pSpecies[i]].estimateModel(bSeq);
            }
        }
        System.err.println();

        mvmt = 0;
    }

    public static void smoothStoch(double[][] M) {
        double epsilon = 0.000001;
        for(int i = 0 ; i < M.length ; i++) {
            for(int j = 0 ; j < M[0].length ; j++) {
                M[i][j] += epsilon;
                M[i][j] /= (1 + M[0].length * epsilon);
            }
        }
    }



    // This function looks for the best permutation to make A1 and A2 closer
    public static Permut permut(double[][] A1, double[][] A2) {

        int[] p = new int[] {0,1,2,3,4};
        int[] minP = new int[] {0,1,2,3,4};
        double err = 0;
        double minErr = 10000;
        
        // Lot of loops to find all the permutations
        for(int i0 = 0 ; i0 < 5 ; i0++) {
            p[0] = i0;
        for(int i1 = 0 ; i1 < 5 ; i1++) {
            if(p[0] == i1)
                continue;
            p[1] = i1;
        for(int i2 = 0 ; i2 < 5 ; i2++) {
            if(p[0] == i2 || p[1] == i2)
                continue;
            p[2] = i2;
        for(int i3 = 0 ; i3 < 5 ; i3++) {
            if(p[0] == i3 || p[1] == i3 || p[2] == i3)
                continue;
            p[3] = i3;
        for(int i4 = 0 ; i4 < 5 ; i4++) {
            if(p[0] == i4 || p[1] == i4 || p[2] == i4 || p[3] == i4)
                continue;
            p[4] = i4;

            // Now we have the permutation, we can compute the error
            err = 0;
            for(int i = 0 ; i < 5 ; i++) {
                for(int j = 0 ; j < 5 ; j++) {
                    err += (A1[i][j] - A2[p[i]][p[j]]) * (A1[i][j] - A2[p[i]][p[j]]);
                }
            }

            if(err < minErr) {
                minErr = err;
                for(int i = 0 ; i < 5 ; i++) {
                    minP[i] = p[i];
                }
            }
        }
        }
        }
        }
        }

        // Now let's apply the permutation on A2
        double[][] A2_permut = new double[5][5];
        for(int i = 0 ; i < 5 ; i++) {
            for(int j = 0 ; j < 5 ; j++) {
                A2_permut[i][j] = A2[minP[i]][minP[j]];
            }
        }

        Permut permut = new Permut(A2_permut, minErr);
        return permut;
    }

    // Compute the mean of A1 & A2
    public static double[][] meanMat(double[][] A1, double[][] A2) {
        double[][] M = new double[5][5];
        for(int i = 0 ; i < 5 ; i++) {
            for(int j = 0 ; j < 5 ; j++) {
                M[i][j] = (A1[i][j] + A2[i][j]) / 2;
            }
        }

        return M;
    }

    // Just print the matrix
    public static void printMat(double[][] M) {
        double sum = 0.0;
        int nL = M.length;
        int nR = M[0].length;

        for(int i = 0 ; i < nL ; i++) {
            sum = 0.0;
            for(int j = 0 ; j < nR ; j++) {
                System.err.printf("%.4f ", M[i][j]);
                sum += M[i][j];
            }
            System.err.printf(" = %.8f\n", sum);
        }
        System.err.println();
    }

    // Generate random stochastic matrices
    // with nL lines and nR rows
    public static double[][] randStochMat(int nL, int nR) {
        if(nL < 1) {
            nL = 1;
        }
        if(nR < 1) {
            nR = 1;
        }

        Random randomGenerator = new Random();

        double[][] M = new double[nL][nR];
        int sumNorm = 0;
        for(int i = 0 ; i < nL ; i++) {
            for(int j = 0 ; j < nR ; j++) {
                M[i][j] = (double)randomGenerator.nextInt(100);
                sumNorm += M[i][j];
            }
            for(int j = 0 ; j < nR ; j++) {
                M[i][j] = M[i][j] / (double)sumNorm;
            }
        }

        return M;
    }

    // This method adjusts M to make it a stochastic matix
    // It divides all lines by its sum so that they sum to 1
    // Assuming all coefficients are positives and the sum is non null
    public static void makeStochastic(double[][] M) {
        double sum;
        for(int i = 0 ; i < M.length ; i++) {
            sum = 0.0;
            for(int j = 0 ; j < M[0].length ; j++) {
                sum += M[i][j];
            }
            for(int j = 0 ; j < M[0].length ; j++) {
                M[i][j] /= sum;
            }
        }
    }

    public static final Action cDontShoot = new Action(-1, -1);

    public static double[][][] As = new double[6][5][5];
    public static double[][][] Bs = new double[6][5][9];
    public static double[][] pis = new double[6][5];

    private static class Permut {
        double[][] M;
        double err;

        Permut(double[][] M, double err) {
            this.M = M;
            this.err = err;
        }
    }
}
