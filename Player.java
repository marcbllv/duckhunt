import java.util.Vector;
import java.util.Random;
import java.util.ArrayList;

class Player {

    private static final int NBR_HMM = 1;
    private static final int NBR_OBS = 1000;

    public ArrayList<HMM> roundHmms;
    public HMM[][] speciesHmm = new HMM[6][NBR_HMM];
    public int mvmt = 0;
    public boolean[] speciesSeen = new boolean[] {false, false, false, false, false, false};

    public int[][] obsSpecies = new int[6][NBR_OBS];
    public int[] obsSpeciesSize = new int[] {0, 0, 0, 0, 0, 0};

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

            // Init new HMM randomly for each bird :
            for(int i = 0 ; i < pState.getNumBirds() ; i++) {
                HMM h = randHMM();
                hmms.add(h);
            }

            roundHmms = new ArrayList<HMM>(hmms);
        }

        mvmt++;

        // This line chooses not to shoot.
        // return new Action(0, Constants.MOVE_RIGHT);

        // This line would predict that bird 0 will move right and shoot at it.
        return new Action(-1, -1);
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
            Random randomGenerator = new Random();
            for (int i = 0; i < pState.getNumBirds(); ++i) {
                lGuess[i] = randomGenerator.nextInt(6);
            }
        } else {
            for (int i = 0; i < pState.getNumBirds(); ++i) {

                if(pState.getBird(i).isDead()) {
                    continue;
                }
                int[] bSeq = new int[mvmt];
                for(int j = 0 ; j < bSeq.length ; j++) {
                    bSeq[j] = pState.getBird(i).getObservation(j);
                }
                Player.smoothStoch(roundHmms.get(i).A);
                Player.smoothStoch(roundHmms.get(i).B);
                roundHmms.get(i).estimateModel(bSeq);

                HMM birdHmm = roundHmms.get(i);
                double[] proba = new double[NBR_HMM];
                double pr;
                double maxProba = -1;
                int mostProbSp = -1;

                for(int k = 0 ; k < 6 ; k++) {
                    if(speciesSeen[k]) {
                        // Getting all the probabilities of all hmm from this species
                        for(int l = 0 ; l < NBR_HMM ; l++) {
                            proba[l] = (speciesHmm[k][l]).estimateProbabilityOfEmissionSequence(bSeq);
                        }
                        // Getting the best probas using : mean -> 97 on kattis
                        pr = 0;
                        for(int l = 0 ; l < NBR_HMM ; l++) {
                            pr += proba[l];
                        }
                        pr /= NBR_HMM;

                        // /!\
                        // NB : NBR_NMM must be > 3 in the following selectors !!
                        // /!\
                        
                        /*
                        // Sorting probas (decr) in array, needed for median & mean best 3
                        for(int idx1 = 0 ; idx1 < NBR_HMM - 1 ; idx1++) {
                            double max = proba[idx1];
                            int maxIdx = idx1;
                            for(int idx2 = idx1 + 1 ; idx2 < NBR_HMM ; idx2++) {
                                if(max < proba[idx2]) {
                                    max = proba[idx2];
                                    maxIdx = idx2;
                                }
                            }
                            double c = proba[maxIdx];
                            proba[maxIdx] = proba[idx1];
                            proba[idx1] = c;
                        }
                        // Getting the best proba using : best median -> 77 on kattis
                        //pr = proba[NBR_HMM / 2];

                        // Getting the best probas using : mean -> 97 on kattis
                        //pr = 0;
                        //for(int l = 0 ; l < NBR_HMM ; l++) {
                        //    pr += proba[l];
                        //}
                        //pr /= NBR_HMM;
                        
                        // Getting the best probas using : mean among 3 bests -> 75 on kattis
                        pr = 0;
                        for(int l = 0 ; l < 3 ; l++) {
                            pr += proba[l];
                        }
                        pr /= NBR_HMM;
                        */

                        if(pr > maxProba) {
                            maxProba = pr;
                            mostProbSp = k;
                        }
                    }
                }

                lGuess[i] = mostProbSp;
            }
        }
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

        for(int i = 0 ; i < pState.getNumBirds() ; i++) {
            if(pState.getBird(i).isDead()) {
                continue;
            }
            // If new species : training new HMMs and store them
            if(!speciesSeen[pSpecies[i]]) {
                speciesSeen[pSpecies[i]] = true;
                for(int l = 0 ; l < NBR_HMM ; l++) {
                    speciesHmm[pSpecies[i]][l] = randHMM();
                }
            }

            // Adding the observations to observations collection
            int s = obsSpeciesSize[pSpecies[i]]; // current number of obs for this species
            for(int j = 0 ; j < mvmt ; j++) {
                if(s + j < NBR_OBS) {
                    obsSpecies[pSpecies[i]][s + j] = pState.getBird(i).getObservation(j);
                }
            }
            s += mvmt;
            if(s < NBR_OBS) {
                obsSpeciesSize[pSpecies[i]] += s;
            } else {
                obsSpeciesSize[pSpecies[i]] = NBR_OBS;
            }

            // Clearing hmm with smoothing function & reestimate model
            for(int l = 0; l < NBR_HMM ; l++) {
                Player.smoothStoch(speciesHmm[pSpecies[i]][l].A);
                Player.smoothStoch(speciesHmm[pSpecies[i]][l].B);

                speciesHmm[pSpecies[i]][l].estimateModel( obsSpecies[pSpecies[i]] );
            }

        }
        mvmt = 0;
    }

    public static void smoothStoch(double[][] M) {
        double epsilon = 0.00000001;
        for(int i = 0 ; i < M.length ; i++) {
            for(int j = 0 ; j < M[0].length ; j++) {
                M[i][j] += epsilon;
                M[i][j] /= (1 + M[0].length * epsilon);
            }
        }
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
                M[i][j] = (double)randomGenerator.nextInt(100) + 1;
                if(i == j) {
                    M[i][j] += 100.0;
                }
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

    public static HMM randHMM() {
        double[][] A_init = randStochMat(5,5);
        double[][] B_init = randStochMat(5,9);
        double[][] pi_init_2D = randStochMat(1,5);
        double[] pi_init = new double[5];
        for(int j = 0 ; j < 5 ; j++) {
            pi_init[j] = pi_init_2D[0][j];
        }

        return new HMM(A_init, B_init, pi_init);
    }

    public static final Action cDontShoot = new Action(-1, -1);

}
