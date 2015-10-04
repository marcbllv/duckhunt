import java.util.Vector;
import java.util.Random;

class Player {

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

        // This line chooses not to shoot.
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
            System.err.println("Round " + pState.getRound() + " " + pState.getNumBirds() + " birds");
            for(int i = 0 ; i < pState.getNumBirds() ; i++) {
                Bird b = pState.getBird(i);

                // Get the sequence of observations
                int[] bSeq = new int[99];
                for(int j = 0 ; j < bSeq.length ; j++) {
                    bSeq[j] = b.getObservation(j);
                }

                double[][] A_init = randStochMat(5,5);
                double[][] B_init = randStochMat(5,9);
                double[][] pi_init_2D = randStochMat(1,5);
                double[] pi_init = new double[5];
                for(int j = 0 ; j < 5 ; j++) {
                    pi_init[j] = pi_init_2D[0][j];
                }

                HMM h = new HMM(A_init, B_init, pi_init);

                // Train hmm
                h.estimateModel(bSeq);

                // Comparing matrices with all species to get the closer
                double minErr = 10000;
                Double err = new Double(0.0);
                int probableSp = -1;
                for(int sp = 0 ; sp < 7 ; sp++) {
                    Player.permut(As[sp], h.A, err);
                    if(err < minErr) {
                        minErr = err;
                        probableSp = sp;
                    }
                }

                lGuess[i] = probableSp;
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

        System.err.print("Reveal part - Species seen : ");
        for(int i = 0 ; i < pSpecies.length ; i++) {
            System.err.print(pSpecies[i] + " ");
        }
        System.err.println();

        for(int i = 0 ; i < pState.getNumBirds() ; i++) {
            Bird b = pState.getBird(i);

            // Get the sequence of observations
            int[] bSeq = new int[99];
            for(int j = 0 ; j < bSeq.length ; j++) {
                bSeq[j] = b.getObservation(j);
            }

            // Init Matrices & Hmm
            double[][] A_init;
            double[][] B_init;
            double[] pi_init;
            if(!speciesSeen[pSpecies[i]]) {
                // If we havn't seen this species before : random init
                A_init = randStochMat(5,5);
                B_init = randStochMat(5,9);
                double[][] pi_init_2D = randStochMat(1,5);
                pi_init = new double[5];
                for(int j = 0 ; j < 5 ; j++) {
                    pi_init[j] = pi_init_2D[0][j];
                }
            } else {
                // Else : we can start with the species matrix
                A_init = As[pSpecies[i]];
                B_init = Bs[pSpecies[i]];
                pi_init = pis[pSpecies[i]];
            }

            HMM h = new HMM(A_init, B_init, pi_init);

            // Train hmm
            h.estimateModel(bSeq);

            // First time we see this species
            if(!speciesSeen[pSpecies[i]]) {
                speciesSeen[pSpecies[i]] = true;

                As[pSpecies[i]] = h.A;
                Bs[pSpecies[i]] = h.B;
                pis[pSpecies[i]] = h.pi;
            } else {
                // Apply most probable permutation to h.A
                double[][] A_p = Player.permut(As[pSpecies[i]], h.A);
                // Saving new matrix for the current species
                // TODO : not that good to make the mean each time, change that line:
                As[pSpecies[i]] = Player.meanMat(As[pSpecies[i]], A_p);
            }
        }

        for(int i = 0 ; i < 7 ; i++) {
            if(speciesSeen[i]) {
                System.err.println("Species " + i + " A matrix:");
                Player.printMat(As[i]);
            }
        }
        System.err.println();
    }

    // This function looks for the best permutation to make A1 and A2 closer
    public static double[][] permut(double[][] A1, double[][] A2) {
        Double d = new Double(0.0);
        return permut(A1, A2, d);
    }

    // Same function but the error |A1 - A2|^2 is stored in the error parameter
    public static double[][] permut(double[][] A1, double[][] A2, Double error) {

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

        error = new Double(minErr);

        // Now let's apply the permutation on A2
        double[][] A2_permut = new double[5][5];
        for(int i = 0 ; i < 5 ; i++) {
            for(int j = 0 ; j < 5 ; j++) {
                A2_permut[i][j] = A2[minP[i]][minP[j]];
            }
        }
        return A2_permut;
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
        for(int i = 0 ; i < 5 ; i++) {
            for(int j = 0 ; j < 5 ; j++) {
                System.err.printf("%.4f ", M[i][j]);
            }
            System.err.println();
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

    public static final Action cDontShoot = new Action(-1, -1);

    public static boolean[] speciesSeen = new boolean[] 
        {false, false, false, false, false, false, false};
    public static double[][][] As = new double[7][5][5];
    public static double[][][] Bs = new double[7][5][9];
    public static double[][] pis = new double[7][5];
}
