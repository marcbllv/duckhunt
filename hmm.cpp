#include "hmm.hpp"
#include <iostream>
#include <cmath>
#include <limits>

namespace ducks {

Hmm::Hmm(int numberofstates, int numberofemissions) {
    /**
     * This constructor just sets everything to zero (why is this wrong?).
     * Write a clever way to initialize the HMM!
     */
    numberofstates_ = numberofstates;
    numberofemissions_ = numberofemissions;

    A_ = new double*[numberofstates_];
    for (int i = 0; i < numberofstates_; ++i) {
		A_[i] = new double[numberofstates_];
	}
    B_ = new double*[numberofstates_];
    for (int i = 0; i < numberofstates_; ++i) {
		B_[i] = new double[numberofemissions_];
	}
    pi_ = new double[numberofstates_];

    for (int i = 0; i < numberofstates_; ++i) {
      for (int j = 0; j < numberofstates_; ++j) {
        if(j==i){A_[i][j] = 2.0/(numberofstates_+1);}
        else{A_[i][j] = (1.0/(numberofstates_+1));};
      }
    }

    for (int i = 0; i < numberofstates_; ++i) {
      for (int j = 0; j < numberofemissions_; ++j) {
        if(j==i){B_[i][j] = 2.0/(numberofemissions_+1);}
        else{B_[i][j] = (1.0/(numberofemissions_+1));};
      }
    }

    for (int i = 0; i < numberofstates_; ++i) {
      pi_[i] = 0.0;
    }
    pi_[0]=1.0;
}

Hmm::Hmm(int numberofstates, int numberofemissions, double** A, double** B, double* pi) {
    numberofstates_ = numberofstates;
    numberofemissions_ = numberofemissions;

    A_ = new double*[numberofstates_];
    for (int i = 0; i < numberofstates_; ++i) {
		A_[i] = new double[numberofstates_];
	}
    B_ = new double*[numberofstates_];
    for (int i = 0; i < numberofstates_; ++i) {
		B_[i] = new double[numberofemissions_];
	}
    pi_ = new double[numberofstates_];

    for (int i = 0; i < numberofstates_; ++i) {
      for (int j = 0; j < numberofstates_; ++j) {
        A_[i][j] = A[i][j];
      }
    }

    for (int i = 0; i < numberofstates_; ++i) {
      for (int j = 0; j < numberofemissions_; ++j) {
        B_[i][j] = B[i][j];
      }
    }

    for (int i = 0; i < numberofstates_; ++i) {
      pi_[i] = pi[i];
    }
}

void Hmm::info()
{
    std::cerr << "Vecteur :";
    for(int i=0;i<numberofstates_;++i)
    {
        std::cerr << pi_[i] << " ";
    }
    std::cerr << std::endl;
}
  /**
   * Estimates the probability distribution of the next emission, given the
   * current state probability distribution.
   *
   * Note that this method solves the preparatory exercise HMM1.
   *
   */
  double* Hmm::estimateProbabilityDistributionOfNextEmission(
      double* currentstateprobabilitydistribution){
    double* probabilityofmoves = new double[numberofemissions_];
    for(int i = 0; i < numberofstates_; i++){
      for(int j = 0; j < numberofstates_; j++){
        for(int k = 0; k < numberofemissions_; k++){
          probabilityofmoves[k] +=
              currentstateprobabilitydistribution[j]*A_[j][i]*B_[i][k];
        }
      }
    }
    return probabilityofmoves;
  }

  /**
   * Estimates the probability of a sequence of observed emissions, assuming
   * this HMM.
   *
   * Note that this method solves the preparatory exercise HMM2.
   */
  double Hmm::estimateProbabilityOfEmissionSequence(int* O, int numberofobservations){
    double probability = 0.0;
    double** alpha = new double*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		alpha[i] = new double[numberofstates_];
	}

    for (int i = 0; i < numberofstates_; ++i) {
      alpha[0][i]=pi_[i]*B_[i][O[0]];
    }

    for (int t = 1; t < numberofobservations; ++t) {
      for (int i = 0; i < numberofstates_; ++i) {
        for (int j = 0; j < numberofstates_; ++j) {
          alpha[t][i] += alpha[(t-1)][j]*A_[j][i]*B_[i][O[t]];
        }
      }
    }

    for (int i = 0; i < numberofstates_; ++i) {
      probability += alpha[(numberofobservations-1)][i];
    }
    return probability;
  }

  /**
   * Estimates the hidden states from which a sequence of emissions were
   * observed.
   *
   * Note that this method solves the preparatory exercise HMM3.
   */
  int* Hmm::estimateStateSequence(int* O, int numberofobservations){

    double** alpha = new double*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		alpha[i] = new double[numberofstates_];
	}
    int** path = new int*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		path[i] = new int[numberofstates_];
	}

    int* finalpath = new int[numberofobservations];
    double maxprobability;

    for (int i = 0; i < numberofstates_; ++i) {
      alpha[0][i]=pi_[i]*B_[i][O[0]];
    }

    for (int t = 1; t < numberofobservations; ++t) {
      for (int i = 0; i < numberofstates_; ++i) {
        maxprobability = 0.0;
        for (int j = 0; j < numberofstates_; ++j) {
          if (maxprobability < alpha[t-1][j]*A_[j][i]){
            maxprobability = alpha[t-1][j]*A_[j][i];
            path[t][i] = j;
          }
        }
        alpha[t][i] = maxprobability*B_[i][O[t]];
      }
    }

    maxprobability = 0.0;
    for (int i = 0; i < numberofstates_; ++i) {
      if (maxprobability < alpha[numberofobservations-1][i]){
        maxprobability = alpha[numberofobservations-1][i];
        finalpath[numberofobservations-1] = i;
      }
    }

    for (int t = numberofobservations-2; t >= 0; --t) {
      finalpath[t] = path[t+1][finalpath[t+1]];
    }
    return finalpath;
  }

/**
   * Re-estimates this HMM from a sequence of observed emissions.
   *
   * Note that this method solves the preparatory exercise HMM4.
   */
  void Hmm::estimateModel(int* O, int numberofobservations){
    double*** xi = new double**[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		xi[i] = new double*[numberofstates_];
		for (int j = 0; j < numberofstates_; ++j) {
			xi[i][j] = new double[numberofstates_];
		}
	}
    double** alpha = new double*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		alpha[i] = new double[numberofstates_];
	}

    double** beta = new double*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		beta[i] = new double[numberofstates_];
	}
    double** gamma = new double*[numberofobservations];
    for (int i = 0; i < numberofobservations; ++i) {
		gamma[i] = new double[numberofstates_];
	}

    double* C = new double[numberofobservations]; // A scaling factor
    double numer; // A temporary variable for holding a numerator
    double denom; // A temporary variable for holding a denominator

    /* Iteration-related stuff */
    double oldlogprob = -std::numeric_limits<double>::max();
    int iters = 0;
    double logprob;
    bool finished = false;

    while(!finished && iters < kMaxIters){
      /* Computation of alpha */

      C[0]=0.0;
      for (int i = 0; i < numberofstates_; ++i) {
        alpha[0][i]=pi_[i]*B_[i][O[0]];
        C[0] += alpha[0][i];
      }

      if (C[0] != 0){
        C[0] = 1.0/C[0];
      }
      for (int i = 0; i < numberofstates_; ++i) {
        alpha[0][i]=C[0]*alpha[0][i];
      }

      for (int t = 1; t < numberofobservations; ++t) {
        C[t]=0.0;
        for (int i = 0; i < numberofstates_; ++i) {
          alpha[t][i] = 0.0;
          for (int j = 0; j < numberofstates_; ++j) {
            alpha[t][i] += alpha[(t-1)][j]*A_[j][i];
          }
          alpha[t][i]=alpha[t][i]*B_[i][O[t]];
          C[t] += alpha[t][i];
        }

        if (C[t] != 0){
          C[t] = 1.0/C[t];
        }
        for (int i = 0; i < numberofstates_; ++i) {
          alpha[t][i]=C[t]*alpha[t][i];
        }
      }

      /* Computation of beta */

      for (int i = 0; i < numberofstates_; ++i) {
        beta[(numberofobservations-1)][i] = C[numberofobservations-1];
      }

      for (int t = numberofobservations-2; t >= 0; --t){
        for (int i = 0; i < numberofstates_; ++i) {
          beta[t][i] = 0;
          for (int j = 0; j < numberofstates_; ++j) {
            beta[t][i] += A_[i][j]*B_[j][O[t+1]]*beta[t+1][j];
          }
          beta[t][i] = C[t]*beta[t][i];
        }
      }

      /* Computation of gamma and xi */

      for (int t = 0; t < numberofobservations-1; ++t) {
        denom = 0;
        for (int i = 0; i < numberofstates_; ++i) {
          for (int j = 0; j < numberofstates_; ++j){
            /* Eq. 37 in Rabiner89 */
            denom += alpha[t][i]*A_[i][j]*B_[j][O[t+1]]*beta[(t+1)][j];
          }
        }
        for (int i = 0; i < numberofstates_; ++i) {
          gamma[t][i] = 0.0;
          for (int j = 0; j < numberofstates_; ++j){
            if (denom != 0){
              xi[t][i][j] = (alpha[t][i]*A_[i][j]*B_[j][O[t+1]]*beta[(t+1)][j])/denom;
            } else {
              xi[t][i][j] = 0.0;
            }
            /* Eq. 38 in Rabiner89 */
            gamma[t][i] += xi[t][i][j];
          }
        }
      }
      /* We must also calculate gamma for the last step. This is given by Eq. 27
       * in Rabiner89. */
      denom = 0;
      for (int i = 0; i < numberofstates_; ++i) {
          denom += alpha[numberofobservations-1][i]*beta[numberofobservations-1][i];
      }
      for (int i = 0; i < numberofstates_; ++i) {
        gamma[numberofobservations-1][i] = 0.0;
        gamma[numberofobservations-1][i] += (alpha[numberofobservations-1][i]*beta[numberofobservations-1][i])/denom;
      }

      /* Re-estimate A,B and pi */

      //Pi
      for (int i = 0; i < numberofstates_; ++i) {
        pi_[i] = gamma[0][i];
      }

      //A
      for (int i = 0; i < numberofstates_; ++i) {
        for (int j = 0; j < numberofstates_; ++j) {
          numer = 0.0;
          denom = 0.0;
          for (int t = 0; t < numberofobservations-1; ++t) {
            numer += xi[t][i][j];
            denom += gamma[t][i];
          }
          if (denom != 0){
            A_[i][j] = numer/denom;
          } else {
            A_[i][j] = 0;
          }

        }
      }

      //B
      for (int i = 0; i < numberofstates_; ++i) {
        for (int j = 0; j < numberofemissions_; ++j) {
          numer = 0.0;
          denom = 0.0;
          for (int t = 0; t < numberofobservations; ++t) {
            if (j == O[t]){
              numer += gamma[t][i];
            }
            denom += gamma[t][i];
          }
          if (denom != 0){
            B_[i][j] = numer/denom;
          } else {
            B_[i][j] = 0;
          }
        }
      }

      /* Compute log probability for model generating observed sequence */

      logprob = 0.0;
      for (int t = 0; t < numberofobservations; ++t) {
        if (C[t] != 0){
          logprob += log(C[t]);
        }
      }
      logprob = -logprob;
      if (logprob > oldlogprob){
        iters += 1;
        oldlogprob = logprob;
      } else {
        finished = true;
      }
    }
  }

  /**
   * Copies a HMM onto this one.
   */
  void Hmm::copyHmm(Hmm* hmm){
    for (int i = 0; i < numberofstates_; ++i){
      for (int j = 0; j < numberofstates_; ++j){
        A_[i][j] = hmm->A_[i][j];
      }
      pi_[i] = hmm->pi_[i];
    }
    for (int i = 0; i < numberofstates_; ++i){
      for (int j = 0; j < numberofemissions_; ++j){
        B_[i][j] = hmm->B_[i][j];
      }
    }
  }

  /**
   * Divides all the entries of A, B and pi by a divisor.
   */
  void Hmm::divide(int divisor) {
    for (int i = 0; i < numberofstates_; i++) {
      for (int j = 0; j < numberofstates_; j++) {
        A_[i][j] = A_[i][j] / divisor;
      }
    }
    for (int i = 0; i < numberofstates_; i++) {
      for (int j = 0; j < numberofemissions_; j++) {
        B_[i][j] = B_[i][j] / divisor;
      }
    }
    for (int i = 0; i < numberofstates_; i++) {
      pi_[i] = pi_[i] / divisor;
    }
}

  /**
   * Multiplies all the entries of A, B and pi with a factor.
   */
  void Hmm::multiply(int factor) {
    for (int i = 0; i < numberofstates_; i++) {
      for (int j = 0; j < numberofstates_; j++) {
        A_[i][j] = A_[i][j] * factor;
      }
    }
    for (int i = 0; i < numberofstates_; i++) {
      for (int j = 0; j < numberofemissions_; j++) {
        B_[i][j] = B_[i][j] * factor;
      }
    }
    for (int i = 0; i < numberofstates_; i++) {
      pi_[i] = pi_[i] * factor;
    }
  }

  /**
   * Adds all the entries of A, B and pi of a HMM to the corresponding entries
   * of this HMM.
   */
  void Hmm::add(Hmm* hmm){
    for (int i = 0; i < numberofstates_; ++i){
      for (int j = 0; j < numberofstates_; ++j){
        A_[i][j] += hmm->A_[i][j];
      }
      pi_[i] += hmm->pi_[i];
    }

    for (int i = 0; i < numberofstates_; ++i){
      for (int j = 0; j < numberofemissions_; ++j){
        B_[i][j] += hmm->B_[i][j];
      }
    }
  }


  void Hmm::infoAB(){
        std::cerr << "Matrice A : "<<std::endl;
        for(int i=0;i<numberofstates_;++i)
        {
            for(int j=0;j<numberofstates_;++j)
            {
                if(A_[i][j]<0.0001)
                {
                    std::cerr << 0 << " ";
                }
                else
                {
                    std::cerr << A_[i][j] << " ";
                }
            }
            std::cerr << std::endl;
        }/*
        std::cerr << "Matrice B : "<<std::endl;
        for(int i=0;i<numberofstates_;++i)
        {
            for(int j=0;j<numberofemissions_;++j)
            {
                std::cerr << B_[i][j] << " ";
            }
            std::cerr << std::endl;
        }*/
    }




} /* namespace ducks */
