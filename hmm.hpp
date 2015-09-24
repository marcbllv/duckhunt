#ifndef _DUCKS_HMM_HPP_
#define _DUCKS_HMM_HPP_

namespace ducks {

/**
 * Represents a HMM.
 *
 * Much of the terminology is borrowed from a very pedagogic paper by Rabiner 
 * and Juang: http://www.cs.ubc.ca/~murphyk/Bayes/rabiner.pdf
 *
 * We recommend that you read the article as an introduction to HMMs.
 */
class Hmm {
public:
	Hmm(int numberofstates, int numberofemissions);
	Hmm(int numberofstates, int numberofemissions, double** A, double** B, double* pi);
	
	double* estimateProbabilityDistributionOfNextEmission(
      double* currentStateProbabilityDistribution); //HMM1
    double estimateProbabilityOfEmissionSequence(int* O, int numberofobservations); //HMM2
    int* estimateStateSequence(int* O, int numberofobservations); //HMM3
    void estimateModel(int* O, int numberofobservations); //HMM4
    
    void copyHmm(Hmm* hmm);
    void divide(int divisor);
    void multiply(int factor);
    void add(Hmm* hmm);

	//Add your own functions here...

private:

	//Nothing here yet...
	
public:

    const int kMaxIters = 30; // Max iterations when estimating a new model.
  
	int numberofstates_; // The number of states in the HMM.
	int numberofemissions_; // The number of emissions in the HMM.

	double** A_; // The transition matrix of the HMM.
	double** B_; // The emission matrix of the HMM.
	double* pi_; // The initial state distribution of the HMM.
	
	//Add your own class variables and constants here...
	
};

} /* namespace ducks */

#endif
