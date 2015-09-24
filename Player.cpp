#include "Player.hpp"
#include <cstdlib>
#include <iostream>

namespace ducks
{

Player::Player()
{
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */

    // This line choose not to shoot
    ///pState.info();
    int NbBirds = (int)pState.getNumBirds();
    if(NbBirds!=0)
    {
        int mvt = pState.getBird(0).getSeqLength();
        if(mvt>20)
        {
            for(int i=0;i<NbBirds;++i)
            {
                Hmm* h = new Hmm(5,9);
                std::cout << "lol";
                ///h.infoAB();
                int l = pState.getBird(i).getSeqLength();
                int *p = new int[l];
                for(int j=0;j<l;++j)
                {
                    int u = toInt(pState.getBird(i).getObservation(j));
                    p[j] = u;
                }
                h->estimateModel(p, l);
                ///h.infoAB();
                double* pi = new double[9];
                for(int j=0;j<5;++j)
                {
                    pi[j] = 0;
                }
                pi[0]=1;
/*
                double* d = h.estimateProbabilityDistributionOfNextEmission(pi);
                std::cerr << "Vecteur :";
                for(int j=0;j<9;++j)
                {
                    std::cerr << d[j];
                }
                std::cerr << std::endl;*/
                ///h.info();
                ///Recupération du vecteur de données


            }
        }
    }


    //This line would predict that bird 0 will move right and shoot at it
    //return Action(0, MOVE_RIGHT);
    return cDontShoot;
}

std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */

    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
    return lGuesses;
}

void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
{
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
     pState.info();
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
}


} /*namespace ducks*/
