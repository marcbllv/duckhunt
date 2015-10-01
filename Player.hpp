#ifndef _DUCKS_PLAYER_HPP_
#define _DUCKS_PLAYER_HPP_

#include "Deadline.hpp"
#include "GameState.hpp"
#include "Action.hpp"
#include "hmm.hpp"
#include <vector>

namespace ducks
{

class Player
{
private:
    int NStates;

public:
    /**
     * Constructor
     * There is no data in the beginning, so not much should be done here.
     */
    Player();

    int toInt(EMovement e) const
    {
        switch(e)
        {
            case MOVE_DOWN:
                return 0;
            case MOVE_DOWN_LEFT:
                return 1;
            case MOVE_DOWN_RIGHT:
                return 2;
            case MOVE_LEFT:
                return 3;
            case MOVE_RIGHT:
                return 4;
            case MOVE_UP:
                return 5;
            case MOVE_UP_LEFT:
                return 6;
            case MOVE_UP_RIGHT:
                return 7;
            case MOVE_STOPPED:
                return 8;
            case MOVE_DEAD:
                return 9;
                default:
                    return -1;
        }
    }

    EMovement toEMove(int e) const
    {
        switch(e)
        {
            case 0:
                return MOVE_DOWN;
            case 1:
                return MOVE_DOWN_LEFT;
            case 2:
                return MOVE_DOWN_RIGHT;
            case 3:
                return MOVE_LEFT;
            case 4:
                return MOVE_RIGHT;
            case 5:
                return MOVE_UP;
            case 6:
                return MOVE_UP_LEFT;
            case 7:
                return MOVE_UP_RIGHT;
            case 8:
                return MOVE_STOPPED;
            default:
                return MOVE_DEAD;
        }
    }


    ///Fonctions utiles pour la manipulation des matrices (int**)
    void swapTab(int*, int, int);
    void permutation(int*, int, int, double**, double**, double&, double**);
    void afficherMatrice(double**, int);
    double** creationMatricePermutation(int*, int);
    double** transpose(double **, int);
    double** produit(double**, double**, int);
    double distanceMatriceN1(double**, double**, int);
    double distanceMatriceN2(double**, double**, int);

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each birds contains all past actions.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    Action shoot(const GameState &pState, const Deadline &pDue);

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
    std::vector<ESpecies> guess(const GameState &pState, const Deadline &pDue);

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    void hit(const GameState &pState, int pBird, const Deadline &pDue);

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    void reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue);
};

} /*namespace ducks*/

#endif
