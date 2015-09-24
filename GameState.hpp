#ifndef _DUCKS_GAMESTATE_HPP_
#define _DUCKS_GAMESTATE_HPP_

#include "Bird.hpp"
#include <vector>
#include "iostream"

namespace ducks
{

/**
 * Represents a game state
 */
class GameState
{
public:

    void info() const
    {
        if((int)getNumBirds()==0)
        {
            std::cerr << "No bird on the game" << std::endl;
        }
        else
        {
            for(unsigned int i=0;i<getNumBirds();++i)
            {
                std::cerr << " Bird : " << i;
                for(int j=0;j<getBird(i).getSeqLength();++j)
                {
                    std::cerr << toString(getBird(i).getObservation(j)) << " ";
                }
                std::cerr << std::endl;
            }
        }
    }

    std::string toString(EMovement e) const
    {
        switch(e)
        {
            case MOVE_DOWN:
                return "MOVE_DOWN";
                case MOVE_DOWN_LEFT:
                return "MOVE_DOWN_LEFT";
                case MOVE_DOWN_RIGHT:
                return "MOVE_DOWN_RIGHT";
                case MOVE_LEFT:
                return "MOVE_LEFT";
                case MOVE_RIGHT:
                return "MOVE_RIGHT";
                case MOVE_UP:
                return "MOVE_UP";
                case MOVE_UP_LEFT:
                return "MOVE_UP_LEFT";
                case MOVE_UP_RIGHT:
                return "MOVE_UP_RIGHT";
                case MOVE_STOPPED:
                return "MOVE_STOPPED";
                case MOVE_DEAD:
                return " >MOVE_DEAD< ";
                default:
                    return "MVT NON IDENTIFIE !";
        }
    }




    ///returns what round we are currently playing
    int getRound() const
    {
        return mRound;
    }

    ///returns the number of birds
    size_t getNumBirds() const
    {
        return mBirds.size();
    }

    ///returns a reference to the i-th bird
    const Bird &getBird(int i) const
    {
        return mBirds[i];
    }

    ///returns the index of your player among all players
    int whoAmI() const
    {
        return mWhoIAm;
    }

    ///returns the number of players
    int getNumPlayers() const
    {
        return mScores.size();
    }

    ///returns your current score
    int myScore() const
    {
        return mScores[mWhoIAm];
    }

    ///returns the score of the i-th player
    int getScore(int i) const
    {
        return mScores[i];
    }

    ///returns the number of turns elapsed since last time Shoot was called.
    ///this is the amount of new data available for each bird
    int getNumNewTurns() const
    {
        return mNumNewTurns;
    }

    /**
     * The following methods are used by the Client class.
     * Don't use them yourself!
     */
    GameState(int pWhoIAm, int pNumPlayers)
        : mRound(-1)
        , mWhoIAm(pWhoIAm)
        , mNumNewTurns(0)
        , mScores(pNumPlayers, 0)
    {
    }

    void newRound(int pRound, int pNumBirds)
    {
        // Clear the sky and fill it with new birds
        mRound = pRound;
        mBirds.assign(pNumBirds, Bird());
        mNumNewTurns = 0;
    }

    void addMoves(std::vector<EMovement> pMoves)
    {
        for (size_t b = 0; b < mBirds.size(); ++b)
            mBirds[b].addObservation(pMoves[b]);
        mNumNewTurns += 1;
    }

    void setScores(std::vector<int> pScores)
    {
        mScores = pScores;
    }

    void resetNumNewTurns()
    {
        mNumNewTurns = 0;
    }

private:
    int mRound;
    int mWhoIAm;
    int mNumNewTurns;
    std::vector<Bird> mBirds;
    std::vector<int> mScores;
};

} /*namespace ducks*/

#endif
