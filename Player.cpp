#include "Player.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>

namespace ducks
{

Player::Player()
{
    NStates = 5;
}

Action Player::shoot(const GameState &pState, const Deadline &pDue)
{
    /*
     * Here you should write your clever algorithms to get the best action.
     * This skeleton never shoots.
     */
    int NEmissions = 9;
    int NbBirds = (int)pState.getNumBirds();
    double *tabMax = new double[NbBirds];
    double *tabMaxCoup = new double[NbBirds];
    for(int i=0;i<NbBirds;++i)
    {
        tabMax[i]=0;
        tabMaxCoup[i]=0;
    }
    if(NbBirds!=0)
    {
        int mvt = pState.getBird(0).getSeqLength();
        std::vector<Hmm*> listHMM;
        if(mvt>20)
        {
            for(int i=0;i<NbBirds;++i)
            {
                if(pState.getBird(i).isAlive())
                {
                    //std::cerr << "---------------------------------------Oiseau "<<i<<"----------------------------"<<std::endl;
                    Hmm *h = new Hmm(NStates, NEmissions);
                    int l = pState.getBird(i).getSeqLength();
                    int *p = new int[l];
                    int *tuit = new int[3];

                    for(int j=0;j<l;++j)
                    {
                        int u = toInt(pState.getBird(i).getObservation(j));
                        //std::cerr<<u;
                        p[j] = u;
                    }
                    //std::cerr<<std::endl;
                    /*tuit[0] = p[l-4];
                    tuit[1] = p[l-3];*/
                    tuit[0] = p[l-2];
                    tuit[1] = p[l-1];
                    tuit[2] = 0;
                    h->estimateModel(p, l);
                    //h->infoAB();
                    double maxi = 0.0;
                    int maxiInd = 0;
                    for(int j=0;j<9;++j)
                    {
                        tuit[2] = j;
                        double test =h->estimateProbabilityOfEmissionSequence(tuit, 3);
                        if(test>maxi)
                        {
                            maxiInd = j;
                            maxi = test;
                        }
                    }
                    tabMax[i] = maxi;
                    tabMaxCoup[i] = maxiInd;

                    listHMM.push_back(h);
                }
            }
            //std::cout << "loiujl";


            ///Récupération du maximum
            int maximaxiInd = 0;
            int indBird = 0;
            double maximaxi = 0;
            for(int i=0;i<NbBirds;++i)
            {
                if(maximaxi < tabMax[i])
                {
                    maximaxi = tabMax[i];
                    maximaxiInd = tabMaxCoup[i];
                    indBird = i;
                }
            }
            //std::cerr << maximaxi << std::endl;
            if(maximaxi>0.7)
            {
                std::cerr << "TIR :: " << maximaxi << " : " << indBird << " : " << maximaxiInd<< std::endl;
                return Action(indBird, toEMove(maximaxiInd));
            }

            ///TEst des interactions entre les différentes hmm
            /*
            listHMM[0]->infoAB();
            listHMM[1]->infoAB();
            listHMM[2]->infoAB();
            listHMM[3]->infoAB();
            listHMM[4]->infoAB();
            listHMM[5]->infoAB();*/
            /*
            double maxiCor = 0.0;
            int Ind1 = 0;
            int Ind2 = 0;
            double** matCorelation =  new double*[listHMM.size()];
            std::cerr << "Matrice de corrélation des matrices A : " << std::endl;
            for(unsigned int i=0;i<listHMM.size();++i)
            {
                matCorelation[i] = new double[listHMM.size()];
                for(unsigned int j=0;j<listHMM.size();++j)
                {
                    matCorelation[i][j] = prodScalaire(listHMM[i]->getA(), listHMM[j]->getA(), NStates, NStates);
                    if(maxiCor < matCorelation[i][j]&&i!=j)
                    {
                        maxiCor = matCorelation[i][j];
                        Ind1 = i;
                        Ind2 = j;
                    }
                    std::cerr << matCorelation[i][j] << " ";
                }
                std::cerr << std::endl;
            }
            std::cerr << maxiCor << " " << Ind1 << " " << Ind2 << std::endl;
            listHMM[Ind1]->infoAB();
            listHMM[Ind2]->infoAB();
            for(int i=0;i<pState.getBird(Ind1).getSeqLength();++i)
            {
                std::cerr << toInt(pState.getBird(Ind1).getObservation(i)) << " ";
            }
            std::cerr << std::endl;
            for(int i=0;i<pState.getBird(Ind2).getSeqLength();++i)
            {
                std::cerr << toInt(pState.getBird(Ind2).getObservation(i)) << " ";
            }

            ///Mettre les résultats dans une matrice NBirds²
            std::cout << "Fin";
            ///Fin test décommenter ligne précédente
            */
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
    std::cerr << "HIT BIRD!!!" << std::endl;
}

void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
{
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */

}



///----------------------Fonctions utiles pour la manipulation de matrice sous forme de tableau (int**) -------------------------

double Player::distanceMatriceN2(double** a, double** b, int taille)
{
    double result = 0;
    for(int i=0;i<taille;++i)
    {
        for(int j=0;j<taille;++j)
        {
            result += (a[i][j]-b[i][j])*(a[i][j]-b[i][j]);
        }
    }
    return sqrt(result);
}

double Player::distanceMatriceN1(double** a, double** b, int taille)
{
    double result = 0;
    for(int i=0;i<taille;++i)
    {
        for(int j=0;j<taille;++j)
        {
            result += abs(a[i][j]-b[i][j]);
        }
    }
    return result;
}

double** Player::produit(double** a, double** b, int taille)
{
    double** result = new double*[taille];
    for(int i=0;i<taille;++i)
    {
        result[i] = new double[taille];
        for(int j=0;j<taille;++j)
        {
            result[i][j] = 0;
        }
    }
    for(int i=0;i<taille;++i)
    {
        for(int j=0;j<taille;++j)
        {
            for(int k=0;k<taille;++k)
            {
                result[i][j] += a[i][k]*b[k][j];
            }
        }
    }
    return result;
}

double** Player::transpose(double **a, int taille)
{
    double** result = new double*[taille];
    for(int i=0;i<taille;++i)
    {
        result[i] = new double[taille];
        for(int j=0;j<taille;++j)
        {
            result[i][j] = a[j][i];
        }
    }
    return result;
}

double** Player::creationMatricePermutation(int* tabPermutation, int taille)
{
    double** result = new double*[taille];
    for(int i=0;i<taille;++i)
    {
        result[i] = new double[taille];
        for(int j=0;j<taille;++j)
        {
            result[i][j] = 0.0;
        }
    }
    for(int i=0;i<taille;++i)
    {
        result[i][tabPermutation[i]] = 1.0;
    }
    return result;
}

void Player::afficherMatrice(double** a, int taille)
{
    std::cerr << std::endl;
    for(int i=0;i<taille;++i)
    {
        for(int j=0;j<taille;++j)
        {
            std::cerr << a[i][j] << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

void Player::swapTab(int* a, int i, int j)
{
    int temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}

void Player::permutation(int* tab, int ind, int taille, double** mat1, double** mat2, double &dist, double** closerMatrix)
{
    if(ind==taille-1)
    {
        double** p1 = creationMatricePermutation(tab, taille);
        double** p2 = transpose(p1, taille);
        double** product = produit(p2, produit(mat1, p1, taille), taille);
        double distanceProvisoireN2 = distanceMatriceN2(product, mat2, taille);
        if(dist>distanceProvisoireN2)
        {
            dist = distanceProvisoireN2;
            for(int i=0;i<taille;++i)
            {
                closerMatrix[i] = product[i];
            }
        }
        distanceProvisoireN2 += 0.2;

    }
    else
    {
        for(int i=ind;i<taille;++i)
        {
            swapTab(tab, ind, i);
            permutation(tab, ind+1, taille, mat1, mat2, dist, closerMatrix);
            swapTab(tab, ind, i);
        }
    }
}


} /*namespace ducks*/
