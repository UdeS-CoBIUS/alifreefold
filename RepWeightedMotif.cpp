/*
 * RepWeightedMotif.cpp
 *
 *  Created on: Jul 28, 2014
 *      Author: Séhi
 */

#include "RepWeightedMotif.h"

/*********************************************************************************/

RepWeightedMotif::RepWeightedMotif() {
	// TODO Auto-generated constructor stub
}

/*********************************************************************************/

RepWeightedMotif::RepWeightedMotif(map<string,float> allStructureFeature, vector< map<string,float> > nmotifsForEachStructure,vector<string> headers,vector<string> sequences,vector<string> structures,int minNumberOfSubOpt,int computeConservationIndex)
{
	// TODO Auto-generated constructor stub
    this->allStructureFeatureForWeigthOfMotifs=allStructureFeature;
    this->nmotifsForEachStructure=nmotifsForEachStructure;
    //this->sequenceOfRepresentative="";
    //this->structureOfRepresentative="";
    this->headers=headers;
    this->sequences=sequences;
    this->structures=structures;

    unsigned int nrows=nmotifsForEachStructure.size();
    unsigned int k=0;
    Eigen::RowVectorXi suboptidx(nrows);
    int nbUniqueSeq=nrows/minNumberOfSubOpt;

    for (int i=0;i<nbUniqueSeq;i++)
    {
        for (int j=0;j<minNumberOfSubOpt;j++)
        {
            suboptidx(k)=i;
            k++;
        }
    }


    this->subOptIdx=suboptidx;

    this->matALLMotifWeigthed=filterMotifs(suboptidx,computeConservationIndex);

    RepWeightedMotif::computeRepresentative(matALLMotifWeigthed,headers,sequences,structures);

   // RepWeightedMotif::computeRepresentative(this->matALLMotifWeigthed,this->headers,this->sequences,this->structures);

}


/*********************************************************************************/

RepWeightedMotif::~RepWeightedMotif() {
	// TODO Auto-generated destructor stub
}

/*******************************************************************************************************/
Eigen::ArrayXXf RepWeightedMotif::filterMotifs(Eigen::RowVectorXi suboptidx,int computeConservationIndex){


    Eigen::ArrayXXf matS_nmotifs = Eigen::ArrayXXf::Zero(nmotifsForEachStructure.size(), allStructureFeatureForWeigthOfMotifs.size());
    std::map<string, float>::iterator it;
    int ind = 0;

    for (unsigned int i = 0; i < nmotifsForEachStructure.size(); i++)
    {
        for (std::map<string, float>::iterator it2 =nmotifsForEachStructure[i].begin();it2 != nmotifsForEachStructure[i].end(); ++it2)
        {
            it = allStructureFeatureForWeigthOfMotifs.find(it2->first);
            ind = distance(allStructureFeatureForWeigthOfMotifs.begin(), it);
            matS_nmotifs(i, ind) = it2->second;
        }
    }

    this->matALLMotif=matS_nmotifs;

    Eigen::ArrayXXf matS_weightdNmotifs = Eigen::ArrayXXf::Zero(nmotifsForEachStructure.size(), allStructureFeatureForWeigthOfMotifs.size());
    Eigen::RowVectorXf allWeight(matS_nmotifs.cols());

    if (computeConservationIndex==1)
    {

        vector<int> countNmotifVect;
        vector<int> uniquecountNmotifVect;

        Eigen::RowVectorXf vectFreq;
        Eigen::RowVectorXf vectFreqLog;
        Eigen::RowVectorXf currentCountCol;

        float currentSum;
        float entCount;
        float entLabel;
        float entY;
        float entXY;

        float currentWeight;
        Eigen::RowVectorXf allEntCount(matS_nmotifs.cols());
        Eigen::RowVectorXf allEntLabel(matS_nmotifs.cols());

        vector<int> currentLabelY;
        vector<int> currentUniqueLabelY;

        vector<int> currentLabelXY;
        vector<int> currentUniqueLabelXY;

        vector<int> currentLabelCol;
        vector<int> currentUniqueLabelCol;

        vector<int> currentCol;
        vector<int> currentUniqueCol;

        int tempTest;
        for (unsigned int j = 0; j < matS_nmotifs.cols(); j++)
        {
                for (unsigned int i = 0; i < matS_nmotifs.matrix().col(j).size(); i++)
                {
                    if (matS_nmotifs(i,j)>0)
                    {currentLabelCol.push_back(suboptidx(i)); }
                    currentCol.push_back(int(matS_nmotifs(i,j)));
                    currentLabelY.push_back(suboptidx(i));

                    tempTest= stoi(std::to_string(int(matS_nmotifs(i,j)))+std::to_string(suboptidx(i)));
                    currentLabelXY.push_back(tempTest);


                }

                //Entropy Count
                 std::sort (currentCol.begin(), currentCol.end());
                currentUniqueCol=currentCol;
                std::vector<int>::iterator it;
                it = std::unique (currentUniqueCol.begin(), currentUniqueCol.end());
                currentUniqueCol.resize( std::distance(currentUniqueCol.begin(),it));

                entCount=0;
                for (unsigned int i = 0; i < currentUniqueCol.size(); i++)
                {
                    currentSum=0;
                    for (unsigned int k = 0; k < currentCol.size(); k++)
                        {
                            if(currentCol[k]==currentUniqueCol[i])
                            {currentSum++;}
                        }
                    entCount=entCount+(currentSum/currentCol.size())*log(currentSum/currentCol.size());
                }
                entCount=entCount*-1;
                currentCol.clear();

                currentWeight=1/exp(entCount);//inverse diversity-->nmfold1

                allWeight(j)=currentWeight;
                matS_weightdNmotifs.col(j)=(matS_nmotifs.col(j))*currentWeight;

        }

    }

    else
    {

        Eigen::RowVectorXf allWeight=Eigen::RowVectorXf::Ones(matS_nmotifs.cols());
        matS_weightdNmotifs=matS_nmotifs;
    }

    this->allWeight=allWeight;

   return matS_weightdNmotifs;
}

/*******************************************************************************************************/
void RepWeightedMotif::computeRepresentative(Eigen::ArrayXXf matALLMotifWeigthed,vector<string> headers,vector<string> sequences,vector<string> structures)
{

    Eigen::VectorXf allNorm;
    int n=matALLMotifWeigthed.rows();

    //Compute Centroid
    Eigen::RowVectorXf centroid;
    centroid=matALLMotifWeigthed.colwise().mean();

	//compute distance between centroids and structures
    Eigen::RowVectorXf distToCentroid(n);

	for ( int i=0; i<n;i++)
	{
	    distToCentroid(i) = (centroid.array() - matALLMotifWeigthed.row(i).array()).matrix().lpNorm<2>();
    }

    //Get structure close to the centroid
    Eigen::ArrayXXf::Index minRow, minCol;
    distToCentroid.array().minCoeff(&minRow, &minCol);

    this->headerOfRepresentative=headers[minCol];
    this->sequenceOfRepresentative=sequences[minCol];
    this->structureOfRepresentative=structures[minCol];

    this->distToCentroid=distToCentroid;

}
