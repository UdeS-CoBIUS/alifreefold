//============================================================================
// Name        : alifreefold program
// Author      : Jean-Pierre Sehi Glouzon
// Copyright   : GNU/GPL
// Description : alifreefold algorithm in C++, Ansi-style
//============================================================================

#include<string>
#include<stdio.h>
#include <iostream>
#include <iostream>
#include <string>
#include <fstream>
#include <cstring>
#include <vector>
#include <ctype.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <queue>
#include <stdexcept>
#include "RepMotif.h"
#include "RepWeightedMotif.h"
#include "Eigen/Dense"
#include <cstdio>

using namespace std;


//global variable
extern int minNumberOfSubOpt=25;

string filenameTemp;
clock_t t1,t2,t3,t4;

//Input parameters
string i_PathFileInput;//input file
string o_PathDirOutput;//Output directory
string c_PathFileInput;
string f_outputOption;
int computeConservationIndex;
int maxLevelNmotifs;//level of nmotifs

void initParameters(string&, string&, int&, string&,int&,int&);
void setParameters(int, char*[] ,string&,string&,int&, string&,int&,int&);
void computeNwriteMatDissimS2S(const Eigen::ArrayXXf, const string ,const string);
void computeNWriteMatDissimS2NM(const Eigen::ArrayXXf,const Eigen::ArrayXXf, const string ,const string);

void writeRepStructure(const string, const string, const string, const string, const string );
void writeMatALLMotif(const vector<string> ,const Eigen::ArrayXXf ,map<string,float> ,const string, const string );
void writeMatALLMotifWeighted(const vector<string>,const Eigen::ArrayXXf, map<string,float>,const string ,const string );
void writeWeight(const Eigen::RowVectorXf,const string,const string);
void writeDistToCentroid(const Eigen::RowVectorXf,const string,const string);
void writeSubOpt(const vector<string> ,const vector<string> ,const vector<string> ,const string,const string);

void writeMatNmotifsPositionsInSS(const vector<string> , vector< map<string,vector<int>> > , map<string,float>  ,const string , const string );

string help();

int main(int argc, char *argv[])
{
    t1=clock();t2=clock();t3=clock();t4=clock();

    initParameters(i_PathFileInput,o_PathDirOutput,maxLevelNmotifs, f_outputOption,minNumberOfSubOpt,computeConservationIndex);
    setParameters(argc, argv, i_PathFileInput, o_PathDirOutput, maxLevelNmotifs,f_outputOption,minNumberOfSubOpt,computeConservationIndex);

    cout<<endl<<"Running the alifreefold program..."<<endl;

    //Extract n-motifs.
    cout<<endl<<"Compute the n-motifs representation..."<<endl;
    RepNmotif repnMotif(i_PathFileInput,maxLevelNmotifs);
    t1=clock()-t1;
    cout << "Took:" << ((float)t1)/CLOCKS_PER_SEC << " sec."<<endl;

    //Filter and weight n-motifs : n-motifs representation.
    cout<<"Compute the weighted n-motifs and extract representative structure..."<<endl;
    RepWeightedMotif repFilterWeightNmotif(repnMotif.getNmotifsAllStructure(),repnMotif.getNmotifsForEachStructure(),repnMotif.getHeaders(),repnMotif.getSequences(),repnMotif.getStructures(),minNumberOfSubOpt,computeConservationIndex);
    t2=clock()-t2;
    cout << "Took:" << ((float)t2)/CLOCKS_PER_SEC << " sec."<<endl;

    //The super-n-motifs representation
    cout<<"Write output files ..."<<endl;
    for (unsigned int i=0;i<f_outputOption.size();i++){
        switch (f_outputOption[i])
        {
        case 'a':
        cout<<"Write the representative structure in vienna format. (rep.db) ..."<<endl;
        filenameTemp="rep.db";
        writeRepStructure(repFilterWeightNmotif.getHeaderOfRepresentative(),repFilterWeightNmotif.getSequenceOfRepresentative() , repFilterWeightNmotif.getStructureOfRepresentative(),filenameTemp,o_PathDirOutput);
        break;

        case 'b':
        cout<<"Write the raw n-motif representation of SS.(r_nmRepOfSS.csv) ..."<<endl;
        filenameTemp="r_nmRepOfSS.csv";
        writeMatALLMotif(repnMotif.getHeaders(), repFilterWeightNmotif.getMatALLMotif() ,repnMotif.getNmotifsAllStructure(),o_PathDirOutput,filenameTemp);
        break;

        case 'c':
        cout<<"write the weigthed n-motif representation of SS. (w_nmRepOfSS.csv) ..."<<endl;
        filenameTemp="w_nmRepOfSS.csv";
        writeMatALLMotifWeighted(repnMotif.getHeaders(), repFilterWeightNmotif.getMatAllMotifWeigthed(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs(), o_PathDirOutput,filenameTemp);
        break;

        case 'd':
         //write the nucleotide position of n-motifs in each SS
        cout<<"write the nucleotide position of n-motifs in each SS (posOfnmotifs.csv) ..."<<endl;
        filenameTemp="posOfnmotifs.csv";
        writeMatNmotifsPositionsInSS(repnMotif.getHeaders(),repnMotif.getNmotifsForEachStructureWithPosNucOfnmotifs(), repFilterWeightNmotif.getAllStructureFeatureForWeigthOfMotifs() ,o_PathDirOutput, filenameTemp);
        break;

        case 'e':
         //write the vector of weight
        cout<<"write vector of weight (weights.csv) ..."<<endl;
        filenameTemp="weights.csv";
        writeWeight(repFilterWeightNmotif.getAllWeight(),filenameTemp,o_PathDirOutput);
        break;

        case 'f':
         //write the vector of distance of structures to centroid
        cout<<"write vector of distance of structures to centroid (distToCentroid.csv) ..."<<endl;
        filenameTemp="distToCentroid.csv";
        writeDistToCentroid(repFilterWeightNmotif.getDistToCentroid(),filenameTemp,o_PathDirOutput);
        break;

        case 'g':
         //write the vector of distance of structures to centroid
        cout<<"write generated suboptimal structures (subopt.db) ..."<<endl;
        filenameTemp="subopt.db";
        writeSubOpt(repnMotif.getHeaders() ,repnMotif.getSequences(),repnMotif.getStructures() ,o_PathDirOutput,filenameTemp);
        break;

        default:
        break;
        }
    }
    t3=clock()-t3;
    cout << "Took:" << ((float)t3)/CLOCKS_PER_SEC << " sec."<<endl<<endl;


    cout<<"Execution of aliFreeFold program completed"<<endl;
     t4=clock()-t4;
    cout << "Took:" << ((float)t4)/CLOCKS_PER_SEC << " sec."<<endl<<endl;


    return 0;
}

/******************************************************************************/
void initParameters(string& i_PathFileInput,string& o_PathDirOutput,int& maxLevelNmotifs, string& f_outputOption, int& minNumberOfSubOpt, int& computeConservationIndex){
    i_PathFileInput=""; //input file
    o_PathDirOutput=""; //Output directory
    maxLevelNmotifs=1; //n-motifs parameters
    f_outputOption="a";
    minNumberOfSubOpt=25;
    computeConservationIndex=1;
}

/******************************************************************************/
void setParameters(int argc,char* argv[],string& i_PathFileInput,string& o_PathDirOutput,int& maxLevelNmotifs, string& f_outputOption, int& minNumberOfSubOpt, int& computeConservationIndex){

    ifstream viennaFile;
    ifstream viennaFile1;

    string lineHelpFile;
    bool paramRequiredinput=false;
    bool paramRequiredoutput=false;

    struct stat sb;
    for (int i=0;i<argc;i++){
        if (argv[i][0]=='-'){
            switch ( argv[i][1] ) {
            case 'i':
                if (argv[i+1]!=nullptr)
                {
                    i_PathFileInput=string(argv[i+1]);
                    viennaFile.open(i_PathFileInput.c_str());
                      if (!(viennaFile.is_open()))
                      {throw invalid_argument("Unable to open fasta file. Please check the input file path.");}
                    viennaFile.close();
                    paramRequiredinput=true;


                }
                else
                {   throw invalid_argument("Empty value for parameter -i. Please enter the input file path.");}
            break;
            case 'o':
                if (argv[i+1]!=nullptr)
                {
                 o_PathDirOutput=string(argv[i+1]);
                 if (!(stat(o_PathDirOutput.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
                  {
                    cerr<<"The ouput path directory doesn't exits. It will be created."<<endl;
                    #if defined(_WIN32)
                        mkdir(o_PathDirOutput.c_str());
                         #else
                        mkdir(o_PathDirOutput.c_str(), 0700);
                         #endif
                  }
                 paramRequiredoutput=true;
                }
                else
                {   throw invalid_argument("Empty value for parameter -o. Please enter the output directory path.");}
              break;
            case 'm':
                if (argv[i+1]!=nullptr)
                {
                    minNumberOfSubOpt=stoi(argv[i+1]);
                    if (!(minNumberOfSubOpt>=5))
                    {throw invalid_argument("Maximum number of Sub-optimal structures for each sequence (-m). It must be greater or equal to 5. ");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -n. "
                                           "Maximum number of Sub-optimal structures for each sequence (-m). It must be greater or equal to 5. ");};
              break;
            case 'c':
                if (argv[i+1]!=nullptr)
                {
                    computeConservationIndex=stoi(argv[i+1]);
                    if (!(computeConservationIndex==0)&&!(computeConservationIndex==1))
                    {throw invalid_argument("Activate (1) or deactivate (0) computation of the conservation index. It must be 0 or 1. ");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -c. \n"
                                           "It must be 1 or 0 to activate or deactivate the computation of conservation index .");};
              break;
              case 'f':
                if (argv[i+1]!=nullptr)
                {
                    f_outputOption=argv[i+1];
                    string characterAllowed = "abcdefg";
                    if ( f_outputOption.find_first_not_of(characterAllowed)!= string::npos)
                    {throw invalid_argument("Output file option (-f) must be a combination of the following characters : a,b,c,d,e and f.\n"
                                            "   a : the representative structure. (rep.db)\n"
                                            "   b : the raw n-motif representation of SS.(r_nmRepOfSS.csv)\n"
                                            "   c : the weigthed n-motif representation of SS.(w_nmRepOfSS.csv)\n"
                                            "   d : the nucleotide position of n-motifs in each SS (posOfnmotifs.csv)\n"
                                            "   e : the vector of weight (weights.csv)\n"
                                            "   f : the vector of distance of structures to centroid (distToCentroid.csv)\n"
                                            "   g : the generated suboptimal structures  (subopt.csv)\n"
                                            );
                    }
                    else
                    {
                        //sort and remove duplicate characters in f_outputOption
                        sort(f_outputOption.begin(), f_outputOption.end());
                        string::iterator it;
                        it = unique (f_outputOption.begin(), f_outputOption.end());
                        f_outputOption.resize( distance(f_outputOption.begin(),it) );
                    }
                }
                else
                {   throw invalid_argument("Empty value for parameter -f. "
                                           "Output file option (-f) must be a combination of the following characters : a,b,c,d,e and f.\n"
                                            "   a : the representative structure. (rep.db)\n"
                                            "   b : the raw n-motif representation of SS.(r_nmRepOfSS.csv)\n"
                                            "   c : the weigthed n-motif representation of SS.(w_nmRepOfSS.csv)\n"
                                            "   d : the nucleotide position of n-motifs in each SS (posOfnmotifs.csv)\n"
                                            "   e : the vector of weight (weights.csv)\n"
                                            "   f : the vector of distance of structures to centroid (distToCentroid.csv)\n"
                                            "   g : the generated suboptimal structures  (subopt.csv)\n"
                                           );
                }
              break;
            case 'n':
                if (argv[i+1]!=nullptr)
                {
                    maxLevelNmotifs=stoi(argv[i+1]);
                    if (!(maxLevelNmotifs==0)&&!(maxLevelNmotifs==1)&& !(maxLevelNmotifs==2))
                    {throw invalid_argument("Maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2. "
                                            "For instance: when it is set to 0, 0-motifs will be extracted. "
                                            "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");}
                }
                else
                {   throw invalid_argument("Empty value for parameter -n. "
                                           "The maximum level of n-motifs parameter (-n) is used to extract the n-motifs. It must be 0, 1 or 2."
                                           "For instance: when it is set to 0, 0-motifs will be extracted. "
                                           "when it is set to 1, 0-motifs and 1-motifs will be extracted. etc.");};
              break;
            case 'h':
                  cout<< help(); exit(EXIT_SUCCESS);
              break;
            default:
                cerr<<"No parameters found or wrong parameters. Please check parameter list in the help (-h)."<<endl;
                exit(EXIT_FAILURE);
            break;
            }
            if ((argv[i][1]!='h')&&(argv[i][1]!='i')&&(argv[i][1]!='o')&&(argv[i][1]!='f')&&(argv[i][1]!='n')&&(argv[i][1]!='c')&&(argv[i][1]!='m'))
            {cerr<<"Wrong parameters. Please check parameter list in the help (-h)."; exit(EXIT_FAILURE);}
        }
    }
    if((paramRequiredinput==false || paramRequiredoutput==false))
   // if((paramRequiredinput==false|| paramRequiredoutput==false))
        {cerr<<"The program requires the input fasta file (-i),  and the output folder (-o)."<<endl; exit(EXIT_FAILURE);}
}

/*********************************************************************************************/
string help()
{
    string help;
    help=
"\n"
"  ### The aliFreeFold program -- predict RNA secondary structure. ###\n"
"\n"
"  Usage:\n"
"\n"
"  /path_to_alifreefold_program/alifreefold -i path_to_fileInDb\n"
"  -o path_to_folderOfResults\n"
"\n"
"  The alifreefold program takes as input a file of sequences\n"
"  in fasta format (-i parameter):\n"
"\n"
"  >RNA1\n"
"  GCCCCGCUGAUGAGGUCAGGGAAAACCGAAAGUGUCGACUCUACGGGGC\n"
"\n"
"  It ouputs the representative structure by default\n"
"  (-o parameter). For further options:\n"
"  '/path_to_alifreefold_program/alifreefold -h'\n"
"\n"
"  ### parameters ###\n"
"\n"
"  -h"
"\n"
"    Help regarding the parameters of aliFreeFold program.\n"
"\n"
"  -i [input vienna file]\n"
"\n"
"    The Input file is in fasta format (required).\n"
"    >Y08502.1-137669_137741\n"
"    ACCUACUUGACUCAGCGGUUAGAGUAUCGCUUUCAUACGGCGAGAGUCAUUGGUUCAAAUCCAAUAGUAGGUA\n"
"    >AF070678.1-91_163\n"
"    GGGGCCUUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCAGCGGUUCGAUCCCGCUAGGCUCCA\n"
"    >AJ271079.2-114727_114656\n"
"    UCCUCAGUAGCUCAGUGGUAGAGCGGUCGGCUGUUAACCGAUUGGUCGUAGGUUCGAAUCCUACUUGGGGAG\n"
"    >X61698.1-1470_1542\n"
"    ACCUACUUAACUCAGUGGUUAGAGUACUGCUUUCAUACGGCGGGAGGCAUUGGUUCAAAUCCAAUAGUAGGUA\n"
"    >V00654.1-12038_12108\n"
"    ACUUUUAAAGGAUAGUAGUUUAUCCGUUGGUCUUAGGAACCAAAAAAUUGGUGCAACUCCAAAUAAAAGUA\n"
"\n"
"  -o [output folder]\n"
"\n"
"  -c [1|0]\n"
"\n"
"    Activation (1) or deactivation (0) of computation of \n"
"    conservation index.\n"
"\n"
"    The result folder (required).\n"
"\n"
"  -n [0|1|2]\n"
"\n"
"    The maximum level of n-motifs used to\n"
"    extract the n-motifs. It must be 0, 1 or 2. For instance: when\n"
"    it is set to 0, 0-motifs will be extracted. While it is set\n"
"    to 1, 0-motifs and 1-motifs will be extracted. etc.\n"
"    (default : -n 1)\n"
"\n"
"  -m [>=5]\n"
"\n"
"    The maximum number of sub-optimal structures for each sequence.\n"
"    It must be greater or equal to 5.\n"
"    (default : -m 25)\n"
"\n"
"  -f [abcdefg]\n"
"\n"
"    The output file option.\n"
"    (default: -f a)\n"
"   a : the representative structure. (rep.db)\n"
"   b : the raw n-motif representation of SS.(r_nmRepOfSS.csv)\n"
"   c : the weigthed n-motif representation of SS.(w_nmRepOfSS.csv)\n"
"   d : the nucleotide position of n-motifs in each SS (posOfnmotifs.csv)\n"
"   e : the vector of weight (weights.csv)\n"
"   f : the vector of distance of structures to centroid (distToCentroid.csv)\n"
"   g : the generated suboptimal structures  (subopt.csv)\n"
"\n"
"  ### Licence ###\n"
"  The aliFreeFold program is released under the terms of the GNU GPL licence.\n\n";


return help;
}

/******************************************************************************/
void writeRepStructure(const string headerOfRepresentative, const string sequenceOfRepresentative, const string structureOfRepresentative, const string filename, const string o_PathDirOutput ){

    string tempfilename;
    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    fprintf(outputFile, ">%s\n", headerOfRepresentative.c_str());
    fprintf(outputFile, "%s\n", sequenceOfRepresentative.c_str());
    fprintf(outputFile, "%s\n", structureOfRepresentative.c_str());

    fclose(outputFile);
}

/******************************************************************************/
void writeMatALLMotif(const vector<string> headers , const Eigen::ArrayXXf matS_nmotifs, map<string,float> allStructureFeature  ,const string o_PathDirOutput, const string filename)
{
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write labels
    map<string,float>::iterator it1;
    fprintf(outputFile, " ");

    for(it1=allStructureFeature.begin();it1!=allStructureFeature.end();++it1)
    {fprintf(outputFile, ",%s", it1->first.c_str());}
    fprintf(outputFile, "\n");

    //write data
    for (unsigned int i=0;i<matS_nmotifs.rows();i++)
    {
     fprintf(outputFile, "%s,", headers[i].c_str());
        for (unsigned int j=0;j<matS_nmotifs.cols();j++)
        {fprintf(outputFile, "%f,", matS_nmotifs(i,j));}
        fprintf(outputFile, "\n");
    }

    fclose(outputFile);
}

/******************************************************************************/
void writeMatALLMotifWeighted(const vector<string> headers, const Eigen::ArrayXXf matALLMotifWeighted, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){

    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write n-motifs labels
    map<string,float>::iterator it;
    fprintf(outputFile, " ");

    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {fprintf(outputFile, ",%s", it->first.c_str());}
    fprintf(outputFile, "\n");

    //write data
    for (unsigned int i=0;i<matALLMotifWeighted.rows();i++)
    {
     fprintf(outputFile, "%s,", headers[i].c_str());
        for (unsigned int j=0;j<matALLMotifWeighted.cols();j++)
        {fprintf(outputFile, "%f,", matALLMotifWeighted(i,j));}
        fprintf(outputFile, "\n");
    }

    fclose(outputFile);
}

//Write all Weight
/******************************************************************************/
void writeWeight(const Eigen::RowVectorXf allWeight, const string filename,const string o_PathDirOutput){
    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write data
    for (unsigned int i=0;i<allWeight.size();i++)
    {fprintf(outputFile, "%f,", allWeight[i]);}
    fprintf(outputFile, "\n");

    fclose(outputFile);
}

//Write distCentroid
void writeDistToCentroid(const Eigen::RowVectorXf distToCentroid, const string filename,const string o_PathDirOutput){

    string tempfilename;

    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write data
    for (unsigned int i=0;i<distToCentroid.size();i++)
    {fprintf(outputFile, "%f,", distToCentroid[i]);}
    fprintf(outputFile, "\n");

    fclose(outputFile);

}

/******************************************************************************/
void writeMatNmotifsPositionsInSS(const vector<string> headers, vector< map<string,vector<int>> > FeatureForEachStructureWithPosNuc, map<string,float> allStructureFeatureForWeigthOfMotifs ,const string o_PathDirOutput, const string filename){

    string tempfilename;
    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write n-motifs labels
    map<string,float>::iterator it;
    fprintf(outputFile, " ");
    for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
    {fprintf(outputFile, ",%s",it->first.c_str());}
    fprintf(outputFile, "\n");

    //write data
    map<string,vector<int>>::iterator it2;
    for(unsigned int i=0;i<headers.size();i++)
    {
        fprintf(outputFile, "%s,",headers[i].c_str());
        for(it=allStructureFeatureForWeigthOfMotifs.begin();it!=allStructureFeatureForWeigthOfMotifs.end();++it)
        {
            it2=FeatureForEachStructureWithPosNuc[i].find(it->first);
            if (it2!=FeatureForEachStructureWithPosNuc[i].end())
            {
                for(unsigned int j=0;j<it2->second.size();j++)
                {fprintf(outputFile, "%i|",it2->second[j]+1);}
            }
            else
            {fprintf(outputFile, "x");}
            fprintf(outputFile, ",");
        }
        fprintf(outputFile, "\n");
    }
    fclose(outputFile);
}

/******************************************************************************/
void writeSubOpt(const vector<string> headers,const vector<string> sequences,const vector<string> structures,const string o_PathDirOutput, const string filename){

   string tempfilename;
    tempfilename=o_PathDirOutput+"/"+filename;
    FILE* outputFile = fopen(tempfilename.c_str(), "wb");

    //write data
    map<string,vector<int>>::iterator it2;
    for(unsigned int i=0;i<headers.size();i++)
    {
        fprintf(outputFile, ">%s\n",headers[i].c_str());
        fprintf(outputFile, "%s\n",sequences[i].c_str());
        fprintf(outputFile, "%s\n",structures[i].c_str());
    }
    fclose(outputFile);
}

