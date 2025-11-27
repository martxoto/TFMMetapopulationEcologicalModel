//First experiment

#include "model.h"

int main(void)
{
    double h, t;
    ofstream fichp, fichv;
    
    cout << "Introduce the step: " << endl;
    cin >> h;
    
    //Gamma
    
    map<string, int> plantIndex;
    map<string, int> insectIndex;
    
    int plantCount = 0, insectCount = 0, numPatch = 0;    
    
    vector<vector<vector<double>>> gamma;
    
    loadGamma("interactions_Penhale_Sands_patches.txt", plantIndex,insectIndex,plantCount,insectCount,numPatch,gamma);
    
    cout << "Number of plants: " << plantCount << endl;
    cout << "Number of insects: " << insectCount << endl;
    
    //Initialization
    
    vector <vector<double>> p(plantCount,vector<double>(numPatch));
    vector <vector<double>> v(insectCount,vector<double>(numPatch));
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            if (plantExistsInPatch(i, site, insectCount, gamma)) 
            {
                p[i][site] = 100.0;
            } 
            else 
            {
                p[i][site] = 0.0; 
            } 
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            if (insectExistsInPatch(i, site, plantCount, gamma)) 
            {
                v[i][site] = 500.0;
            } 
            else 
            {
                v[i][site] = 0.0;
            }
        }
    }   
    
    //Run
    fichp.open("evolutionp.txt");
    fichv.open("evolutionv.txt");

    t = 0.0;
    
    findSteadyState(t, p, v, fichp, fichv, plantCount, insectCount, numPatch, gamma, h);
    
    runExtinctionExperiment(p,v,gamma,h,plantCount,insectCount,numPatch);
    
    fichp.close();
    fichv.close();
    
    return 0;
}
