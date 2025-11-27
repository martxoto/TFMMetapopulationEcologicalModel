//This program simulates a metapopulation model of pollinators-plants ecosystem

#include "model.h"

bool plantExistsInPatch(int plantID, int site, int insectCount, const vector<vector<vector<double>>>& gamma)
{
    for (int j = 0; j < insectCount; ++j) {
        if (gamma[site][plantID][j] > 0.0) {
            return true;
        }
    }
    return false;
}

bool insectExistsInPatch(int insectID, int site, int plantCount, const vector<vector<vector<double>>>& gamma)
{
    for (int i = 0; i < plantCount; ++i) {
        if (gamma[site][i][insectID] > 0.0) {
            return true;
        }
    }
    return false; 
}

void evaluaFp(double p, const vector<vector<double>>& v, double &fp, const vector<vector<vector<double>>>& gamma, int pindex, int site, int insectCount)
{
    fp = 0.0;
    double sum = 0.0;
    for( int j=0 ; j<insectCount ; j++)
        sum += (gamma[site][pindex][j]*v[j][site]) / (1.0 + (ha*gamma[site][pindex][j]*v[j][site])); 
    fp = (r*p*(1.0-(p/Kp))) + (p*sum); 
    return;
}

void evaluaFv(const vector<vector<double>>& p, const vector<vector<double>>& v, double &fv, const vector<vector<vector<double>>>& gamma, int vindex, int site, int plantCount, int numPatch)
{
    fv = 0.0;
    double sum = 0.0;
    double sumD = 0.0;
    for( int i=0 ; i<plantCount ; i++)
        sum += (gamma[site][i][vindex]*p[i][site]*v[vindex][site]) / (1.0 + (ha*gamma[site][i][vindex]*p[i][site]));
    for( int i=0 ; i<numPatch ; i++)
        if ( i!=site)
            sumD += (v[vindex][i]-v[vindex][site]); 
    fv = (-1.0*d*v[vindex][site]*(1.0+(v[vindex][site]/Kv))) + sum + (D*sumD); 
    return;
}

void rungekutta(vector<vector<double>>& p, vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch)
{   
    double sumav, sumap;
    
    vector<vector<double>> k1v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k2v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k3v(insectCount, vector<double>(numPatch));
    vector<vector<double>> k4v(insectCount, vector<double>(numPatch));
    
    vector<vector<double>> k1p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k2p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k3p(plantCount, vector<double>(numPatch));
    vector<vector<double>> k4p(plantCount, vector<double>(numPatch)); 
    
    vector<vector<double>> fv(insectCount, vector<double>(numPatch));
    vector<vector<double>> fp(plantCount, vector<double>(numPatch));
    
    vector<vector<double>> auxp(plantCount, vector<double>(numPatch));
    vector<vector<double>> auxv(insectCount, vector<double>(numPatch));
    
    //k1p k1v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k1p[i][site] = 0.0;
            evaluaFp(p[i][site],v,fp[i][site],gamma,i,site,insectCount);
            k1p[i][site] = h*fp[i][site];
        }
    }
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k1v[i][site] = 0.0;
            evaluaFv(p,v,fv[i][site],gamma,i,site,plantCount,numPatch);
            k1v[i][site] = h*fv[i][site];
        }
    }
    //k2p k2v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k1p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = 0.0; 
            auxv[i][site] = v[i][site]+(0.5*k1v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k2p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k2p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k2v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch);
            k2v[i][site] = h*fv[i][site];
        }
    }
    //k3p k3v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site]+(0.5*k2p[i][site]);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = 0.0;
            auxv[i][site] = v[i][site]+(0.5*k2v[i][site]);
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k3p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k3p[i][site] = h*fp[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k3v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch);
            k3v[i][site] = h*fv[i][site];
        }
    }
    //k4p k4v
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxp[i][site] = 0.0;
            auxp[i][site] = p[i][site] + k3p[i][site];
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            auxv[i][site] = v[i][site] + k3v[i][site];
        }
    }
    
    for(int i=0 ; i<plantCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k4p[i][site] = 0.0;
            evaluaFp(auxp[i][site],auxv,fp[i][site],gamma,i,site,insectCount);
            k4p[i][site] = h*fp[i][site];
            sumap=0.0;
            sumap=k1p[i][site]+(2*k2p[i][site])+(2*k3p[i][site])+k4p[i][site];
            p[i][site]+=(sumap/6.0);
        }
    }
    
    for(int i=0 ; i<insectCount ; i++)
    {
        for(int site=0 ; site<numPatch ; site++)
        {
            k4v[i][site] = 0.0;
            evaluaFv(auxp,auxv,fv[i][site],gamma,i,site,plantCount,numPatch);
            k4v[i][site] = h*fv[i][site];
            sumav=0.0;
            sumav=k1v[i][site]+(2*k2v[i][site])+(2*k3v[i][site])+k4v[i][site];
            v[i][site]+=(sumav/6.0);
        }
    }
    return;
}

void findSteadyState(double t, vector<vector<double>>& p, vector<vector<double>>& v, ofstream& fichp, ofstream& fichv, int plantCount, int insectCount, int numPatch,  const vector<vector<vector<double>>>& gamma, double h)
{
    //Stationary state detection
    double max_delta = 1.0;
    const double TOLERANCE = 1e-9;
    int iter_count = 0;
    int max_iter = 100000;
    
    vector<vector<double>> p_prev(plantCount, vector<double>(numPatch));
    vector<vector<double>> v_prev(insectCount, vector<double>(numPatch));
    
    while(iter_count < max_iter)
    {   
        fichp << t << " ";
        fichv << t << " ";
        
        for(int i=0 ; i<plantCount ; i++)
            for(int site=0 ; site<numPatch ; site++)
                fichp << p[i][site] << " ";
        
        for(int i=0 ; i<insectCount ; i++)
            for(int site=0 ; site<numPatch ; site++)
                fichv << v[i][site] << " ";

        fichp << endl;
        fichv << endl;

        p_prev = p;
        v_prev = v;

        rungekutta(p,v,gamma,h,plantCount,insectCount,numPatch); 
        t+=h;   
        iter_count++;
        
        max_delta = 0.0;
        for(int i=0 ; i<plantCount ; i++)
        {
            for(int site=0 ; site<numPatch ; site++)
            {
                double delta = abs(p[i][site] - p_prev[i][site]);
                if(delta > max_delta)
                    max_delta = delta;
            }
        }
        
        for(int i=0 ; i<insectCount ; i++)
        {
            for(int site=0 ; site<numPatch ; site++)
            {
                double delta = abs(v[i][site] - v_prev[i][site]);
                if(delta > max_delta)
                    max_delta = delta;
            }
        }
        
        if (max_delta < TOLERANCE && iter_count > 1000) 
        {
            cout << "Stationary state at t = " << t << " (iter " << iter_count << ")" << endl;
            break;
        }
    }
    
    if (iter_count == max_iter)
        cout << "  No convergence." << endl;
    
    fichp << t << " ";
    fichv << t << " ";
    for(int i=0 ; i<plantCount ; i++)
        for(int site=0 ; site<numPatch ; site++)
            fichp << p[i][site] << " ";
    for(int i=0 ; i<insectCount ; i++)
        for(int site=0 ; site<numPatch ; site++)
            fichv << v[i][site] << " ";
    fichp << endl;
    fichv << endl;
    
    return;
}

void loadGamma(const string &filename, map<string, int>& plantIndex, map<string, int>& insectIndex, int& plantCount, int& insectCount, int& numPatch, vector<vector<vector<double>>>& gamma)
{
    ifstream intfich(filename);

    string plant, insect;
    double weight;
    int site;

    string line;
    
    vector<tuple<int, string, string, double>> data;
    
    while (getline(intfich, line)) 
    {
        istringstream iss(line);
        if (!(iss >> site >> plant >> insect >> weight)) continue;

        if (plantIndex.find(plant) == plantIndex.end())
            plantIndex[plant] = plantCount++;
        if (insectIndex.find(insect) == insectIndex.end())
            insectIndex[insect] = insectCount++;
        
        if (site + 1 > numPatch)
            numPatch = site + 1;

        data.push_back({site, plant, insect, weight});
    }
    
    intfich.close();
    
    cout << endl << numPatch << " patches." << endl << endl;
    
    gamma.resize(numPatch, vector<vector<double>>(plantCount, vector<double>(insectCount, 0.0)));
    
    for (const auto& item : data) 
    {
        int p_id = get<0>(item);
        string pl = get<1>(item);
        string ins = get<2>(item);
        double w = get<3>(item);
        
        int i = plantIndex[pl];
        int j = insectIndex[ins];
        gamma[p_id][i][j] = w;
    }
    
    return;
}

void runExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numPatch)
{
    cout << "\n---Initializing extinction experiment ---" << endl;
    
    //Plant ranking
    vector<pair<double, int>> plantRanking;
    
    for(int i=0 ; i<plantCount ; i++)
    {
        double abundance = 0.0;
        for(int site=0 ; site<numPatch ; site++)
        {
            abundance += p[i][site];
        }
        plantRanking.push_back({abundance,i});
    }
    
    sort(plantRanking.begin(), plantRanking.end());
    
    
    //Preparation
    ofstream experimentFile("results.txt");
    
    experimentFile << "# Num_Extinctions Robustness_Ratio Surv_Plants Surv_Insects Pollination_Service Gini_Plants Gini_Insects" << endl;
    
    ofstream dummy_p("/dev/null");
    ofstream dummy_v("/dev/null");
    
    vector<vector<double>> pCurrent = p;
    vector<vector<double>> vCurrent = v;
    
    double tDummy = 0.0;
    
    //Experiment
    for(int k=0; k<=plantCount ; k++)
    {
        if (k>0)
        {
            int plantToRemove = plantRanking[k-1].second;
            
            for(int site=0 ; site<numPatch ; site++)
                pCurrent[plantToRemove][site] = 0.0;
            
            findSteadyState(tDummy, pCurrent, vCurrent, dummy_p, dummy_v, plantCount, insectCount, numPatch, gamma, h);
        }
        
        //Metrics
        int survPlants = 0;
        int survInsects = 0;
        double plantBiomass = 0.0, insectBiomass = 0.0;
        double sum2p = 0.0, sum2v = 0.0;
        
        for (int i=0 ; i<plantCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += pCurrent[i][site];
            
            if(total > viability) 
            {
                survPlants++;
                plantBiomass += total;
            }            
        }
        
        if (plantBiomass > 0)
        {
            for(int i=0 ; i<plantCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += pCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / plantBiomass;
                    sum2p += (prob * prob);
                }
            }
        }
        
        for (int i=0 ; i<insectCount ; i++)
        {
            double total = 0.0;
            for(int site=0 ; site<numPatch ; site++)
                total += vCurrent[i][site];
            
            if(total > viability) 
            {
                survInsects++;
                insectBiomass += total;  
            }          
        }
        
        if (insectBiomass > 0)
        {
            for(int i=0 ; i<insectCount ; i++)
            {
                double total = 0.0;
                for(int site=0 ; site<numPatch ; site++)
                    total += vCurrent[i][site];
                if(total > viability)
                {
                    double prob = total / insectBiomass;
                    sum2v += (prob * prob);
                }
            }
        }
        
        double robustness = (double)(survPlants + survInsects) / (plantCount + insectCount);
        double pollinationService = insectBiomass;
        double giniP = 1.0 - sum2p;
        double giniV = 1.0 - sum2v;
        
        experimentFile << k << " " << fixed << setprecision(6) << robustness << " " << survPlants << " " << survInsects << " " << pollinationService << " " << giniP << " " << giniV << endl;
    }
    
    experimentFile.close();
    dummy_p.close();
    dummy_v.close();
    
    cout << "---Extinction experiment complete ---" << endl;
    
    return;
}
