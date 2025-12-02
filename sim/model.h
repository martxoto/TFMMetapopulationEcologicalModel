#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iomanip>
#include <random>

using namespace std;

constexpr double m = 0.2;
constexpr double r =  2.0;
constexpr int Kp = 100;
constexpr int Kv = 1000;
constexpr double ha = 1.0;
constexpr double d = 0.3;
//constexpr double D = 2.5;
constexpr double viability = 1e-6;
constexpr double alpha = 1.0;



bool plantExistsInPatch(int plantID, int site, int insectCount, const vector<vector<vector<double>>>& gamma);
bool insectExistsInPatch(int insectID, int site, int plantCount, const vector<vector<vector<double>>>& gamma);

void evaluaFp(double p, const vector<vector<double>>& v, double &fp, const vector<vector<vector<double>>>& gamma, int pindex, int site, int insectCount);

void evaluaFv(const vector<vector<double>>& p, const vector<vector<double>>& v, double &fv, const vector<vector<vector<double>>>& gamma, int vindex, int site, int plantCount, int numpatch, double D);

void rungekutta(vector<vector<double>>& p, vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);

void findSteadyState(double t, vector<vector<double>>& p, vector<vector<double>>& v, ofstream& fichp, ofstream& fichv, int plantCount, int insectCount, int numpatch,  const vector<vector<vector<double>>>& gamma, double h, double D);

void loadGamma(const string &filename, map<string, int>& plantIndex, map<string, int>& insectIndex, int& plantCount, int& insectCount, int& numpatch, vector<vector<vector<double>>>& gamma);

void runExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);


void runRandomExtinctionExperiment(const vector<vector<double>>& p, const vector<vector<double>>& v, const vector<vector<vector<double>>>& gamma, double h, int plantCount, int insectCount, int numpatch, double D);

#endif
