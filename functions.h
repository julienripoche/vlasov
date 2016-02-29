#ifndef DEF_FUNCTIONS
#define DEF_FUNCTIONS

#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

//Phase coordinates generator
double alea();
vector<vector<double> > sphere(double radius_max, unsigned int nbr_points);
void write(string const nomFichier, vector<vector<double> > &data);

//Evolution
void read(string const nomFichier, vector<vector<double> > &data);
double module(vector<double> coords, int i, int j, int k, double l);
double gaussian(double modulus, double sigma);
vector<double> harm_osc_str(vector<double> coords);

#endif
