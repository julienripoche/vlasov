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
double module(vector<double> &r_real, vector<double> &r_box);
double gaussian(vector<double> &r_real, vector<double> &r_box, double sigma);
int key(int x, int y, int z, int N);
void rho_map(vector<double> &rho, vector<vector<double> > &coords, double sigma, double L, int box_size, int N);
double U(double rho);
vector<double> minus_gradU(vector<double> &r, vector<double> &rho_map, double nbr_sigma, int N, double L, double sigma);
vector<double> minus_gradU2(vector<double> &r, vector<double> &rho_map, double nbr_sigma, int N, double L, double sigma);

#endif
