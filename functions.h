#ifndef DEF_FUNCTIONS
#define DEF_FUNCTIONS

#define _A_ 56
#define _N_ 5
#define _NA_ (_N_ * _A_)
#define _SIGMA_ 1

#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

//Phase coordinates generator functions
double alea();
void sphere(vector<vector<double> > &r, double radius_max);

//Write and read functions
void write(vector<vector<double> > &data, string const nomFichier);
void read(vector<vector<double> > &data, string const nomFichier);

//Useful functions
int key(int x, int y, int z, int N);
double module(vector<double> &r_real, vector<double> &r_box);
double gaussian(vector<double> &r_real, vector<double> &r_box);
double U(double rho);

//Evolution functions
void rho(vector<double> &rho_map, vector<vector<double> > &coords, double l0, int box_size);
void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r, int nbr_sigma, int N, double l0);

#endif
