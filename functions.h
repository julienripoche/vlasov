#ifndef DEF_FUNCTIONS
#define DEF_FUNCTIONS

#define _A_ 56
#define _N_ 20
#define _NA_ (_N_ * _A_)
#define _SIGMA_ (double) 1
#define _BOX_NBR_ 30
#define _L_ 20
#define _L0_ (double) _L_ / _BOX_NBR_

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
void rho(vector<double> &rho_map, vector<vector<double> > &coords, double l0, int box_nbr);
void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r, double l0, int box_nbr, int sigma_nbr);
double get_ubar(vector<double> &rho_map, vector<double> &r, double l0, int box_nbr, int sigma_nbr);

void coords_generate(vector<vector<double> > &r, double radius_max);
void momenta_generate(vector<vector<double> > &r, vector<vector<double> > &p);

#endif
