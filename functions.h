#ifndef DEF_FUNCTIONS
#define DEF_FUNCTIONS

#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>

#define _A_ 56
#define _N_ 100
#define _NA_ (_N_ * _A_)
#define _SIGMA_ (double) 0.5
#define _SIGMA_NBR_ 3
#define _BOX_NBR_ 40
#define _L_ 20
#define _BOX_NBR_X_ 40
#define _BOX_NBR_Y_ 40
#define _BOX_NBR_Z_ 80
#define _L_X_ 20
#define _L_Y_ 20
#define _L_Z_ 40
#define _L0_ (double) _L_ / _BOX_NBR_
#define _R0_ 1.12
#define _HBAR_C_ 197.3
#define _M_ ( 938.3 + 939.6) / 2
#define _DT_ 0.5
#define _A_WS_ (double) _SIGMA_

using namespace std;

//Phase coordinates generator functions
double alea();
void sphere(vector<vector<double> > &r, double radius_max);
void coords_generate(vector<vector<double> > &r, double radius_max);
void momenta_generate(vector<vector<double> > &r, vector<vector<double> > &p);
double rho_ws(double r);
double fermi_momentum(double rho);

//Write and read functions
void write(vector<vector<double> > &data, string const nomFichier);
void read(vector<vector<double> > &data, string const nomFichier);

//Useful functions
int key(int x, int y, int z, int N);
double module(vector<double> &r);
double vect_module(vector<double> &r_real, vector<double> &r_box);
double gaussian(vector<double> &r_real, vector<double> &r_box);
double U(double rho);

//Evolution functions
void rho(vector<double> &rho_map, vector<vector<double> > &coords);
void minus_gradU(vector<double> &gradu, vector<double> &rho_map, vector<double> &r);
double get_ubar(vector<double> &rho_map, vector<double> &r);

#endif
