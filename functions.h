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

//Phase coordinates generator functions
double alea();
void sphere(vector<vector<double> > &r, double radius_max);

//Write and read functions
void write(vector<vector<double> > &data, string const nomFichier);
void read(vector<vector<double> > &data, string const nomFichier);

//Useful functions
int key(int x, int y, int z, int N);
double module(vector<double> &r_real, vector<double> &r_box);
double gaussian(vector<double> &r_real, vector<double> &r_box, double sigma);
double U(double rho);

//Evolution functions
void rho_map(vector<double> &rho, vector<vector<double> > &coords, double sigma, double L, int box_size, int N);
void minus_gradU(vector<double> &gradu, vector<double> &r, vector<double> &rho_map, double nbr_sigma, int N, double L, double sigma);

#endif
