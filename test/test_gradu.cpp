#include <iostream>
#include <iomanip> //setprecision
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "functions.h"

using namespace std;

int main()
{
    //Initialize positions and momenta values
    vector<vector<double> > r;
    vector<vector<double> > p;
    read(r, "coords.gnu");
    read(p, "momenta.gnu");

    //Useful variables
    int n_box = 30;
    double L = 20.;
    double l0 = L/n_box;
    int nbr_sigma = 2;

    //Initialize strength
    vector<vector<double> > F(_NA_, vector<double>(3,0)); 
    vector<double> rho_map(n_box*n_box*n_box,0);

    //Get rho and gradu maps
    rho(rho_map, r, l0, n_box);
    minus_gradU(F[0], rho_map, r[0], nbr_sigma, n_box, l0);

    //Display results
    cout << F[0][0] << " " << F[0][1] << " " << F[0][2] << endl;

    return 0;
}

